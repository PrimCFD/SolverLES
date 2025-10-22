// File: PressurePoisson.cpp  (matrix-free Amat + Galerkin MG (PMAT), PETSc 3.20+)
//
// This version keeps a matrix-free MatShell as the fine-grid operator (Amat) and
// assembles a fine-grid AIJ as Pmat; PCMG forms coarse PMATs variationally via Galerkin.
// Transfers use DMDA_Q0 (cell-centered) interpolation. Level smoothers are Jacobi-Chebyshev
// Outer solve is PiPeCG (no communication CG) or MINRES if poissson is singular (all Neumann).
//
// Boundary treatment matches the assembled path:
//   - Neumann: natural BC (zero flux unless s->value specified -> RHS shift)
//   - Dirichlet: face coeff added to diagonal; RHS shifted by coeff*value
//
// References (PETSc 3.20/3.23 docs):
// - PCMG & Galerkin: PCMGSetGalerkin(), PCMGSetOperators()
// - DMDA interpolation: DMDASetInterpolationType(DMDA_Q0) and DMCreateInterpolation()
// - KSP/Amat/Pmat: KSPSetOperators(); Chebyshev smoother: KSPCHEBYSHEV/KSPChebyshevEstEigSet()

#include "PressurePoisson.hpp"
#include "Program.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include "master/Views.hpp"
#include "kernels_fluids.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// PETSc
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>

#include "memory/MpiBox.hpp"

#ifdef HAVE_MPI
static inline MPI_Comm valid_comm_or(MPI_Comm cand, MPI_Comm fallback)
{
    if (cand == MPI_COMM_NULL)
        return fallback;
    int sz = -1, rc = MPI_Comm_size(cand, &sz);
    if (rc == MPI_SUCCESS && sz > 0)
        return cand;
    return fallback;
}
#endif

using core::master::FieldCatalog;
using core::master::MeshTileView;

namespace fluids
{

// --------------------------- small helpers ---------------------------------

static inline double hmean(double a, double b)
{
    const double eps = 1.0e-300;
    a = std::max(a, eps);
    b = std::max(b, eps);
    return 2.0 / (1.0 / a + 1.0 / b);
}

// pressure BC bundle
struct PBC
{
    const BcSpec* W;
    const BcSpec* E;
    const BcSpec* S;
    const BcSpec* N;
    const BcSpec* B;
    const BcSpec* T;
};

static inline BcSpec::Type norm_p_type(BcSpec::Type t)
{
    // For pressure, "mirror" behaves like homogeneous Neumann
    if (t == BcSpec::Type::mirror)
        return BcSpec::Type::neumann;
    return t;
}

// --------------------------- MatShell context -------------------------------

struct LevelBeta
{
    // β = 1/ρ at cell-centers for THIS MG level, interior-only (no ghosts).
    // size = nxi*nyi*nzi for the level.
    std::vector<double> data;
    int nxi = 0, nyi = 0, nzi = 0;
    inline double B(int i, int j, int k) const noexcept
    {
        return data[std::size_t(i) +
                    std::size_t(nxi) * (std::size_t(j) + std::size_t(nyi) * std::size_t(k))];
    }
};

// Per-level face transmissibilities (TPFA conductances).
// Stored as T = β * (A / Δ), i.e. T_x = β * (dy*dz / dx), T_y = β * (dx*dz / dy), T_z = β * (dx*dy
// / dz). The discrete operator uses (1/Vol) * sum_faces( T_face * (p_nb - p_c) ).
struct LevelTrans
{
    int nxi = 0, nyi = 0, nzi = 0; // cell counts
    // Faces: x-faces sit between (i,j,k) and (i+1,j,k), i in [0, nxi-2]
    //        y-faces sit between (i,j,k) and (i,j+1,k), j in [0, nyi-2]
    //        z-faces sit between (i,j,k) and (i,j,k+1), k in [0, nzi-2]
    std::vector<double> Tx, Ty, Tz;
    inline std::size_t idxTx(int i, int j, int k) const
    {
        return std::size_t(i) +
               std::size_t(nxi - 1) * (std::size_t(j) + std::size_t(nyi) * std::size_t(k));
    }
    inline std::size_t idxTy(int i, int j, int k) const
    {
        return std::size_t(j) +
               std::size_t(nyi - 1) * (std::size_t(i) + std::size_t(nxi) * std::size_t(k));
    }
    inline std::size_t idxTz(int i, int j, int k) const
    {
        return std::size_t(k) +
               std::size_t(nzi - 1) * (std::size_t(i) + std::size_t(nxi) * std::size_t(j));
    }
    inline double TX(int i, int j, int k) const { return Tx[idxTx(i, j, k)]; }
    inline double TY(int i, int j, int k) const { return Ty[idxTy(i, j, k)]; }
    inline double TZ(int i, int j, int k) const { return Tz[idxTz(i, j, k)]; }
};

struct ShellCtx
{
    DM da;                         // level DMDA (interior)
    int nxi, nyi, nzi;             // interior extents (no ghosts)
    int ng;                        // original ghost width in user arrays (for indexing RHS shifts)
    int nxc_tot, nyc_tot, nzc_tot; // original center total sizes (fine grid)
    double dx, dy, dz;             // level spacings
    const LevelTrans* trans = nullptr; // face transmissibilities on THIS level (β*A/Δ)
    const LevelBeta* beta = nullptr;   // cell-centered β on THIS level (interior only)
    // For constant ρ path when trans is null:
    bool use_const_beta = false;
    double beta_const = 1.0;

    PBC pbc;

    // cached diagonals to speed repeated GET_DIAGONAL calls (rebuilt when spacing/beta changes)
    Vec diag = nullptr;
    bool diag_built = false;

    // periodicity per axis (derived from p.* BCs)
    bool perX = false, perY = false, perZ = false;
};

// MULT: y = A x  (matrix-free  (1/V) ∑_faces T (x_nb - x_c)  )
static PetscErrorCode ShellMult(Mat A, Vec x, Vec y)
{
    PetscFunctionBegin;
    ShellCtx* ctx = nullptr;
    PetscCall(MatShellGetContext(A, (void**) &ctx));
    const bool useTx = (ctx->trans && !ctx->trans->Tx.empty());
    const bool useTy = (ctx->trans && !ctx->trans->Ty.empty());
    const bool useTz = (ctx->trans && !ctx->trans->Tz.empty());

    const double invV = 1.0 / (ctx->dx * ctx->dy * ctx->dz);

    // Constant-β fallback face conductances (used whenever no per-face Tx/Ty/Tz)
    const double TXc = ctx->beta_const * (ctx->dy * ctx->dz / ctx->dx);
    const double TYc = ctx->beta_const * (ctx->dx * ctx->dz / ctx->dy);
    const double TZc = ctx->beta_const * (ctx->dx * ctx->dy / ctx->dz);
    const bool constBeta = ctx->use_const_beta;

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscCall(DMDAGetCorners(ctx->da, &xs, &ys, &zs, &xm, &ym, &zm));

    Vec xloc;
    PetscCall(DMGetLocalVector(ctx->da, &xloc));
    PetscCall(DMGlobalToLocalBegin(ctx->da, x, INSERT_VALUES, xloc));
    PetscCall(DMGlobalToLocalEnd(ctx->da, x, INSERT_VALUES, xloc));

    const PetscScalar*** xarr;
    PetscScalar*** yarr;
    PetscCall(DMDAVecGetArrayRead(ctx->da, xloc, &xarr));
    PetscCall(DMDAVecGetArray(ctx->da, y, &yarr));

#pragma omp parallel for collapse(3) schedule(static)
    for (PetscInt k = zs; k < zs + zm; ++k)
        for (PetscInt j = ys; j < ys + ym; ++j)
            for (PetscInt i = xs; i < xs + xm; ++i)
            {
                const double xc = xarr[k][j][i];

                // neighbor presence (periodic-aware)
                const bool hasW = (i > 0) || ctx->perX;
                const bool hasE = (i < ctx->nxi - 1) || ctx->perX;
                const bool hasS = (j > 0) || ctx->perY;
                const bool hasN = (j < ctx->nyi - 1) || ctx->perY;
                const bool hasB = (k > 0) || ctx->perZ;
                const bool hasT = (k < ctx->nzi - 1) || ctx->perZ;

                // local neighbor indices (ghosts hold periodic wraps)
                const PetscInt iw = i - 1;
                const PetscInt ie = i + 1;
                const PetscInt js = j - 1;
                const PetscInt jn = j + 1;
                const PetscInt kb = k - 1;
                const PetscInt kt = k + 1;

                // Face transmissibilities (conductances) computed once per cell.
                // For interior faces, use precomputed Tx/Ty/Tz if present, else constant-β TXc/TYc/TZc.
                // For boundary faces, fall back to local β * (A/Δ).
                const double bC = constBeta ? ctx->beta_const : ctx->beta->B(i, j, k);

                double Tw = 0.0, Te = 0.0, Ts = 0.0, Tn = 0.0, Tb = 0.0, Tt = 0.0;
                // X faces
                if (i > 0)        Tw = useTx ? ctx->trans->TX(i - 1, j, k) : TXc;
                else              Tw = bC * (ctx->dy * ctx->dz / ctx->dx);
                if (i < ctx->nxi - 1) Te = useTx ? ctx->trans->TX(i, j, k)     : TXc;
                else                   Te = bC * (ctx->dy * ctx->dz / ctx->dx);
                // Y faces
                if (j > 0)        Ts = useTy ? ctx->trans->TY(i, j - 1, k) : TYc;
                else              Ts = bC * (ctx->dx * ctx->dz / ctx->dy);
                if (j < ctx->nyi - 1) Tn = useTy ? ctx->trans->TY(i, j, k)     : TYc;
                else                   Tn = bC * (ctx->dx * ctx->dz / ctx->dy);
                // Z faces
                if (k > 0)        Tb = useTz ? ctx->trans->TZ(i, j, k - 1) : TZc;
                else              Tb = bC * (ctx->dx * ctx->dy / ctx->dz);
                if (k < ctx->nzi - 1) Tt = useTz ? ctx->trans->TZ(i, j, k)     : TZc;
                else                   Tt = bC * (ctx->dx * ctx->dy / ctx->dz);

                double acc = 0.0;
                if (hasW)
                    acc += Tw * (xarr[k][j][iw] - xc);
                else if (!ctx->perX && ctx->pbc.W &&
                         norm_p_type(ctx->pbc.W->type) == BcSpec::Type::dirichlet)
                    acc += Tw * ((ctx->pbc.W ? ctx->pbc.W->value : 0.0) - xc);

                if (hasE)
                    acc += Te * (xarr[k][j][ie] - xc);
                else if (!ctx->perX && ctx->pbc.E &&
                         norm_p_type(ctx->pbc.E->type) == BcSpec::Type::dirichlet)
                    acc += Te * ((ctx->pbc.E ? ctx->pbc.E->value : 0.0) - xc);

                if (hasS)
                    acc += Ts * (xarr[k][js][i] - xc);
                else if (!ctx->perY && ctx->pbc.S &&
                         norm_p_type(ctx->pbc.S->type) == BcSpec::Type::dirichlet)
                    acc += Ts * ((ctx->pbc.S ? ctx->pbc.S->value : 0.0) - xc);

                if (hasN)
                    acc += Tn * (xarr[k][jn][i] - xc);
                else if (!ctx->perY && ctx->pbc.N &&
                         norm_p_type(ctx->pbc.N->type) == BcSpec::Type::dirichlet)
                    acc += Tn * ((ctx->pbc.N ? ctx->pbc.N->value : 0.0) - xc);

                if (hasB)
                    acc += Tb * (xarr[kb][j][i] - xc);
                else if (!ctx->perZ && ctx->pbc.B &&
                         norm_p_type(ctx->pbc.B->type) == BcSpec::Type::dirichlet)
                    acc += Tb * ((ctx->pbc.B ? ctx->pbc.B->value : 0.0) - xc);

                if (hasT)
                    acc += Tt * (xarr[kt][j][i] - xc);
                else if (!ctx->perZ && ctx->pbc.T &&
                         norm_p_type(ctx->pbc.T->type) == BcSpec::Type::dirichlet)
                    acc += Tt * ((ctx->pbc.T ? ctx->pbc.T->value : 0.0) - xc);

                yarr[k][j][i] = -invV * acc;
            }

    PetscCall(DMDAVecRestoreArrayRead(ctx->da, xloc, &xarr));
    PetscCall(DMDAVecRestoreArray(ctx->da, y, &yarr));
    PetscCall(DMRestoreLocalVector(ctx->da, &xloc));
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode ShellBuildDiagonal(ShellCtx* ctx)
{
    const bool useTx = (ctx->trans && !ctx->trans->Tx.empty());
    const bool useTy = (ctx->trans && !ctx->trans->Ty.empty());
    const bool useTz = (ctx->trans && !ctx->trans->Tz.empty());
    PetscFunctionBegin;

    // Precompute constants for constant-β path
    const double TXc = ctx->beta_const * (ctx->dy * ctx->dz / ctx->dx);
    const double TYc = ctx->beta_const * (ctx->dx * ctx->dz / ctx->dy);
    const double TZc = ctx->beta_const * (ctx->dx * ctx->dy / ctx->dz);
    const bool constBeta = ctx->use_const_beta;

    if (!ctx->diag)
    {
        PetscCall(DMCreateGlobalVector(ctx->da, &ctx->diag));
    }

    const double invV = 1.0 / (ctx->dx * ctx->dy * ctx->dz);

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscCall(DMDAGetCorners(ctx->da, &xs, &ys, &zs, &xm, &ym, &zm));

    PetscScalar*** darr;
    PetscCall(DMDAVecGetArray(ctx->da, ctx->diag, &darr));

    for (PetscInt k = zs; k < zs + zm; ++k)
        for (PetscInt j = ys; j < ys + ym; ++j)
            for (PetscInt i = xs; i < xs + xm; ++i)
            {
                const bool hasW = (i > 0) || ctx->perX;
                const bool hasE = (i < ctx->nxi - 1) || ctx->perX;
                const bool hasS = (j > 0) || ctx->perY;
                const bool hasN = (j < ctx->nyi - 1) || ctx->perY;
                const bool hasB = (k > 0) || ctx->perZ;
                const bool hasT = (k < ctx->nzi - 1) || ctx->perZ;
                const double bC = constBeta ? ctx->beta_const : ctx->beta->B(i, j, k);
                double Tw, Te, Ts, Tn, Tb, Tt;
                // X faces
                if (i > 0)             Tw = useTx ? ctx->trans->TX(i - 1, j, k) : TXc;
                else                   Tw = bC * (ctx->dy * ctx->dz / ctx->dx);
                if (i < ctx->nxi - 1)  Te = useTx ? ctx->trans->TX(i, j, k)     : TXc;
                else                   Te = bC * (ctx->dy * ctx->dz / ctx->dx);
                // Y faces
                if (j > 0)             Ts = useTy ? ctx->trans->TY(i, j - 1, k) : TYc;
                else                   Ts = bC * (ctx->dx * ctx->dz / ctx->dy);
                if (j < ctx->nyi - 1)  Tn = useTy ? ctx->trans->TY(i, j, k)     : TYc;
                else                   Tn = bC * (ctx->dx * ctx->dz / ctx->dy);
                // Z faces
                if (k > 0)             Tb = useTz ? ctx->trans->TZ(i, j, k - 1) : TZc;
                else                   Tb = bC * (ctx->dx * ctx->dy / ctx->dz);
                if (k < ctx->nzi - 1)  Tt = useTz ? ctx->trans->TZ(i, j, k)     : TZc;
                else                   Tt = bC * (ctx->dx * ctx->dy / ctx->dz);

                double diag = 0.0;
                if (hasE) diag += Te;
                if (hasW) diag += Tw;
                if (hasN) diag += Tn;
                if (hasS) diag += Ts;
                if (hasT) diag += Tt;
                if (hasB) diag += Tb;

                // Dirichlet faces contribute to diagonal
                if (!hasW && !ctx->perX && ctx->pbc.W &&
                    norm_p_type(ctx->pbc.W->type) == BcSpec::Type::dirichlet)
                    diag += Tw;
                if (!hasE && !ctx->perX && ctx->pbc.E &&
                    norm_p_type(ctx->pbc.E->type) == BcSpec::Type::dirichlet)
                    diag += Te;
                if (!hasS && !ctx->perY && ctx->pbc.S &&
                    norm_p_type(ctx->pbc.S->type) == BcSpec::Type::dirichlet)
                    diag += Ts;
                if (!hasN && !ctx->perY && ctx->pbc.N &&
                    norm_p_type(ctx->pbc.N->type) == BcSpec::Type::dirichlet)
                    diag += Tn;
                if (!hasB && !ctx->perZ && ctx->pbc.B &&
                    norm_p_type(ctx->pbc.B->type) == BcSpec::Type::dirichlet)
                    diag += Tb;
                if (!hasT && !ctx->perZ && ctx->pbc.T &&
                    norm_p_type(ctx->pbc.T->type) == BcSpec::Type::dirichlet)
                    diag += Tt;

                darr[k][j][i] = invV * diag;
            }

    PetscCall(DMDAVecRestoreArray(ctx->da, ctx->diag, &darr));
    ctx->diag_built = true;
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode ShellGetDiagonal(Mat A, Vec d)
{
    PetscFunctionBegin;
    ShellCtx* ctx = nullptr;
    PetscCall(MatShellGetContext(A, (void**) &ctx));
    if (!ctx->diag_built)
        PetscCall(ShellBuildDiagonal(ctx));
    PetscCall(VecCopy(ctx->diag, d));
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode AssembleAIJFromShell(const ShellCtx& ctx, Mat* Aout)
{
    const bool useTx = (ctx.trans && !ctx.trans->Tx.empty());
    const bool useTy = (ctx.trans && !ctx.trans->Ty.empty());
    const bool useTz = (ctx.trans && !ctx.trans->Tz.empty());
    PetscFunctionBegin;
    DM da = ctx.da;

    AO ao = NULL;
    PetscCall(DMDAGetAO(da, &ao));
    const PetscInt nxi = ctx.nxi, nyi = ctx.nyi, nzi = ctx.nzi;

    // map (i,j,k) -> global row/col using AO
    auto gid = [&](PetscInt i, PetscInt j, PetscInt k) -> PetscInt
    {
        PetscInt app = i + nxi * (j + nyi * k); // application/global logical index
        AOApplicationToPetsc(ao, 1, &app);      // map to PETSc global id for this DMDA layout
        return app;
    };

    // Create a properly preallocated AIJ using the DMDA stencil (STAR, width=1)
    PetscCall(DMSetMatType(da, MATAIJ));
    PetscCall(DMCreateMatrix(da, Aout)); // Aout gets AIJ with 7-pt prealloc
    PetscCall(MatSetOption(*Aout, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

    const double invV = 1.0 / (ctx.dx * ctx.dy * ctx.dz);

    PetscInt xs, ys, zs, xm, ym, zm;
    PetscCall(DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm));

    auto Tw = [&](int i, int j, int k)
    {
        if (i > 0)
            return useTx ? ctx.trans->TX(i - 1, j, k) : ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dy * ctx.dz / ctx.dx);
    };
    auto Te = [&](int i, int j, int k)
    {
        if (i < ctx.nxi - 1)
            return useTx ? ctx.trans->TX(i, j, k) : ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dy * ctx.dz / ctx.dx);
    };
    auto Ts = [&](int i, int j, int k)
    {
        if (j > 0)
            return useTy ? ctx.trans->TY(i, j - 1, k) : ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dx * ctx.dz / ctx.dy);
    };
    auto Tn = [&](int i, int j, int k)
    {
        if (j < ctx.nyi - 1)
            return useTy ? ctx.trans->TY(i, j, k) : ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dx * ctx.dz / ctx.dy);
    };
    auto Tb = [&](int i, int j, int k)
    {
        if (k > 0)
            return useTz ? ctx.trans->TZ(i, j, k - 1) : ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dx * ctx.dy / ctx.dz);
    };
    auto Tt = [&](int i, int j, int k)
    {
        if (k < ctx.nzi - 1)
            return useTz ? ctx.trans->TZ(i, j, k) : ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz);
        const double b = ctx.use_const_beta ? ctx.beta_const : ctx.beta->B(i, j, k);
        return b * (ctx.dx * ctx.dy / ctx.dz);
    };

    for (PetscInt k = zs; k < zs + zm; ++k)
        for (PetscInt j = ys; j < ys + ym; ++j)
            for (PetscInt i = xs; i < xs + xm; ++i)
            {
                // coefficients (same logic as ShellMult / ShellGetDiagonal)
                const bool hasW = (i > 0) || ctx.perX;
                const bool hasE = (i < ctx.nxi - 1) || ctx.perX;
                const bool hasS = (j > 0) || ctx.perY;
                const bool hasN = (j < ctx.nyi - 1) || ctx.perY;
                const bool hasB = (k > 0) || ctx.perZ;
                const bool hasT = (k < ctx.nzi - 1) || ctx.perZ;

                const PetscInt iw = (i > 0) ? (i - 1) : (ctx.nxi - 1);
                const PetscInt ie = (i < ctx.nxi - 1) ? (i + 1) : 0;
                const PetscInt js = (j > 0) ? (j - 1) : (ctx.nyi - 1);
                const PetscInt jn = (j < ctx.nyi - 1) ? (j + 1) : 0;
                const PetscInt kb = (k > 0) ? (k - 1) : (ctx.nzi - 1);
                const PetscInt kt = (k < ctx.nzi - 1) ? (k + 1) : 0;

                const double te = hasE ? Te(i, j, k) : 0.0;
                const double tw = hasW ? Tw(i, j, k) : 0.0;
                const double tn = hasN ? Tn(i, j, k) : 0.0;
                const double ts = hasS ? Ts(i, j, k) : 0.0;
                const double tt = hasT ? Tt(i, j, k) : 0.0;
                const double tb = hasB ? Tb(i, j, k) : 0.0;

                // central diagonal and 6 neighbors (Dirichlet faces add to diagonal;
                // Neumann faces affect only RHS -> not part of the matrix)
                double diag = 0.0;

                PetscInt row = gid(i, j, k);
                PetscInt cols[7];
                PetscScalar vals[7];
                int nz = 0;

                auto add_off = [&](int ii, int jj, int kk, double coef)
                {
                    cols[nz] = gid(ii, jj, kk);
                    vals[nz] = -invV * coef;
                    nz++;
                    diag += coef;
                };

                if (hasW)
                    add_off(iw, j, k, tw);
                else if (!ctx.perX && ctx.pbc.W &&
                         norm_p_type(ctx.pbc.W->type) == BcSpec::Type::dirichlet)
                    diag += Tw(i, j, k);

                if (hasE)
                    add_off(ie, j, k, te);
                else if (!ctx.perX && ctx.pbc.E &&
                         norm_p_type(ctx.pbc.E->type) == BcSpec::Type::dirichlet)
                    diag += Te(i, j, k);

                if (hasS)
                    add_off(i, js, k, ts);
                else if (!ctx.perY && ctx.pbc.S &&
                         norm_p_type(ctx.pbc.S->type) == BcSpec::Type::dirichlet)
                    diag += Ts(i, j, k);

                if (hasN)
                    add_off(i, jn, k, tn);
                else if (!ctx.perY && ctx.pbc.N &&
                         norm_p_type(ctx.pbc.N->type) == BcSpec::Type::dirichlet)
                    diag += Tn(i, j, k);

                if (hasB)
                    add_off(i, j, kb, tb);
                else if (!ctx.perZ && ctx.pbc.B &&
                         norm_p_type(ctx.pbc.B->type) == BcSpec::Type::dirichlet)
                    diag += Tb(i, j, k);

                if (hasT)
                    add_off(i, j, kt, tt);
                else if (!ctx.perZ && ctx.pbc.T &&
                         norm_p_type(ctx.pbc.T->type) == BcSpec::Type::dirichlet)
                    diag += Tt(i, j, k);

                // center
                cols[nz] = row;
                vals[nz] = invV * diag;
                nz++;

                PetscCall(MatSetValues(*Aout, 1, &row, nz, cols, vals, INSERT_VALUES));
            }

    PetscCall(MatAssemblyBegin(*Aout, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*Aout, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
}

// --------------------------- Impl & builder ---------------------------------

struct PPImpl
{
    // geometry / sizes
    int nxc_tot{}, nyc_tot{}, nzc_tot{}, ng{};
    double dx{}, dy{}, dz{};
    bool varrho{};
    double rho_const{};
    bool all_neumann{};

    // Borrowed communicator owned by the main application.
    // Lifetime is managed externally; we must not duplicate or free it.
    MPI_Comm comm{MPI_COMM_NULL};

    // PETSc objects
    DM da_fine{};                      // finest DM (interior)
    std::vector<DM> da_lvl;            // 0=coarsest ... L-1=finest
    std::vector<Mat> A_lvl;            // level operators (MatShell)
    std::vector<LevelBeta> beta_lvl;   // per-level β (only used to build Tx/Ty/Tz each step)
    std::vector<LevelTrans> trans_lvl; // per-level face transmissibilities
    std::vector<std::unique_ptr<ShellCtx>> ctx_lvl;
    KSP ksp{};       // outer CG (MG preconditioner inside)
    Vec x{}, b{};    // fine-level unknown/RHS
    Mat Afine_aij{}; // assembled AIJ on finest level (PMAT for MG + outer KSP)
    // Reusable nullspace for pure Neumann problems (created once per hierarchy)
    MatNullSpace ns_const{};

    // host scratch for p
    std::vector<double> p_host;

    // full teardown (safe on partial initialization)
    void destroy()
    {

        if (Afine_aij)
        {
            MatDestroy(&Afine_aij);
            Afine_aij = nullptr;
        }

        // vectors
        if (x)
        {
            VecDestroy(&x);
            x = nullptr;
        }
        if (b)
        {
            VecDestroy(&b);
            b = nullptr;
        }
        // level mats (MatShell) – PCMG retains borrowed refs; destroy our refs
        for (auto& A : A_lvl)
            if (A)
            {
                MatDestroy(&A);
            }
        A_lvl.clear();
        // cached diagonals on each ctx
        for (auto& c : ctx_lvl)
            if (c && c->diag)
            {
                VecDestroy(&c->diag);
                c->diag = nullptr;
            }
        ctx_lvl.clear();
        // DMs (levels then fine)
        for (auto& d : da_lvl)
            if (d)
            {
                DMDestroy(&d);
            }
        da_lvl.clear();
        if (da_fine)
        {
            DMDestroy(&da_fine);
            da_fine = nullptr;
        }
        // KSP
        if (ksp)
        {
            KSPDestroy(&ksp);
            ksp = nullptr;
        }
        if (ns_const)
        {
            MatNullSpaceDestroy(&ns_const);
            ns_const = nullptr;
        }
        // Borrowed communicator: do NOT free; just clear the handle.
        comm = MPI_COMM_NULL;
    }
};

static void average_beta_coarsen(const LevelBeta& fine, LevelBeta& coarse)
{
    // simple 2x restriction by arithmetic averaging of 2x2x2 children
    for (int K = 0; K < coarse.nzi; ++K)
        for (int J = 0; J < coarse.nyi; ++J)
            for (int I = 0; I < coarse.nxi; ++I)
            {
                const int i0 = 2 * I, j0 = 2 * J, k0 = 2 * K;
                double s = 0.0;
                int cnt = 0;
                for (int dk = 0; dk < 2; ++dk)
                    for (int dj = 0; dj < 2; ++dj)
                        for (int di = 0; di < 2; ++di)
                        {
                            const int ii = i0 + di, jj = j0 + dj, kk = k0 + dk;
                            if (ii < fine.nxi && jj < fine.nyi && kk < fine.nzi)
                            {
                                s += fine.B(ii, jj, kk);
                                ++cnt;
                            }
                        }
                coarse.data[std::size_t(I) +
                            std::size_t(coarse.nxi) *
                                (std::size_t(J) + std::size_t(coarse.nyi) * std::size_t(K))] =
                    s / std::max(1, cnt);
            }
}

// Build fine-level face transmissibilities from β and spacings
static void build_fine_trans_from_beta(const LevelBeta& B, double dx, double dy, double dz,
                                       LevelTrans& T)
{
    T.nxi = B.nxi;
    T.nyi = B.nyi;
    T.nzi = B.nzi;
    T.Tx.assign(std::size_t(std::max(0, B.nxi - 1)) * B.nyi * B.nzi, 0.0);
    T.Ty.assign(std::size_t(B.nxi) * std::max(0, B.nyi - 1) * B.nzi, 0.0);
    T.Tz.assign(std::size_t(B.nxi) * B.nyi * std::max(0, B.nzi - 1), 0.0);

    for (int k = 0; k < B.nzi; ++k)
        for (int j = 0; j < B.nyi; ++j)
            for (int i = 0; i < B.nxi - 1; ++i)
            {
                const double hf = hmean(B.B(i, j, k), B.B(i + 1, j, k));
                T.Tx[T.idxTx(i, j, k)] = hf * (dy * dz / dx);
            }
    for (int k = 0; k < B.nzi; ++k)
        for (int j = 0; j < B.nyi - 1; ++j)
            for (int i = 0; i < B.nxi; ++i)
            {
                const double hf = hmean(B.B(i, j, k), B.B(i, j + 1, k));
                T.Ty[T.idxTy(i, j, k)] = hf * (dx * dz / dy);
            }
    for (int k = 0; k < B.nzi - 1; ++k)
        for (int j = 0; j < B.nyi; ++j)
            for (int i = 0; i < B.nxi; ++i)
            {
                const double hf = hmean(B.B(i, j, k), B.B(i, j, k + 1));
                T.Tz[T.idxTz(i, j, k)] = hf * (dx * dy / dz);
            }
}

// Conservative 2x face-aggregation: coarse-face transmissibility is the SUM of
// fine-face transmissibilities across that coarse face (parallel conductances add).
static void coarsen_trans_2x(const LevelTrans& Tf, LevelTrans& Tc)
{
    // Geometry: coarse cell (I,J,K) covers fine [2I:2I+1]×[2J:2J+1]×[2K:2K+1]
    // X-faces at (I+1/2,J,K) collect fine x-faces at i=2I+1 and j,k in those 2x2 slabs.
    Tc.Tx.assign(std::size_t(std::max(0, Tc.nxi - 1)) * Tc.nyi * Tc.nzi, 0.0);
    Tc.Ty.assign(std::size_t(Tc.nxi) * std::max(0, Tc.nyi - 1) * Tc.nzi, 0.0);
    Tc.Tz.assign(std::size_t(Tc.nxi) * Tc.nyi * std::max(0, Tc.nzi - 1), 0.0);

    // X faces
    for (int K = 0; K < Tc.nzi; ++K)
        for (int J = 0; J < Tc.nyi; ++J)
            for (int I = 0; I < Tc.nxi - 1; ++I)
            {
                const int i = 2 * I + 1;
                double sum = 0.0;
                for (int dj = 0; dj < 2; ++dj)
                    for (int dk = 0; dk < 2; ++dk)
                    {
                        const int j = 2 * J + dj;
                        const int k = 2 * K + dk;
                        if (i >= 0 && i < Tf.nxi - 1 && j < Tf.nyi && k < Tf.nzi)
                            sum += Tf.TX(i, j, k);
                    }
                Tc.Tx[Tc.idxTx(I, J, K)] = sum;
            }
    // Y faces
    for (int K = 0; K < Tc.nzi; ++K)
        for (int J = 0; J < Tc.nyi - 1; ++J)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int j = 2 * J + 1;
                double sum = 0.0;
                for (int di = 0; di < 2; ++di)
                    for (int dk = 0; dk < 2; ++dk)
                    {
                        const int i = 2 * I + di;
                        const int k = 2 * K + dk;
                        if (j >= 0 && j < Tf.nyi - 1 && i < Tf.nxi && k < Tf.nzi)
                            sum += Tf.TY(i, j, k);
                    }
                Tc.Ty[Tc.idxTy(I, J, K)] = sum;
            }
    // Z faces
    for (int K = 0; K < Tc.nzi - 1; ++K)
        for (int J = 0; J < Tc.nyi; ++J)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int k = 2 * K + 1;
                double sum = 0.0;
                for (int di = 0; di < 2; ++di)
                    for (int dj = 0; dj < 2; ++dj)
                    {
                        const int i = 2 * I + di;
                        const int j = 2 * J + dj;
                        if (k >= 0 && k < Tf.nzi - 1 && i < Tf.nxi && j < Tf.nyi)
                            sum += Tf.TZ(i, j, k);
                    }
                Tc.Tz[Tc.idxTz(I, J, K)] = sum;
            }
}

// Build solver hierarchy. Periodicity comes from BOTH the BC table and mesh.periodic.
// meshPerX/Y/Z override acts as a single source of truth when the mesh is periodic.
static void build_hierarchy(PPImpl& impl, MPI_Comm user_comm_in, const PBC& pbc, bool meshPerX,
                            bool meshPerY, bool meshPerZ)
{

    // Use the application's communicator as-is (borrowed; not owned here).
    if (impl.comm == MPI_COMM_NULL)
        impl.comm = user_comm_in;
    MPI_Comm comm = impl.comm;

#ifdef HAVE_MPI
    {
        int sz = -1, rk = -1;
        if (MPI_Comm_size(comm, &sz) != MPI_SUCCESS || sz < 1)
        {
            fprintf(stderr,
                    "FATAL: invalid MPI_Comm passed into PressurePoisson (size=%d). "
                    "Falling back to PETSC_COMM_WORLD.\n",
                    sz);
            comm = PETSC_COMM_WORLD;
            impl.comm = comm;
        }
        MPI_Comm_rank(comm, &rk); // also warms up the comm; ignore return for brevity
    }
#endif

    // ----- finest DMDA over interior cells -----
    const PetscInt nxi = impl.nxc_tot - 2 * impl.ng;
    const PetscInt nyi = impl.nyc_tot - 2 * impl.ng;
    const PetscInt nzi = impl.nzc_tot - 2 * impl.ng;

    // Periodicity from BCs OR mesh flags (mesh flags take precedence)
    auto is_periodic = [](const BcSpec* s) { return s && s->type == BcSpec::Type::periodic; };
    const bool perX = meshPerX || (is_periodic(pbc.W) && is_periodic(pbc.E));
    const bool perY = meshPerY || (is_periodic(pbc.S) && is_periodic(pbc.N));
    const bool perZ = meshPerZ || (is_periodic(pbc.B) && is_periodic(pbc.T));
    DMBoundaryType bx = perX ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType by = perY ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType bz = perZ ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;

    DMDACreate3d(comm, bx, by, bz, DMDA_STENCIL_STAR, nxi, nyi, nzi, PETSC_DECIDE, PETSC_DECIDE,
                 PETSC_DECIDE, 1, 1, nullptr, nullptr, nullptr, &impl.da_fine);

    DMSetUp(impl.da_fine);
    DMDASetInterpolationType(impl.da_fine, DMDA_Q0); // Q0 = cell-centered (FV) transfers

    // coarsen all dim <= 8
    std::vector<DM> fine_to_coarse;
    {
        DM cur = impl.da_fine;
        PetscObjectReference((PetscObject) cur);
        fine_to_coarse.push_back(cur);
        while (true)
        {
            DM c = nullptr;
            if (DMCoarsen(cur, comm, &c) || !c)
                break;
            PetscInt M = 0, N = 0, P = 0;
            DMDAGetInfo(c, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

            // Stop if each dim is small or total DOFs small
            if (M <= 8 && N <= 8 && P <= 8)
            {
                fine_to_coarse.push_back(c);
                break;
            }
            fine_to_coarse.push_back(c);
            cur = c;
        }
    }
    // reorder to PCMG convention: 0 = coarsest ... L-1 = finest
    impl.da_lvl.assign(fine_to_coarse.rbegin(), fine_to_coarse.rend());
    const PetscInt L = (PetscInt) impl.da_lvl.size();

    // ----- nullspace applicability (pure Neumann only) -----
    const bool allNeumann = (!pbc.W || norm_p_type(pbc.W->type) == BcSpec::Type::neumann) &&
                            (!pbc.E || norm_p_type(pbc.E->type) == BcSpec::Type::neumann) &&
                            (!pbc.S || norm_p_type(pbc.S->type) == BcSpec::Type::neumann) &&
                            (!pbc.N || norm_p_type(pbc.N->type) == BcSpec::Type::neumann) &&
                            (!pbc.B || norm_p_type(pbc.B->type) == BcSpec::Type::neumann) &&
                            (!pbc.T || norm_p_type(pbc.T->type) == BcSpec::Type::neumann);

    // Treat fully periodic box as having a constant nullspace too.
    const bool fullyPeriodic = perX && perY && perZ;
    impl.all_neumann = allNeumann || fullyPeriodic;
    // When all faces are (effective) Neumann, ∇·(β∇·) has a constant nullspace.
    // Otherwise (any Dirichlet present), the operator is SPD without a nullspace.
    // We’ll attach MatNullSpace only when allNeumann == true.

    // ----- per-level β (interior) -----
    impl.beta_lvl.clear();
    impl.beta_lvl.resize(L);
    impl.trans_lvl.clear();
    impl.trans_lvl.resize(L);
    // build fine β from host rho if available
    if (impl.varrho)
    {
        LevelBeta fine{};
        fine.nxi = nxi;
        fine.nyi = nyi;
        fine.nzi = nzi;
        fine.data.resize(std::size_t(nxi) * nyi * nzi);
        // read β=1/ρ from user's center array (with ghosts)
        // note: index mapping uses original totals with ng shift
        for (int k = 0; k < nzi; ++k)
            for (int j = 0; j < nyi; ++j)
                for (int i = 0; i < nxi; ++i)
                {
                    const int ic = i + impl.ng;
                    const int jc = j + impl.ng;
                    const int kc = k + impl.ng;
                    const std::size_t c =
                        std::size_t(ic) +
                        std::size_t(impl.nxc_tot) *
                            (std::size_t(jc) + std::size_t(impl.nyc_tot) * std::size_t(kc));
                    fine.data[std::size_t(i) +
                              std::size_t(nxi) *
                                  (std::size_t(j) + std::size_t(nyi) * std::size_t(k))] =
                        1.0; // filled later per-step
                }
        // placeholder; actual values filled each execute() before solve
        impl.beta_lvl[L - 1] = std::move(fine);
        // Set up containers for coarse β and transmissibilities (filled later per-step)
        for (int l = (int) L - 2; l >= 0; --l)
        {
            LevelBeta c{};
            PetscInt ci = 0, cj = 0, ck = 0;
            DMDAGetInfo(impl.da_lvl[l], nullptr, &ci, &cj, &ck, nullptr, nullptr, nullptr, nullptr,
                        nullptr, nullptr, nullptr, nullptr, nullptr);
            c.nxi = (int) ci;
            c.nyi = (int) cj;
            c.nzi = (int) ck;
            c.data.resize(std::size_t(ci) * cj * ck);
            impl.beta_lvl[l] = std::move(c);
            LevelTrans t{};
            t.nxi = (int) ci;
            t.nyi = (int) cj;
            t.nzi = (int) ck;
            impl.trans_lvl[l] = std::move(t);
        }
        // finest trans container as well
        impl.trans_lvl[L - 1].nxi = nxi;
        impl.trans_lvl[L - 1].nyi = nyi;
        impl.trans_lvl[L - 1].nzi = nzi;
    }

    // ----- build MatShell operators and KSP/PCMG -----
    impl.ctx_lvl.clear();
    impl.ctx_lvl.resize(L);
    impl.A_lvl.assign(L, nullptr);

    // outer KSP with MG preconditioner
    KSPCreate(comm, &impl.ksp);
    // For pure-Neumann (singular) use MINRES (left-preconditioned),
    // otherwise CG (PIPECG) with symmetric preconditioning.
    if (impl.all_neumann)
    {
        KSPSetType(impl.ksp, KSPMINRES);
        KSPSetPCSide(impl.ksp, PC_LEFT); // MINRES supports only left PC
        KSPSetNormType(impl.ksp,
                       KSP_NORM_PRECONDITIONED); // MINRES does NOT support UNPRECONDITIONED
    }
    else
    {
        KSPSetType(impl.ksp, KSPPIPECG);
        KSPSetPCSide(impl.ksp, PC_SYMMETRIC);
        KSPSetNormType(impl.ksp, KSP_NORM_UNPRECONDITIONED); // exact control of ||r||_2
    }

    PC pc;
    KSPGetPC(impl.ksp, &pc);
    PCSetType(pc, PCMG); // geometric multigrid
    PCMGSetLevels(pc, L, NULL);

    // set per-level interpolation (coarse -> fine), using Q0 on both DMs
    for (PetscInt l = 1; l < L; ++l)
    {
        // Enforce Q0 (cell-centered) on both levels
        DMDASetInterpolationType(impl.da_lvl[l - 1], DMDA_Q0);
        DMDASetInterpolationType(impl.da_lvl[l], DMDA_Q0);
        Mat P = NULL;
        DMCreateInterpolation(/*coarse*/ impl.da_lvl[l - 1], /*fine*/ impl.da_lvl[l], &P, NULL);
        PCMGSetInterpolation(pc, l, P);
        MatDestroy(&P);
    }

    // create per-level MatShell (for finest Amat) and record contexts
    for (PetscInt l = 0; l < L; ++l)
    {
        // sizes
        PetscInt nxi_l = 0, nyi_l = 0, nzi_l = 0;
        DMDAGetInfo(impl.da_lvl[l], nullptr, &nxi_l, &nyi_l, &nzi_l, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);

        // global sizes for MatShell
        Vec tmp;
        DMCreateGlobalVector(impl.da_lvl[l], &tmp);
        PetscInt M, N;
        VecGetSize(tmp, &M);
        N = M;
        VecDestroy(&tmp);

        // ctx
        auto ctx = std::make_unique<ShellCtx>();
        ctx->da = impl.da_lvl[l];
        ctx->nxi = (int) nxi_l;
        ctx->nyi = (int) nyi_l;
        ctx->nzi = (int) nzi_l;
        ctx->ng = impl.ng;
        ctx->nxc_tot = impl.nxc_tot;
        ctx->nyc_tot = impl.nyc_tot;
        ctx->nzc_tot = impl.nzc_tot;
        const int refine = (int) std::pow(2.0, (double) (L - 1 - l));
        ctx->dx = impl.dx * refine;
        ctx->dy = impl.dy * refine;
        ctx->dz = impl.dz * refine;
        ctx->trans = impl.varrho ? &impl.trans_lvl[l] : nullptr;
        ctx->beta = impl.varrho ? &impl.beta_lvl[l] : nullptr;
        ctx->use_const_beta = !impl.varrho;
        ctx->beta_const = impl.varrho ? 1.0 : (1.0 / impl.rho_const);
        ctx->pbc = pbc;

        // propagate periodicity to shell contexts
        ctx->perX = perX;
        ctx->perY = perY;
        ctx->perZ = perZ;

        Mat A;
        MatCreateShell(comm, PETSC_DECIDE, PETSC_DECIDE, M, N, (void*) ctx.get(), &A);
        MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) ShellMult);
        MatShellSetOperation(A, MATOP_GET_DIAGONAL,
                             (void (*)(void)) ShellGetDiagonal); // for Jacobi/Chebyshev
        // make symmetry explicit to solvers that may query A^T
        MatShellSetOperation(A, MATOP_MULT_TRANSPOSE, (void (*)(void)) ShellMult);
        MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
        if (!allNeumann)
            MatSetOption(A, MAT_SPD, PETSC_TRUE);

        impl.ctx_lvl[l] = std::move(ctx);
        impl.A_lvl[l] = A;
    }

    // MG cycle config (multiplicative V-cycle, 4 pre/post by default)
    PCMGSetType(pc, PC_MG_MULTIPLICATIVE);
    PCMGSetCycleType(pc, PC_MG_CYCLE_V);
    PCMGSetNumberSmooth(pc, 4);

    // connect fine DM to outer KSP (for viewers/options), but keep DM-owned by PCMG
    KSPSetDM(impl.ksp, impl.da_fine);
    KSPSetDMActive(impl.ksp, PETSC_FALSE);
    PCSetUseAmat(pc, PETSC_FALSE);            // use Pmat inside PC/MG
    PCMGSetGalerkin(pc, PC_MG_GALERKIN_PMAT); // build only coarse Pmats

    // Fetch and configure the actual per-level solvers
    {
        // Smoothers on levels 1..L-1 (JACOBI-CHEBYSHEV)

        for (PetscInt l = 1; l < L; ++l)
        {
            KSP kspl = NULL;
            PC pcl = NULL;
            PetscCallAbort(comm, PCMGGetSmoother(pc, l, &kspl));
            if (!kspl)
                continue;

            // Chebyshev + Jacobi (GET_DIAGONAL)
            PetscCallAbort(comm, KSPSetType(kspl, KSPCHEBYSHEV));
            PetscCallAbort(comm, KSPGetPC(kspl, &pcl));
            PetscCallAbort(comm, PCSetType(pcl, PCJACOBI));
            PCJacobiSetType(pcl, PC_JACOBI_DIAGONAL);
            PetscCallAbort(comm, KSPChebyshevEstEigSet(kspl, PETSC_DEFAULT, PETSC_DEFAULT,
                                                       PETSC_DEFAULT, PETSC_DEFAULT));
            PetscCallAbort(comm, KSPSetConvergenceTest(kspl, KSPConvergedSkip, NULL, NULL));
            PetscCallAbort(comm, KSPSetNormType(kspl, KSP_NORM_NONE));
            PetscCallAbort(comm,
                           KSPSetTolerances(kspl, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 3));
        }
    }

    // ---- Coarse level: force a direct solve by default (preonly + LU) ----
    // This pins the L=0 KSP/PC to a direct factorization unless the user overrides
    {
        KSP kspc = NULL;
        PC  pcc  = NULL;
        // Get the PCMG-managed coarse solver (level 0)
        PetscCallAbort(comm, PCMGGetCoarseSolve(pc, &kspc));
        if (kspc) {
            // No iterations on coarse grid: just apply the preconditioner once
            PetscCallAbort(comm, KSPSetType(kspc, KSPPREONLY));
            PetscCallAbort(comm, KSPGetPC(kspc, &pcc));
            // Prefer LU to work for both SPD and non-SPD coarse operators
            PetscCallAbort(comm, PCSetType(pcc, PCLU));
            if (!impl.all_neumann) PetscCallAbort(comm, PCSetType(pcc, PCCHOLESKY));
        }
    }

    // fine vectors
    DMCreateGlobalVector(impl.da_fine, &impl.x);
    DMCreateGlobalVector(impl.da_fine, &impl.b);

    // Assemble fine-level AIJ (PMAT) from the fine ShellCtx and wire operators:
    //  - Outer KSP: Amat = fine MatShell; Pmat = fine AIJ
    //  - PCMG: set only the finest level operators; Galerkin builds coarse PMATs
    if (impl.Afine_aij)
    {
        MatDestroy(&impl.Afine_aij);
        impl.Afine_aij = nullptr;
    }
    PetscCallAbort(comm, AssembleAIJFromShell(*impl.ctx_lvl[L - 1], &impl.Afine_aij));
    MatSetOption(impl.Afine_aij, MAT_SYMMETRIC, PETSC_TRUE);
    if (!impl.all_neumann)
    {
        MatSetOption(impl.Afine_aij, MAT_SPD, PETSC_TRUE);
    }

    // Outer KSP operators (Amat=MatShell, Pmat=AIJ)
    KSPSetOperators(impl.ksp, impl.A_lvl[L - 1], impl.Afine_aij);

    // Finest-level MG operators (Amat=MatShell, Pmat=AIJ); coarser levels via Galerkin(PMAT)
    PC pc_for_ops = NULL;
    KSPGetPC(impl.ksp, &pc_for_ops);
    PCMGSetOperators(pc_for_ops, L - 1, impl.A_lvl[L - 1], impl.Afine_aij);

    // ---------- Proper nullspace propagation (pure Neumann) ----------
    if (impl.all_neumann)
    {
        // One constant nullspace object reused everywhere
        if (!impl.ns_const)
        {
            PetscCallAbort(comm, MatNullSpaceCreate(comm, PETSC_TRUE, 0, NULL, &impl.ns_const));
        }
        // Attach to fine Amat and Pmat so KSP can handle the singular operator.
        PetscCallAbort(comm, MatSetNullSpace(impl.A_lvl[L - 1], impl.ns_const));
        PetscCallAbort(comm, MatSetNullSpace(impl.Afine_aij, impl.ns_const));
        // Help MG/AMG via near-nullspace on PMAT; coarse PMATs built by Galerkin will inherit.
        PetscCallAbort(comm, MatSetNearNullSpace(impl.Afine_aij, impl.ns_const));
    }

    KSPSetFromOptions(impl.ksp); // allow CLI/options to refine everything
}

// --------------------------- IAction impl -----------------------------------

PressurePoisson::~PressurePoisson()
{
    if (impl_)
        impl_->destroy();
}

PressurePoisson::PressurePoisson(double rho, double dx, double dy, double dz, int iters,
                                 void* mpi_comm, const BcTable& bcs)
    : rho_(rho), dx_(dx), dy_(dy), dz_(dz), iters_(iters), mpi_comm_(mpi_comm), bcs_(bcs)
{
    info_.name = "poisson";
    info_.phases = core::master::plugin::Phase::PostExchange;
    info_.access.reads = {{"u", 1}, {"v", 1}, {"w", 1}, {"p", 1}, {"rho", 1}};
    info_.access.writes = {{"p", 1}};
}

void PressurePoisson::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");
    const bool have_rho = fields.contains("rho");
    auto vrho = have_rho ? fields.view("rho") : vp; // shape match

    const int nxc_tot = vp.extents[0];
    const int nyc_tot = vp.extents[1];
    const int nzc_tot = vp.extents[2];

    const int nx = tile.box.hi[0] - tile.box.lo[0];
    const int ng = (nxc_tot - nx) / 2;

// Borrowed communicator: if null, fall back to self for single-rank use.
#ifdef HAVE_MPI
    // Robust unbox + validate. Prefer the app comm; fall back to PETSC_COMM_WORLD.
    MPI_Comm user_comm = PETSC_COMM_WORLD;
    if (mpi_comm_)
    {
        // Our public headers promise "boxed" MPI_Comm via mpi_box(); unbox it if so.
        MPI_Comm cand = mpi_unbox(mpi_comm_);
        user_comm = valid_comm_or(cand, PETSC_COMM_WORLD);
    }
    int sz = 1;
    MPI_Comm_size(user_comm, &sz);
    if (sz == 1)
        user_comm = PETSC_COMM_SELF;
#else
    MPI_Comm user_comm = PETSC_COMM_SELF;
#endif

    // rebuild hierarchy if geometry changed
    if (!impl_ || impl_->nxc_tot != nxc_tot || impl_->nyc_tot != nyc_tot ||
        impl_->nzc_tot != nzc_tot || impl_->ng != ng || impl_->dx != dx_ || impl_->dy != dy_ ||
        impl_->dz != dz_)
    {
        if (impl_)
            impl_->destroy();
        impl_.reset(new PPImpl{});
        impl_->nxc_tot = nxc_tot;
        impl_->nyc_tot = nyc_tot;
        impl_->nzc_tot = nzc_tot;
        impl_->ng = ng;
        impl_->dx = dx_;
        impl_->dy = dy_;
        impl_->dz = dz_;
        impl_->varrho = have_rho;
        impl_->rho_const = rho_;
        impl_->p_host.assign(std::size_t(nxc_tot) * nyc_tot * nzc_tot, 0.0);

        auto find_bc = [&](const char* k) -> const BcSpec*
        {
            auto it = bcs_.find(std::string("p.") + k);
            return (it == bcs_.end() ? nullptr : &it->second);
        };
        // Start from deck BCs, then override with mesh.periodic if requested
        PBC pbc{find_bc("west"),  find_bc("east"),   find_bc("south"),
                find_bc("north"), find_bc("bottom"), find_bc("top")};

        const bool perX = tile.mesh && tile.mesh->periodic[0];
        const bool perY = tile.mesh && tile.mesh->periodic[1];
        const bool perZ = tile.mesh && tile.mesh->periodic[2];
        static const BcSpec kPeriodic{BcSpec::Type::periodic, 0.0};
        if (perX)
        {
            pbc.W = &kPeriodic;
            pbc.E = &kPeriodic;
        }
        if (perY)
        {
            pbc.S = &kPeriodic;
            pbc.N = &kPeriodic;
        }
        if (perZ)
        {
            pbc.B = &kPeriodic;
            pbc.T = &kPeriodic;
        }

        // Gather mesh periodicity (authoritative when true)
        const bool meshPerX = (tile.mesh && tile.mesh->periodic[0]);
        const bool meshPerY = (tile.mesh && tile.mesh->periodic[1]);
        const bool meshPerZ = (tile.mesh && tile.mesh->periodic[2]);

        build_hierarchy(*impl_, user_comm, pbc, meshPerX, meshPerY, meshPerZ);
    }

    // ---------------- RHS = - (1/dt) * div(u*) plus BC contributions ----------------
    std::vector<double> div(std::size_t(nxc_tot) * nyc_tot * nzc_tot, 0.0);
    divergence_mac_c(
        static_cast<const double*>(vu.host_ptr), static_cast<const double*>(vv.host_ptr),
        static_cast<const double*>(vw.host_ptr), vu.extents[0], vu.extents[1], vu.extents[2],
        vv.extents[0], vv.extents[1], vv.extents[2], vw.extents[0], vw.extents[1], vw.extents[2],
        nxc_tot, nyc_tot, nzc_tot, ng, dx_, dy_, dz_, div.data());

    // fill fine-level β from current rho (if present), then rebuild fine T and coarsen T
    if (impl_->varrho)
    {
        LevelBeta& fine = impl_->beta_lvl.back();
        const int nxi = fine.nxi, nyi = fine.nyi, nzi = fine.nzi;
        for (int k = 0; k < nzi; ++k)
            for (int j = 0; j < nyi; ++j)
                for (int i = 0; i < nxi; ++i)
                {
                    const int ic = i + impl_->ng;
                    const int jc = j + impl_->ng;
                    const int kc = k + impl_->ng;
                    const std::size_t c =
                        std::size_t(ic) +
                        std::size_t(nxc_tot) *
                            (std::size_t(jc) + std::size_t(nyc_tot) * std::size_t(kc));
                    fine.data[std::size_t(i) +
                              std::size_t(nxi) *
                                  (std::size_t(j) + std::size_t(nyi) * std::size_t(k))] =
                        1.0 / static_cast<const double*>(vrho.host_ptr)[c];
                }
        // build fine face transmissibilities and coarsen (sum of subfaces)
        build_fine_trans_from_beta(fine, impl_->dx, impl_->dy, impl_->dz, impl_->trans_lvl.back());
        for (int l = (int) impl_->trans_lvl.size() - 2; l >= 0; --l)
        {
            average_beta_coarsen(impl_->beta_lvl[l + 1], impl_->beta_lvl[l]);
            // Conservative transmissibility aggregation
            LevelTrans& Tc = impl_->trans_lvl[l];
            Tc.nxi = impl_->beta_lvl[l].nxi;
            Tc.nyi = impl_->beta_lvl[l].nyi;
            Tc.nzi = impl_->beta_lvl[l].nzi;
            coarsen_trans_2x(impl_->trans_lvl[l + 1], Tc);
        }
        // mark diagonals dirty (β changed)
        for (auto& ctx : impl_->ctx_lvl)
            ctx->diag_built = false;

        // --- Rebuild fine assembled operator (PMAT) and rewire finest level (reflects new β/T) ---
        if (impl_->Afine_aij)
        {
            MatDestroy(&impl_->Afine_aij);
            impl_->Afine_aij = nullptr;
        }
        PetscCallAbort(impl_->comm,
                       AssembleAIJFromShell(*impl_->ctx_lvl.back(), &impl_->Afine_aij));

        // Re-attach the same constant nullspace after the fine Pmat is rebuilt
        if (impl_->all_neumann)
        {
            // Ensure the reusable nullspace exists (should already be set in build_hierarchy)
            if (!impl_->ns_const)
            {
                PetscCallAbort(impl_->comm, MatNullSpaceCreate(impl_->comm, PETSC_TRUE, 0, NULL,
                                                               &impl_->ns_const));
            }
            PetscCallAbort(impl_->comm, MatSetNullSpace(impl_->A_lvl.back(), impl_->ns_const));
            PetscCallAbort(impl_->comm, MatSetNullSpace(impl_->Afine_aij, impl_->ns_const));
            // Keep near-nullspace too so Galerkin coarse Pmats remain aware of the constant mode
            PetscCallAbort(impl_->comm, MatSetNearNullSpace(impl_->Afine_aij, impl_->ns_const));
        }

        // Update operators on outer KSP and MG finest (coarse PMATs will be regen by Galerkin)
        KSPSetOperators(impl_->ksp, impl_->A_lvl.back(), impl_->Afine_aij);
        PC pc_refresh = NULL;
        KSPGetPC(impl_->ksp, &pc_refresh);
        const PetscInt L = (PetscInt) impl_->da_lvl.size();
        PCMGSetOperators(pc_refresh, L - 1, impl_->A_lvl.back(), impl_->Afine_aij);
    }

    // Build RHS on fine DMDA (owned box only) with BC shifts
    {
        auto find_bc = [&](const char* k) -> const BcSpec*
        {
            auto it = bcs_.find(std::string("p.") + k);
            return (it == bcs_.end() ? nullptr : &it->second);
        };
        PBC pbc{find_bc("west"),  find_bc("east"),   find_bc("south"),
                find_bc("north"), find_bc("bottom"), find_bc("top")};

        PetscInt xs, ys, zs, xm, ym, zm;
        DMDAGetCorners(impl_->da_fine, &xs, &ys, &zs, &xm, &ym, &zm);

        PetscScalar*** barr;
        DMDAVecGetArray(impl_->da_fine, impl_->b, &barr);

        // Face-T accessors on fine level
        const LevelTrans* Tf = impl_->varrho ? &impl_->trans_lvl.back() : nullptr;
        const LevelBeta* Bf = impl_->varrho ? &impl_->beta_lvl.back() : nullptr;
        const double invV = 1.0 / (dx_ * dy_ * dz_);
        const double Tcx = (impl_->varrho ? 0.0 : (1.0 / rho_) * (dy_ * dz_ / dx_));
        const double Tcy = (impl_->varrho ? 0.0 : (1.0 / rho_) * (dx_ * dz_ / dy_));
        const double Tcz = (impl_->varrho ? 0.0 : (1.0 / rho_) * (dx_ * dy_ / dz_));
        auto Tw = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (i > 0)
                    return Tf->TX(i - 1, j, k);
                const double b = Bf->B(i, j, k);
                return b * (dy_ * dz_ / dx_);
            }
            else
                return Tcx;
        };
        auto Te = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (i < (nxc_tot - 2 * ng) - 1)
                    return Tf->TX(i, j, k);
                const double b = Bf->B(i, j, k);
                return b * (dy_ * dz_ / dx_);
            }
            else
                return Tcx;
        };
        auto Ts = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (j > 0)
                    return Tf->TY(i, j - 1, k);
                const double b = Bf->B(i, j, k);
                return b * (dx_ * dz_ / dy_);
            }
            else
                return Tcy;
        };
        auto Tn = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (j < (nyc_tot - 2 * ng) - 1)
                    return Tf->TY(i, j, k);
                const double b = Bf->B(i, j, k);
                return b * (dx_ * dz_ / dy_);
            }
            else
                return Tcy;
        };
        auto Tb = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (k > 0)
                    return Tf->TZ(i, j, k - 1);
                const double b = Bf->B(i, j, k);
                return b * (dx_ * dy_ / dz_);
            }
            else
                return Tcz;
        };
        auto Tt = [&](int i, int j, int k)
        {
            if (Tf)
            {
                if (k < (nzc_tot - 2 * ng) - 1)
                    return Tf->TZ(i, j, k);
                const double b = Bf->B(i, j, k);
                return b * (dx_ * dy_ / dz_);
            }
            else
                return Tcz;
        };

        for (int k = zs; k < zs + zm; ++k)
            for (int j = ys; j < ys + ym; ++j)
                for (int i = xs; i < xs + xm; ++i)
                {
                    const std::size_t c =
                        std::size_t(i + ng) +
                        std::size_t(nxc_tot) *
                            (std::size_t(j + ng) + std::size_t(nyc_tot) * std::size_t(k + ng));
                    // A is -div(β∇p) (SPD). Projection gives -div(β∇p) = -(1/dt) div(u*).
                    double rhs = -div[c] / dt;

                    // Use the finest-level periodic flags
                    const ShellCtx* fctx = impl_->ctx_lvl.back().get();
                    const bool perX = fctx->perX, perY = fctx->perY, perZ = fctx->perZ;
                    const bool hasW = (i > 0) || perX;
                    const bool hasE = (i < (nxc_tot - 2 * ng) - 1) || perX;
                    const bool hasS = (j > 0) || perY;
                    const bool hasN = (j < (nyc_tot - 2 * ng) - 1) || perY;
                    const bool hasB = (k > 0) || perZ;
                    const bool hasT = (k < (nzc_tot - 2 * ng) - 1) || perZ;

                    const double te = hasE ? Te(i, j, k) : 0.0;
                    const double tw = hasW ? Tw(i, j, k) : 0.0;
                    const double tn = hasN ? Tn(i, j, k) : 0.0;
                    const double ts = hasS ? Ts(i, j, k) : 0.0;
                    const double tt = hasT ? Tt(i, j, k) : 0.0;
                    const double tb = hasB ? Tb(i, j, k) : 0.0;

                    // Dirichlet: RHS += (1/V) * T_face * value (honor axis periodicity)
                    auto add_dirichlet = [&](const BcSpec* s, bool per_axis, double coef, double value)
                    {
                        if (!per_axis && s && norm_p_type(s->type) == BcSpec::Type::dirichlet)
                            rhs += invV * coef * value;
                    };
                    // Neumann: flux = -β A * (∂p/∂n) ⇒ RHS += (1/V) * (± β A * g) = (± invV *
                    // (T_face * h) * g)
                    auto add_neumann = [&](const BcSpec* s, double Tface, double h, int sgn)
                    {
                        if (s && norm_p_type(s->type) == BcSpec::Type::neumann)
                            rhs += sgn * invV * (Tface * h) * s->value;
                    };
                    // WEST/EAST
                    if (!hasW && !perX)
                    {
                        add_dirichlet(pbc.W, perX, Tw(i,j,k), pbc.W ? pbc.W->value : 0.0);
                        add_neumann(pbc.W, Tw(i, j, k), dx_, -1);
                    }
                    if (!hasE && !perX)
                    {
                        add_dirichlet(pbc.E, perX, Te(i,j,k), pbc.E ? pbc.E->value : 0.0);
                        add_neumann(pbc.E, Te(i, j, k), dx_, +1);
                    }
                    // SOUTH/NORTH
                    if (!hasS && !perY)
                    {
                        add_dirichlet(pbc.S, perY, Ts(i,j,k), pbc.S ? pbc.S->value : 0.0);
                        add_neumann(pbc.S, Ts(i, j, k), dy_, -1);
                    }
                    if (!hasN && !perY)
                    {
                        add_dirichlet(pbc.N, perY, Tn(i,j,k), pbc.N ? pbc.N->value : 0.0);
                        add_neumann(pbc.N, Tn(i, j, k), dy_, +1);
                    }
                    // BOTTOM/TOP
                    if (!hasB && !perZ)
                    {
                        add_dirichlet(pbc.B, perZ, Tb(i,j,k), pbc.B ? pbc.B->value : 0.0);
                        add_neumann(pbc.B, Tb(i, j, k), dz_, -1);
                    }
                    if (!hasT && !perZ)
                    {
                        add_dirichlet(pbc.T, perZ, Tt(i,j,k), pbc.T ? pbc.T->value : 0.0);
                        add_neumann(pbc.T, Tt(i, j, k), dz_, +1);
                    }

                    barr[k][j][i] = rhs;
                }
        DMDAVecRestoreArray(impl_->da_fine, impl_->b, &barr);
    }

    // initial guess: previous p
    {
        PetscInt xs, ys, zs, xm, ym, zm;
        DMDAGetCorners(impl_->da_fine, &xs, &ys, &zs, &xm, &ym, &zm);

        PetscScalar*** xarr;
        DMDAVecGetArray(impl_->da_fine, impl_->x, &xarr);
        for (int k = zs; k < zs + zm; ++k)
            for (int j = ys; j < ys + ym; ++j)
                for (int i = xs; i < xs + xm; ++i)
                {
                    const std::size_t c =
                        std::size_t(i + ng) +
                        std::size_t(nxc_tot) *
                            (std::size_t(j + ng) + std::size_t(nyc_tot) * std::size_t(k + ng));
                    xarr[k][j][i] = impl_->p_host[c];
                }
        DMDAVecRestoreArray(impl_->da_fine, impl_->x, &xarr);
    }

    // Remove constant component from BOTH b and x when a nullspace is present (pure Neumann).
    // This keeps the singular system well-posed for Krylov and avoids stalls.
    {
        // Use the vectors' communicator (same as impl.comm)
        MPI_Comm vcomm = impl_->comm;
        PetscCallAbort(vcomm, PetscObjectGetComm((PetscObject) impl_->b, &vcomm));

        MatNullSpace nsp = NULL;
        PetscCallAbort(vcomm, MatGetNullSpace(impl_->A_lvl.back(), &nsp)); // may be NULL

        if (nsp)
        {
            // Ensure b ⟂ 1 and (if used) x0 ⟂ 1.
            PetscCallAbort(vcomm, MatNullSpaceRemove(nsp, impl_->b));

            PetscBool guessNonzero = PETSC_FALSE;
            PetscCallAbort(vcomm, KSPGetInitialGuessNonzero(impl_->ksp, &guessNonzero));
            if (guessNonzero)
            {
                PetscCallAbort(vcomm, MatNullSpaceRemove(nsp, impl_->x));
            }
        }
        else
        {
            // Not a pure-Neumann problem: safest is a zero initial guess to avoid
            // measuring an odd initial residual against a tiny RHS.
            PetscCallAbort(vcomm, VecZeroEntries(impl_->x));
            PetscCallAbort(vcomm, KSPSetInitialGuessNonzero(impl_->ksp, PETSC_FALSE));
        }
    }

    // If the projected RHS is (numerically) zero, the solution is trivially zero.
    // Short-circuit to avoid DTOL-at-0 and pointless work.
    {
        PetscReal nb = 0.0;
        MPI_Comm vcomm = impl_->comm;
        PetscCallAbort(vcomm, PetscObjectGetComm((PetscObject) impl_->b, &vcomm));
        PetscCallAbort(vcomm, VecNorm(impl_->b, NORM_2, &nb));
        // Absolute floor
        if (nb < 1e-12)
        {
            PetscCallAbort(vcomm, VecZeroEntries(impl_->x));
        }
    }

    // ---- set outer KSP tolerances from div_tol & dt ----
    // Stop when ||r||_2 <= div_tol / dt  (since div^{n+1} = dt * r)
    {
        PetscReal nb = 0.0; // nb = ||b||_2
        VecNorm(impl_->b, NORM_2, &nb);
        const PetscReal atol = user_div_tol_ / dt; // absolute target: ||r||_2 <= div_tol/dt
        PetscReal rtol = 0.0;
        if (nb > 0.0)
        {
            const PetscReal cand = atol / nb; // desired relative tol
            // If cand >= 1, relative test is meaningless; rely on absolute tol instead.
            rtol = (cand < 1.0) ? cand : 0.0;
        }
        // Large dtol to avoid false divergence on tiny RHS; cap iters reasonably.
        KSPSetTolerances(impl_->ksp, rtol, atol, 1e50, 200);
    }
    // solve
    MPI_Comm pcomm;
    PetscObjectGetComm((PetscObject) impl_->ksp, &pcomm);
    PetscCallAbort(pcomm, KSPSolve(impl_->ksp, impl_->b, impl_->x));

    // Gauge fix for pure Neumann (zero-mean)
    if (impl_->all_neumann)
    {

        // --- enforce zero-mean gauge on interior cells ---
        PetscReal sum;
        PetscInt N;
        MPI_Comm vcomm = impl_->comm;
        PetscCallAbort(vcomm, PetscObjectGetComm((PetscObject) impl_->x, &vcomm));
        PetscCallAbort(vcomm, VecSum(impl_->x, &sum));
        PetscCallAbort(vcomm, VecGetSize(impl_->x, &N));
        const PetscScalar mean = (N ? sum / (PetscReal) N : 0.0);
        PetscCallAbort(vcomm, VecShift(impl_->x, -mean)); // subtract global mean
    }

    // copy solution back to p_host (owned box)
    {
        PetscInt xs, ys, zs, xm, ym, zm;
        DMDAGetCorners(impl_->da_fine, &xs, &ys, &zs, &xm, &ym, &zm);
        const PetscScalar*** xarr;
        DMDAVecGetArrayRead(impl_->da_fine, impl_->x, &xarr);
        for (int k = zs; k < zs + zm; ++k)
            for (int j = ys; j < ys + ym; ++j)
                for (int i = xs; i < xs + xm; ++i)
                {
                    const std::size_t c =
                        std::size_t(i + ng) +
                        std::size_t(nxc_tot) *
                            (std::size_t(j + ng) + std::size_t(nyc_tot) * std::size_t(k + ng));
                    impl_->p_host[c] = xarr[k][j][i];
                }
        DMDAVecRestoreArrayRead(impl_->da_fine, impl_->x, &xarr);
    }

    // SYNC: copy the solved pressure from our host scratch into the actual
    // field catalog storage so downstream stages (halo exchange + corrector)
    // see the updated 'p'.
    {
        auto vp_sync = fields.view("p");
        const std::size_t nTot = std::size_t(nxc_tot) * nyc_tot * nzc_tot;
        std::memcpy(vp_sync.host_ptr, impl_->p_host.data(), nTot * sizeof(double));
    }
}

// factory
std::shared_ptr<core::master::plugin::IAction> make_poisson(const core::master::plugin::KV& kv,
                                                            const core::master::RunContext& rc)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };
    const double rho = std::stod(get("rho", "1.0"));
    const double dx = std::stod(get("dx", "1.0"));
    const double dy = std::stod(get("dy", "1.0"));
    const double dz = std::stod(get("dz", "1.0"));
    const int iters = std::stoi(get("iters", "50")); // retained for compatibility
    const Params p = parse_params(kv);
    auto obj = std::make_shared<PressurePoisson>(rho, dx, dy, dz, iters, rc.mpi_comm, p.bcs);
    // propagate div_tol from user params if present (fallback 1e-7)
    obj->set_div_tol(std::stod(get("div_tol", "1e-7")));
    return obj;
}

} // namespace fluids
