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
#include "master/Log.hpp"
#include "master/Views.hpp"
#include "MacOps.hpp"

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

static inline MPI_Comm valid_comm_or(MPI_Comm cand, MPI_Comm fallback)
{
    if (cand == MPI_COMM_NULL)
        return fallback;
    int sz = -1, rc = MPI_Comm_size(cand, &sz);
    if (rc == MPI_SUCCESS && sz > 0)
        return cand;
    return fallback;
}

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
    // Boundary face conductances (physical boundaries only; periodic axes ignore these):
    //  WEST/EAST slabs: size = nyi*nzi
    //  SOUTH/NORTH slabs: size = nxi*nzi
    //  BOTTOM/TOP slabs: size = nxi*nyi
    std::vector<double> Tw_bnd, Te_bnd;
    std::vector<double> Ts_bnd, Tn_bnd;
    std::vector<double> Tb_bnd, Tt_bnd;
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

    // boundary accessors
    inline std::size_t idxWE(int j, int k) const
    {
        return std::size_t(j) + std::size_t(nyi) * std::size_t(k);
    }
    inline std::size_t idxSN(int i, int k) const
    {
        return std::size_t(i) + std::size_t(nxi) * std::size_t(k);
    }
    inline std::size_t idxBT(int i, int j) const
    {
        return std::size_t(i) + std::size_t(nxi) * std::size_t(j);
    }
    inline double TWb(int j, int k) const { return Tw_bnd[idxWE(j, k)]; }
    inline double TEb(int j, int k) const { return Te_bnd[idxWE(j, k)]; }
    inline double TSb(int i, int k) const { return Ts_bnd[idxSN(i, k)]; }
    inline double TNb(int i, int k) const { return Tn_bnd[idxSN(i, k)]; }
    inline double TBb(int i, int j) const { return Tb_bnd[idxBT(i, j)]; }
    inline double TTb(int i, int j) const { return Tt_bnd[idxBT(i, j)]; }
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
    const bool constBeta = ctx->use_const_beta || (ctx->beta == nullptr);

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

                // Decide neighbors against *global* domain edges.
                // DMDA ghosts provide inter-rank neighbors automatically.
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
                // For interior faces, use precomputed Tx/Ty/Tz if present, else constant-β
                // TXc/TYc/TZc. For boundary faces, fall back to local β * (A/Δ).
                const double bC = constBeta ? ctx->beta_const : ctx->beta->B(i, j, k);

                double Tw = 0.0, Te = 0.0, Ts = 0.0, Tn = 0.0, Tb = 0.0, Tt = 0.0;
                // X faces (periodic seams use harmonic mean across wrap; physical boundaries use
                // pre-coarsened boundary T)
                if (i > 0)
                {
                    if (useTx)
                        Tw = ctx->trans->TX(i - 1, j, k);
                    else
                        Tw = constBeta ? TXc
                                       : (hmean(ctx->beta->B(i - 1, j, k), bC) *
                                          (ctx->dy * ctx->dz / ctx->dx));
                }
                else if (ctx->perX)
                {
                    if (constBeta)
                        Tw = TXc;
                    else
                    {
                        const double bw = ctx->beta->B(ctx->nxi - 1, j, k);
                        const double bh = hmean(bw, bC);
                        Tw = bh * (ctx->dy * ctx->dz / ctx->dx);
                    }
                }
                else
                {
                    Tw = useTx ? ctx->trans->TWb(j, k) : (bC * (ctx->dy * ctx->dz / ctx->dx));
                }
                if (i < ctx->nxi - 1)
                {
                    if (useTx)
                        Te = ctx->trans->TX(i, j, k);
                    else
                        Te = constBeta ? TXc
                                       : (hmean(ctx->beta->B(i + 1, j, k), bC) *
                                          (ctx->dy * ctx->dz / ctx->dx));
                }
                else if (ctx->perX)
                {
                    if (constBeta)
                        Te = TXc;
                    else
                    {
                        const double be = ctx->beta->B(0, j, k);
                        const double bh = hmean(be, bC);
                        Te = bh * (ctx->dy * ctx->dz / ctx->dx);
                    }
                }
                else
                {
                    Te = useTx ? ctx->trans->TEb(j, k) : (bC * (ctx->dy * ctx->dz / ctx->dx));
                }
                // Y faces
                if (j > 0)
                {
                    if (useTy)
                        Ts = ctx->trans->TY(i, j - 1, k);
                    else
                        Ts = constBeta ? TYc
                                       : (hmean(ctx->beta->B(i, j - 1, k), bC) *
                                          (ctx->dx * ctx->dz / ctx->dy));
                }
                else if (ctx->perY)
                {
                    if (constBeta)
                        Ts = TYc;
                    else
                    {
                        const double bs = ctx->beta->B(i, ctx->nyi - 1, k);
                        const double bh = hmean(bs, bC);
                        Ts = bh * (ctx->dx * ctx->dz / ctx->dy);
                    }
                }
                else
                {
                    Ts = useTy ? ctx->trans->TSb(i, k) : (bC * (ctx->dx * ctx->dz / ctx->dy));
                }
                if (j < ctx->nyi - 1)
                {
                    if (useTy)
                        Tn = ctx->trans->TY(i, j, k);
                    else
                        Tn = constBeta ? TYc
                                       : (hmean(ctx->beta->B(i, j + 1, k), bC) *
                                          (ctx->dx * ctx->dz / ctx->dy));
                }
                else if (ctx->perY)
                {
                    if (constBeta)
                        Tn = TYc;
                    else
                    {
                        const double bn = ctx->beta->B(i, 0, k);
                        const double bh = hmean(bn, bC);
                        Tn = bh * (ctx->dx * ctx->dz / ctx->dy);
                    }
                }
                else
                {
                    Tn = useTy ? ctx->trans->TNb(i, k) : (bC * (ctx->dx * ctx->dz / ctx->dy));
                }
                // Z faces
                if (k > 0)
                {
                    if (useTz)
                        Tb = ctx->trans->TZ(i, j, k - 1);
                    else
                        Tb = constBeta ? TZc
                                       : (hmean(ctx->beta->B(i, j, k - 1), bC) *
                                          (ctx->dx * ctx->dy / ctx->dz));
                }
                else if (ctx->perZ)
                {
                    if (constBeta)
                        Tb = TZc;
                    else
                    {
                        const double bb = ctx->beta->B(i, j, ctx->nzi - 1);
                        const double bh = hmean(bb, bC);
                        Tb = bh * (ctx->dx * ctx->dy / ctx->dz);
                    }
                }
                else
                {
                    Tb = useTz ? ctx->trans->TBb(i, j) : (bC * (ctx->dx * ctx->dy / ctx->dz));
                }
                if (k < ctx->nzi - 1)
                {
                    if (useTz)
                        Tt = ctx->trans->TZ(i, j, k);
                    else
                        Tt = constBeta ? TZc
                                       : (hmean(ctx->beta->B(i, j, k + 1), bC) *
                                          (ctx->dx * ctx->dy / ctx->dz));
                }
                else if (ctx->perZ)
                {
                    if (constBeta)
                        Tt = TZc;
                    else
                    {
                        const double bt = ctx->beta->B(i, j, 0);
                        const double bh = hmean(bt, bC);
                        Tt = bh * (ctx->dx * ctx->dy / ctx->dz);
                    }
                }
                else
                {
                    Tt = useTz ? ctx->trans->TTb(i, j) : (bC * (ctx->dx * ctx->dy / ctx->dz));
                }

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
    const bool constBeta = ctx->use_const_beta || (ctx->beta == nullptr);

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
                // Use *global* tests; inter-rank neighbors arrive via DMDA ghosts.
                const bool hasW = (i > 0) || ctx->perX;
                const bool hasE = (i < ctx->nxi - 1) || ctx->perX;
                const bool hasS = (j > 0) || ctx->perY;
                const bool hasN = (j < ctx->nyi - 1) || ctx->perY;
                const bool hasB = (k > 0) || ctx->perZ;
                const bool hasT = (k < ctx->nzi - 1) || ctx->perZ;
                const double bC = constBeta ? ctx->beta_const : ctx->beta->B(i, j, k);
                double Tw, Te, Ts, Tn, Tb, Tt;
                // X faces (match ShellMult logic)
                if (i > 0)
                    Tw = useTx ? ctx->trans->TX(i - 1, j, k)
                               : (constBeta ? TXc
                                            : (hmean(ctx->beta->B(i - 1, j, k), bC) *
                                               (ctx->dy * ctx->dz / ctx->dx)));
                else if (ctx->perX)
                    Tw = constBeta ? TXc
                                   : (hmean(ctx->beta->B(ctx->nxi - 1, j, k), bC) *
                                      (ctx->dy * ctx->dz / ctx->dx));
                else
                    Tw = useTx ? ctx->trans->TWb(j, k) : (bC * (ctx->dy * ctx->dz / ctx->dx));
                if (i < ctx->nxi - 1)
                    Te = useTx ? ctx->trans->TX(i, j, k)
                               : (constBeta ? TXc
                                            : (hmean(ctx->beta->B(i + 1, j, k), bC) *
                                               (ctx->dy * ctx->dz / ctx->dx)));
                else if (ctx->perX)
                    Te = constBeta
                             ? TXc
                             : (hmean(ctx->beta->B(0, j, k), bC) * (ctx->dy * ctx->dz / ctx->dx));
                else
                    Te = useTx ? ctx->trans->TEb(j, k) : (bC * (ctx->dy * ctx->dz / ctx->dx));
                // Y faces
                if (j > 0)
                    Ts = useTy ? ctx->trans->TY(i, j - 1, k)
                               : (constBeta ? TYc
                                            : (hmean(ctx->beta->B(i, j - 1, k), bC) *
                                               (ctx->dx * ctx->dz / ctx->dy)));
                else if (ctx->perY)
                    Ts = constBeta ? TYc
                                   : (hmean(ctx->beta->B(i, ctx->nyi - 1, k), bC) *
                                      (ctx->dx * ctx->dz / ctx->dy));
                else
                    Ts = useTy ? ctx->trans->TSb(i, k) : (bC * (ctx->dx * ctx->dz / ctx->dy));
                if (j < ctx->nyi - 1)
                    Tn = useTy ? ctx->trans->TY(i, j, k)
                               : (constBeta ? TYc
                                            : (hmean(ctx->beta->B(i, j + 1, k), bC) *
                                               (ctx->dx * ctx->dz / ctx->dy)));
                else if (ctx->perY)
                    Tn = constBeta
                             ? TYc
                             : (hmean(ctx->beta->B(i, 0, k), bC) * (ctx->dx * ctx->dz / ctx->dy));
                else
                    Tn = useTy ? ctx->trans->TNb(i, k) : (bC * (ctx->dx * ctx->dz / ctx->dy));
                // Z faces
                if (k > 0)
                    Tb = useTz ? ctx->trans->TZ(i, j, k - 1)
                               : (constBeta ? TZc
                                            : (hmean(ctx->beta->B(i, j, k - 1), bC) *
                                               (ctx->dx * ctx->dy / ctx->dz)));
                else if (ctx->perZ)
                    Tb = constBeta ? TZc
                                   : (hmean(ctx->beta->B(i, j, ctx->nzi - 1), bC) *
                                      (ctx->dx * ctx->dy / ctx->dz));
                else
                    Tb = useTz ? ctx->trans->TBb(i, j) : (bC * (ctx->dx * ctx->dy / ctx->dz));
                if (k < ctx->nzi - 1)
                    Tt = useTz ? ctx->trans->TZ(i, j, k)
                               : (constBeta ? TZc
                                            : (hmean(ctx->beta->B(i, j, k + 1), bC) *
                                               (ctx->dx * ctx->dy / ctx->dz)));
                else if (ctx->perZ)
                    Tt = constBeta
                             ? TZc
                             : (hmean(ctx->beta->B(i, j, 0), bC) * (ctx->dx * ctx->dy / ctx->dz));
                else
                    Tt = useTz ? ctx->trans->TTb(i, j) : (bC * (ctx->dx * ctx->dy / ctx->dz));

                double diag = 0.0;
                if (hasE)
                    diag += Te;
                if (hasW)
                    diag += Tw;
                if (hasN)
                    diag += Tn;
                if (hasS)
                    diag += Ts;
                if (hasT)
                    diag += Tt;
                if (hasB)
                    diag += Tb;

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

// Ensure Chebyshev smoothers are (re)tuned after operator/PMAT updates.
// This re-runs PETSc's eigenvalue estimator so the Chebyshev bounds track
// changes in the preconditioned operator (crucial when coefficients/RHS vary).
static void retune_chebyshev_on_levels(KSP outerKSP)
{
    MPI_Comm comm = PETSC_COMM_SELF;
    if (outerKSP)
        PetscObjectGetComm((PetscObject) outerKSP, &comm);
    PC pc = NULL;
    PetscCallAbort(comm, KSPGetPC(outerKSP, &pc));
    PetscInt nlev = 0;
    PetscCallAbort(comm, PCMGGetLevels(pc, &nlev));
    for (PetscInt l = 1; l < nlev; ++l)
    {
        KSP kspl = NULL;
        PetscCallAbort(comm, PCMGGetSmoother(pc, l, &kspl));
        if (!kspl)
            continue;
        // Keep Jacobi PC as set elsewhere; just (re)enable Chebyshev + auto eig est.
        PetscCallAbort(comm, KSPSetType(kspl, KSPCHEBYSHEV));
        PetscCallAbort(comm, KSPChebyshevEstEigSet(kspl, PETSC_DEFAULT, PETSC_DEFAULT,
                                                   PETSC_DEFAULT, PETSC_DEFAULT));
        PetscCallAbort(comm,
                       KSPSetTolerances(kspl, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 2));
        PetscCallAbort(comm, KSPSetConvergenceTest(kspl, KSPConvergedSkip, NULL, NULL));
        PetscCallAbort(comm, KSPSetNormType(kspl, KSP_NORM_NONE));
    }
}

// If *Aio is nullptr, create & preallocate. Otherwise, reuse the existing pattern and just
// overwrite values.
static PetscErrorCode AssembleAIJFromShell(const ShellCtx& ctx, Mat* Aio)
{
    const bool useTx = (ctx.trans && !ctx.trans->Tx.empty());
    const bool useTy = (ctx.trans && !ctx.trans->Ty.empty());
    const bool useTz = (ctx.trans && !ctx.trans->Tz.empty());
    PetscFunctionBegin;
    DM da = ctx.da;

    const PetscInt nxi = ctx.nxi, nyi = ctx.nyi, nzi = ctx.nzi;

    // Create on first use; otherwise reuse the existing sparsity pattern
    if (!*Aio)
    {
        // Create a properly preallocated AIJ using the DMDA stencil (STAR, width=1)
        PetscCall(DMSetMatType(da, MATAIJ));
        PetscCall(DMCreateMatrix(da, Aio)); // *Aio gets AIJ with 7-pt prealloc
        PetscCall(MatSetOption(*Aio, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
    }
    else
    {
        // Zero existing entries; keep structure
        PetscCall(MatZeroEntries(*Aio));
    }

    // From this point on, always insert into *Aio (created or reused)
    Mat Afill = *Aio;

    const double invV = 1.0 / (ctx.dx * ctx.dy * ctx.dz);

    // Work with both owned corners and ghosted corners so we can insert using *local* stencil
    // (avoids ISLocalToGlobalMappingApply() overflow on periodic wraps).
    PetscInt xs, ys, zs, xm, ym, zm;
    PetscCall(DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm));
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    PetscCall(DMDAGetGhostCorners(da, &gxs, &gys, &gzs, &gxm, &gym, &gzm));

    // Accept subset/off-proc inserts and allow non-prealloc'd slots (safe with DMDA pattern).
    PetscCall(MatSetOption(Afill, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE));
    PetscCall(MatSetOption(Afill, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE));
    const bool constBeta = ctx.use_const_beta || (ctx.beta == nullptr);
    auto Tw = [&](int i, int j, int k)
    {
        if (i > 0)
            return useTx ? ctx.trans->TX(i - 1, j, k)
                         : (constBeta ? ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx)
                                      : (hmean(ctx.beta->B(i - 1, j, k), ctx.beta->B(i, j, k)) *
                                         (ctx.dy * ctx.dz / ctx.dx)));
        if (ctx.perX)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx);
            const double bC = ctx.beta->B(i, j, k);
            const double bw = ctx.beta->B(ctx.nxi - 1, j, k);
            return hmean(bw, bC) * (ctx.dy * ctx.dz / ctx.dx);
        }
        return useTx ? ctx.trans->TWb(j, k)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dy * ctx.dz / ctx.dx));
    };
    auto Te = [&](int i, int j, int k)
    {
        if (i < ctx.nxi - 1)
            return useTx ? ctx.trans->TX(i, j, k)
                         : (constBeta ? ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx)
                                      : (hmean(ctx.beta->B(i + 1, j, k), ctx.beta->B(i, j, k)) *
                                         (ctx.dy * ctx.dz / ctx.dx)));
        if (ctx.perX)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dy * ctx.dz / ctx.dx);
            const double bC = ctx.beta->B(i, j, k);
            const double be = ctx.beta->B(0, j, k);
            return hmean(be, bC) * (ctx.dy * ctx.dz / ctx.dx);
        }
        return useTx ? ctx.trans->TEb(j, k)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dy * ctx.dz / ctx.dx));
    };
    auto Ts = [&](int i, int j, int k)
    {
        if (j > 0)
            return useTy ? ctx.trans->TY(i, j - 1, k)
                         : (constBeta ? ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy)
                                      : (hmean(ctx.beta->B(i, j - 1, k), ctx.beta->B(i, j, k)) *
                                         (ctx.dx * ctx.dz / ctx.dy)));
        if (ctx.perY)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy);
            const double bC = ctx.beta->B(i, j, k);
            const double bs = ctx.beta->B(i, ctx.nyi - 1, k);
            return hmean(bs, bC) * (ctx.dx * ctx.dz / ctx.dy);
        }
        return useTy ? ctx.trans->TSb(i, k)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dx * ctx.dz / ctx.dy));
    };
    auto Tn = [&](int i, int j, int k)
    {
        if (j < ctx.nyi - 1)
            return useTy ? ctx.trans->TY(i, j, k)
                         : (constBeta ? ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy)
                                      : (hmean(ctx.beta->B(i, j + 1, k), ctx.beta->B(i, j, k)) *
                                         (ctx.dx * ctx.dz / ctx.dy)));
        if (ctx.perY)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dx * ctx.dz / ctx.dy);
            const double bC = ctx.beta->B(i, j, k);
            const double bn = ctx.beta->B(i, 0, k);
            return hmean(bn, bC) * (ctx.dx * ctx.dz / ctx.dy);
        }
        return useTy ? ctx.trans->TNb(i, k)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dx * ctx.dz / ctx.dy));
    };
    auto Tb = [&](int i, int j, int k)
    {
        if (k > 0)
            return useTz ? ctx.trans->TZ(i, j, k - 1)
                         : (constBeta ? ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz)
                                      : (hmean(ctx.beta->B(i, j, k - 1), ctx.beta->B(i, j, k)) *
                                         (ctx.dx * ctx.dy / ctx.dz)));
        if (ctx.perZ)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz);
            const double bC = ctx.beta->B(i, j, k);
            const double bb = ctx.beta->B(i, j, ctx.nzi - 1);
            return hmean(bb, bC) * (ctx.dx * ctx.dy / ctx.dz);
        }
        return useTz ? ctx.trans->TBb(i, j)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dx * ctx.dy / ctx.dz));
    };
    auto Tt = [&](int i, int j, int k)
    {
        if (k < ctx.nzi - 1)
            return useTz ? ctx.trans->TZ(i, j, k)
                         : (constBeta ? ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz)
                                      : (hmean(ctx.beta->B(i, j, k + 1), ctx.beta->B(i, j, k)) *
                                         (ctx.dx * ctx.dy / ctx.dz)));
        if (ctx.perZ)
        {
            if (constBeta)
                return ctx.beta_const * (ctx.dx * ctx.dy / ctx.dz);
            const double bC = ctx.beta->B(i, j, k);
            const double bt = ctx.beta->B(i, j, 0);
            return hmean(bt, bC) * (ctx.dx * ctx.dy / ctx.dz);
        }
        return useTz ? ctx.trans->TTb(i, j)
                     : ((constBeta ? ctx.beta_const : ctx.beta->B(i, j, k)) *
                        (ctx.dx * ctx.dy / ctx.dz));
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

                // local neighbor indices (ghosts hold periodic wraps)
                const PetscInt iw = i - 1;
                const PetscInt ie = i + 1;
                const PetscInt js = j - 1;
                const PetscInt jn = j + 1;
                const PetscInt kb = k - 1;
                const PetscInt kt = k + 1;

                const double te = hasE ? Te(i, j, k) : 0.0;
                const double tw = hasW ? Tw(i, j, k) : 0.0;
                const double tn = hasN ? Tn(i, j, k) : 0.0;
                const double ts = hasS ? Ts(i, j, k) : 0.0;
                const double tt = hasT ? Tt(i, j, k) : 0.0;
                const double tb = hasB ? Tb(i, j, k) : 0.0;

                // central diagonal and 6 neighbors (Dirichlet faces add to diagonal;
                // Neumann faces affect only RHS -> not part of the matrix)
                double diag = 0.0;

                MatStencil row, cols[7];
                PetscScalar vals[7];
                PetscInt nset = 0;
                // Use global logical stencil indices (DMDA handles periodic wraps).
                row.k = k;
                row.j = j;
                row.i = i;
                row.c = 0;

                if (hasW)
                {
                    cols[nset] = MatStencil{.k = k, .j = j, .i = iw, .c = 0};
                    vals[nset] = -invV * tw;
                    diag += tw;
                    ++nset;
                }
                else if (!ctx.perX && ctx.pbc.W &&
                         norm_p_type(ctx.pbc.W->type) == BcSpec::Type::dirichlet)
                    diag += Tw(i, j, k);

                if (hasE)
                {
                    cols[nset] = MatStencil{.k = k, .j = j, .i = ie, .c = 0};
                    vals[nset] = -invV * te;
                    diag += te;
                    ++nset;
                }
                else if (!ctx.perX && ctx.pbc.E &&
                         norm_p_type(ctx.pbc.E->type) == BcSpec::Type::dirichlet)
                    diag += Te(i, j, k);

                if (hasS)
                {
                    cols[nset] = MatStencil{.k = k, .j = js, .i = i, .c = 0};
                    vals[nset] = -invV * ts;
                    diag += ts;
                    ++nset;
                }
                else if (!ctx.perY && ctx.pbc.S &&
                         norm_p_type(ctx.pbc.S->type) == BcSpec::Type::dirichlet)
                    diag += Ts(i, j, k);

                if (hasN)
                {
                    cols[nset] = MatStencil{.k = k, .j = jn, .i = i, .c = 0};
                    vals[nset] = -invV * tn;
                    diag += tn;
                    ++nset;
                }
                else if (!ctx.perY && ctx.pbc.N &&
                         norm_p_type(ctx.pbc.N->type) == BcSpec::Type::dirichlet)
                    diag += Tn(i, j, k);

                if (hasB)
                {
                    cols[nset] = MatStencil{.k = kb, .j = j, .i = i, .c = 0};
                    vals[nset] = -invV * tb;
                    diag += tb;
                    ++nset;
                }
                else if (!ctx.perZ && ctx.pbc.B &&
                         norm_p_type(ctx.pbc.B->type) == BcSpec::Type::dirichlet)
                    diag += Tb(i, j, k);

                if (hasT)
                {
                    cols[nset] = MatStencil{.k = kt, .j = j, .i = i, .c = 0};
                    vals[nset] = -invV * tt;
                    diag += tt;
                    ++nset;
                }
                else if (!ctx.perZ && ctx.pbc.T &&
                         norm_p_type(ctx.pbc.T->type) == BcSpec::Type::dirichlet)
                    diag += Tt(i, j, k);

                // center
                cols[nset] = row;
                vals[nset] = invV * diag;
                ++nset;

                PetscCall(MatSetValuesStencil(Afill, 1, &row, nset, cols, vals, INSERT_VALUES));
            }

    PetscCall(MatAssemblyBegin(Afill, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Afill, MAT_FINAL_ASSEMBLY));
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
    int nxi_glob{}, nyi_glob{}, nzi_glob{};

    // Borrowed communicator owned by the main application.
    // Lifetime is managed externally; we must not duplicate or free it.
    MPI_Comm comm{MPI_COMM_NULL};
    bool comm_owned{false}; // true only if we create a child/reordered communicator

    // PETSc objects
    DM da_fine{};           // finest DM (interior)
    std::vector<DM> da_lvl; // 0=coarsest ... L-1=finest
    Mat A_shell = nullptr;  // finest-level operator (MatShell)
    std::vector<std::unique_ptr<ShellCtx>> ctx_lvl;
    KSP ksp{};       // outer CG (MG preconditioner inside)
    Vec x{}, b{};    // fine-level unknown/RHS
    Mat Afine_aij{}; // assembled AIJ on finest level (PMAT for MG + outer KSP)
    // Reusable nullspace for pure Neumann problems (created once per hierarchy)
    MatNullSpace ns_const{};

    // host scratch for p
    std::vector<double> p_host;

    // Finest-level β field (β = 1/ρ) stored without ghosts for the MatShell & PMAT assembly.
    // Populated in execute() when variable density is active.
    LevelBeta beta_fine;

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
        // PCMG retains borrowed refs; destroy our ref
        if (A_shell)
        {
            MatDestroy(&A_shell);
            A_shell = nullptr;
        }
        // cached diagonal on finest ctx if present
        for (auto& c : ctx_lvl)
        {
            if (c && c->diag)
            {
                VecDestroy(&c->diag);
                c->diag = nullptr;
            }
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
        // Only free if we created a child communicator; never free the caller's.
        if (comm_owned && comm != MPI_COMM_NULL)
        {
            MPI_Comm_free(&comm);
        }
        comm_owned = false;
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

    // allocate boundary slabs (physical boundaries use these; periodic axes can ignore them)
    T.Tw_bnd.assign(std::size_t(B.nyi) * B.nzi, 0.0);
    T.Te_bnd.assign(std::size_t(B.nyi) * B.nzi, 0.0);
    T.Ts_bnd.assign(std::size_t(B.nxi) * B.nzi, 0.0);
    T.Tn_bnd.assign(std::size_t(B.nxi) * B.nzi, 0.0);
    T.Tb_bnd.assign(std::size_t(B.nxi) * B.nyi, 0.0);
    T.Tt_bnd.assign(std::size_t(B.nxi) * B.nyi, 0.0);

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

    // Fine-level boundary face conductances (one-sided: use local cell β)
    // WEST/EAST
    for (int k = 0; k < B.nzi; ++k)
        for (int j = 0; j < B.nyi; ++j)
        {
            T.Tw_bnd[T.idxWE(j, k)] = B.B(0, j, k) * (dy * dz / dx);
            T.Te_bnd[T.idxWE(j, k)] = B.B(B.nxi - 1, j, k) * (dy * dz / dx);
        }
    // SOUTH/NORTH
    for (int k = 0; k < B.nzi; ++k)
        for (int i = 0; i < B.nxi; ++i)
        {
            T.Ts_bnd[T.idxSN(i, k)] = B.B(i, 0, k) * (dx * dz / dy);
            T.Tn_bnd[T.idxSN(i, k)] = B.B(i, B.nyi - 1, k) * (dx * dz / dy);
        }
    // BOTTOM/TOP
    for (int j = 0; j < B.nyi; ++j)
        for (int i = 0; i < B.nxi; ++i)
        {
            T.Tb_bnd[T.idxBT(i, j)] = B.B(i, j, 0) * (dx * dy / dz);
            T.Tt_bnd[T.idxBT(i, j)] = B.B(i, j, B.nzi - 1) * (dx * dy / dz);
        }
}

// Conservative face aggregation for semi-coarsening:
// (rx,ry,rz) ∈ {1,2}^3 ; sum only across axes actually coarsened; copy-through when factor==1.
static void coarsen_trans(const LevelTrans& Tf, LevelTrans& Tc, int rx, int ry, int rz, bool perX,
                          bool perY, bool perZ)
{
    Tc.Tx.assign(std::size_t(std::max(0, Tc.nxi - 1)) * Tc.nyi * Tc.nzi, 0.0);
    Tc.Ty.assign(std::size_t(Tc.nxi) * std::max(0, Tc.nyi - 1) * Tc.nzi, 0.0);
    Tc.Tz.assign(std::size_t(Tc.nxi) * Tc.nyi * std::max(0, Tc.nzi - 1), 0.0);
    // allocate boundary slabs on coarse (sizes depend on coarse extents)
    Tc.Tw_bnd.assign(std::size_t(Tc.nyi) * Tc.nzi, 0.0);
    Tc.Te_bnd.assign(std::size_t(Tc.nyi) * Tc.nzi, 0.0);
    Tc.Ts_bnd.assign(std::size_t(Tc.nxi) * Tc.nzi, 0.0);
    Tc.Tn_bnd.assign(std::size_t(Tc.nxi) * Tc.nzi, 0.0);
    Tc.Tb_bnd.assign(std::size_t(Tc.nxi) * Tc.nyi, 0.0);
    Tc.Tt_bnd.assign(std::size_t(Tc.nxi) * Tc.nyi, 0.0);

    // X faces
    for (int K = 0; K < Tc.nzi; ++K)
        for (int J = 0; J < Tc.nyi; ++J)
            for (int I = 0; I < Tc.nxi - 1; ++I)
            {
                const int iF = (rx == 2) ? (2 * I + 1) : I;
                const int j0 = (ry == 2) ? (2 * J) : J;
                const int k0 = (rz == 2) ? (2 * K) : K;
                const int jc = (ry == 2) ? 2 : 1, kc = (rz == 2) ? 2 : 1;
                double s = 0.0;
                for (int dj = 0; dj < jc; ++dj)
                    for (int dk = 0; dk < kc; ++dk)
                    {
                        const int jF = j0 + dj, kF = k0 + dk;
                        if (iF < Tf.nxi - 1 && jF < Tf.nyi && kF < Tf.nzi)
                            s += Tf.TX(iF, jF, kF);
                    }
                Tc.Tx[Tc.idxTx(I, J, K)] = s;
            }
    // Y faces
    for (int K = 0; K < Tc.nzi; ++K)
        for (int J = 0; J < Tc.nyi - 1; ++J)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int jF = (ry == 2) ? (2 * J + 1) : J;
                const int i0 = (rx == 2) ? (2 * I) : I;
                const int k0 = (rz == 2) ? (2 * K) : K;
                const int ic = (rx == 2) ? 2 : 1, kc = (rz == 2) ? 2 : 1;
                double s = 0.0;
                for (int di = 0; di < ic; ++di)
                    for (int dk = 0; dk < kc; ++dk)
                    {
                        const int iF = i0 + di, kF = k0 + dk;
                        if (jF < Tf.nyi - 1 && iF < Tf.nxi && kF < Tf.nzi)
                            s += Tf.TY(iF, jF, kF);
                    }
                Tc.Ty[Tc.idxTy(I, J, K)] = s;
            }
    // Z faces
    for (int K = 0; K < Tc.nzi - 1; ++K)
        for (int J = 0; J < Tc.nyi; ++J)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int kF = (rz == 2) ? (2 * K + 1) : K;
                const int i0 = (rx == 2) ? (2 * I) : I;
                const int j0 = (ry == 2) ? (2 * J) : J;
                const int ic = (rx == 2) ? 2 : 1, jc = (ry == 2) ? 2 : 1;
                double s = 0.0;
                for (int di = 0; di < ic; ++di)
                    for (int dj = 0; dj < jc; ++dj)
                    {
                        const int iF = i0 + di, jF = j0 + dj;
                        if (kF < Tf.nzi - 1 && iF < Tf.nxi && jF < Tf.nyi)
                            s += Tf.TZ(iF, jF, kF);
                    }
                Tc.Tz[Tc.idxTz(I, J, K)] = s;
            }

    // --- Boundary slabs (only meaningful on non-periodic axes) ---
    // WEST/EAST: aggregate across coarsened Y/Z
    if (!perX)
    {
        for (int K = 0; K < Tc.nzi; ++K)
            for (int J = 0; J < Tc.nyi; ++J)
            {
                const int j0 = (ry == 2) ? 2 * J : J;
                const int k0 = (rz == 2) ? 2 * K : K;
                const int jc = (ry == 2) ? 2 : 1;
                const int kc = (rz == 2) ? 2 : 1;
                double sw = 0.0, se = 0.0;
                for (int dj = 0; dj < jc; ++dj)
                    for (int dk = 0; dk < kc; ++dk)
                    {
                        const int jF = j0 + dj, kF = k0 + dk;
                        if (jF < Tf.nyi && kF < Tf.nzi)
                        {
                            sw += Tf.TWb(jF, kF);
                            se += Tf.TEb(jF, kF);
                        }
                    }
                Tc.Tw_bnd[Tc.idxWE(J, K)] = sw;
                Tc.Te_bnd[Tc.idxWE(J, K)] = se;
            }
    }
    // SOUTH/NORTH: aggregate across coarsened X/Z
    if (!perY)
    {
        for (int K = 0; K < Tc.nzi; ++K)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int i0 = (rx == 2) ? 2 * I : I;
                const int k0 = (rz == 2) ? 2 * K : K;
                const int ic = (rx == 2) ? 2 : 1;
                const int kc = (rz == 2) ? 2 : 1;
                double ss = 0.0, sn = 0.0;
                for (int di = 0; di < ic; ++di)
                    for (int dk = 0; dk < kc; ++dk)
                    {
                        const int iF = i0 + di, kF = k0 + dk;
                        if (iF < Tf.nxi && kF < Tf.nzi)
                        {
                            ss += Tf.TSb(iF, kF);
                            sn += Tf.TNb(iF, kF);
                        }
                    }
                Tc.Ts_bnd[Tc.idxSN(I, K)] = ss;
                Tc.Tn_bnd[Tc.idxSN(I, K)] = sn;
            }
    }
    // BOTTOM/TOP: aggregate across coarsened X/Y
    if (!perZ)
    {
        for (int J = 0; J < Tc.nyi; ++J)
            for (int I = 0; I < Tc.nxi; ++I)
            {
                const int i0 = (rx == 2) ? 2 * I : I;
                const int j0 = (ry == 2) ? 2 * J : J;
                const int ic = (rx == 2) ? 2 : 1;
                const int jc = (ry == 2) ? 2 : 1;
                double sb = 0.0, st = 0.0;
                for (int di = 0; di < ic; ++di)
                    for (int dj = 0; dj < jc; ++dj)
                    {
                        const int iF = i0 + di, jF = j0 + dj;
                        if (iF < Tf.nxi && jF < Tf.nyi)
                        {
                            sb += Tf.TBb(iF, jF);
                            st += Tf.TTb(iF, jF);
                        }
                    }
                Tc.Tb_bnd[Tc.idxBT(I, J)] = sb;
                Tc.Tt_bnd[Tc.idxBT(I, J)] = st;
            }
    }
}

// Build solver hierarchy. Periodicity comes from both the BC table and mesh.periodic
// meshPerX/Y/Z override acts as a single source of truth when the mesh is periodic
// Pass the desired process grid (m,n,p). If any are 0, we'll auto-compute from MPI size
static void build_hierarchy(PPImpl& impl, MPI_Comm user_comm_in, const PBC& pbc, bool meshPerX,
                            bool meshPerY, bool meshPerZ, std::array<int, 3> proc_grid)
{

    // Start from the application's communicator (borrowed).
    if (impl.comm == MPI_COMM_NULL)
        impl.comm = user_comm_in;
    MPI_Comm comm = impl.comm;

    {
        int sz = -1, rk = -1;
        if (MPI_Comm_size(comm, &sz) != MPI_SUCCESS || sz < 1)
        {
            SETERRABORT(comm, PETSC_ERR_ARG_INCOMP,
                        "invalid MPI_Comm passed into PressurePoisson (size=%d) — must be a live "
                        "communicator",
                        sz);
        }
        MPI_Comm_rank(comm, &rk); // warm up the comm
    }

    // --- If the app's communicator is Cartesian, derive a *reordered child communicator*
    //     whose rank order is lexicographic in (x,y,z). This matches DMDA's process order.
    //     We keep the same process set (color=0), only the ordering (key) changes.
    {
        int topo_type = MPI_UNDEFINED;
        MPI_Topo_test(comm, &topo_type);
        if (topo_type == MPI_CART)
        {
            int ndims = 0;
            MPI_Cartdim_get(comm, &ndims);
            if (ndims == 3)
            {
                int dims[3] = {0, 0, 0}, periods[3] = {0, 0, 0}, coords[3] = {0, 0, 0};
                MPI_Cart_get(comm, 3, dims, periods, coords);
                // Lexicographic key: z-major over (x,y,z) -> (k * Ny * Nx) + j * Nx + i
                const int key = coords[2] * (dims[0] * dims[1]) + coords[1] * dims[0] + coords[0];
                MPI_Comm child = MPI_COMM_NULL;
                // color=0 keeps the same group; key sets rank order inside it.
                MPI_Comm_split(comm, /*color*/ 0, /*key*/ key, &child);
                if (child != MPI_COMM_NULL)
                {
                    // If we previously owned a child, free it before replacing.
                    if (impl.comm_owned && impl.comm != MPI_COMM_NULL && impl.comm != user_comm_in)
                    {
                        MPI_Comm_free(&impl.comm);
                    }
                    impl.comm = child;
                    impl.comm_owned = true;
                    comm = impl.comm;
                }
            }
        }
    }

    // ----- finest DMDA over interior cells -----
    const PetscInt nxi = impl.nxi_glob; // use mesh global interior sizes
    const PetscInt nyi = impl.nyi_glob;
    const PetscInt nzi = impl.nzi_glob;

    // Periodicity from BCs OR mesh flags (mesh flags take precedence)
    auto is_periodic = [](const BcSpec* s) { return s && s->type == BcSpec::Type::periodic; };
    const bool perX = meshPerX || (is_periodic(pbc.W) && is_periodic(pbc.E));
    const bool perY = meshPerY || (is_periodic(pbc.S) && is_periodic(pbc.N));
    const bool perZ = meshPerZ || (is_periodic(pbc.B) && is_periodic(pbc.T));
    DMBoundaryType bx = perX ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType by = perY ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType bz = perZ ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;

    // ---------------- Processor grid ----------------
    // If the communicator is (still) Cartesian, read its dims; otherwise fall back
    // to MPI_Dims_create using the mesh hint. With the child comm (if created),
    // the rank order now matches DMDA's lexicographic assumptions.
    int dims3[3] = {proc_grid[0], proc_grid[1], proc_grid[2]};
    {
        int topo_type = MPI_UNDEFINED;
        MPI_Topo_test(comm, &topo_type);
        if (topo_type == MPI_CART)
        {
            int cdims[3] = {0, 0, 0};
            int cper[3] = {0, 0, 0};
            int coords_dummy[3] = {0, 0, 0};
            int ndims = 0;
            MPI_Cartdim_get(comm, &ndims);
            if (ndims == 3)
            {
                MPI_Cart_get(comm, 3, cdims, cper, coords_dummy);
                // Keep (x,y,z) order consistent with MPI_Cart dims and PETSc DMDA:
                dims3[0] = cdims[0]; // x
                dims3[1] = cdims[1]; // y
                dims3[2] = cdims[2]; // z
            }
        }
        else
        {
            int size_world = 1;
            MPI_Comm_size(comm, &size_world);
            bool user_fixed = (dims3[0] > 0 && dims3[1] > 0 && dims3[2] > 0);
            if (user_fixed)
            {
                const long long prod = 1LL * dims3[0] * dims3[1] * dims3[2];
                if (prod != size_world)
                {
                    SETERRABORT(comm, PETSC_ERR_ARG_INCOMP,
                                "mesh.proc_grid=%dx%dx%d must multiply to communicator size %d",
                                dims3[0], dims3[1], dims3[2], size_world);
                }
            }
            else
            {
                int tmp[3] = {dims3[0], dims3[1], dims3[2]};
                int sz;
                MPI_Comm_size(comm, &sz);
                int prod = 1;
                for (int a = 0; a < 3; ++a)
                    if (tmp[a] > 0)
                        prod *= tmp[a];
                if (sz % prod != 0)
                {
                    SETERRABORT(comm, PETSC_ERR_ARG_INCOMP,
                                "mesh.proc_grid contains fixed axes that do not divide MPI size");
                }
                int need = sz / prod;
                int fill[3] = {0, 0, 0};
                MPI_Dims_create(need, 3, fill);
                for (int a = 0, t = 0; a < 3; ++a)
                    if (tmp[a] == 0)
                        tmp[a] = (fill[t++] ? fill[t - 1] : 1);
                dims3[0] = tmp[0];
                dims3[1] = tmp[1];
                dims3[2] = tmp[2];
            }
        }
    }

    // Compute explicit ownership ranges so DMDA matches our mesh partition exactly.
    auto make_ownership = [](PetscInt n, int p)
    {
        std::vector<PetscInt> l(std::max(1, p), 0);
        PetscInt q = n / p, r = n % p;
        for (int i = 0; i < p; ++i)
            l[i] = q + (i < r ? 1 : 0);
        return l;
    };
    const int px = dims3[0], py = dims3[1], pz = dims3[2];
    std::vector<PetscInt> lx = make_ownership(nxi, px);
    std::vector<PetscInt> ly = make_ownership(nyi, py);
    std::vector<PetscInt> lz = make_ownership(nzi, pz);

    // Create the finest DMDA with an explicit (px,py,pz) that multiplies to size.
    DMDACreate3d(comm, bx, by, bz, DMDA_STENCIL_STAR, nxi, nyi, nzi, px, py, pz,
                 /*dof*/ 1, /*stencil width*/ 1,
                 /*lx*/ lx.data(), /*ly*/ ly.data(), /*lz*/ lz.data(), &impl.da_fine);
    // Allow -da_* options at runtime for debugging/overrides
    DMSetFromOptions(impl.da_fine);
    DMSetUp(impl.da_fine);

    // Sanity: the process grid must not exceed the per-axis cell counts.
    // (Catches impossible user hints early and avoids hard-to-trace segfaults later.)
    if (!(dims3[0] <= (int) nxi && dims3[1] <= (int) nyi && dims3[2] <= (int) nzi))
    {
        SETERRABORT(comm, PETSC_ERR_ARG_OUTOFRANGE, "proc grid %dx%dx%d exceeds cell grid %dx%dx%d",
                    dims3[0], dims3[1], dims3[2], (int) nxi, (int) nyi, (int) nzi);
    }

    DMDASetInterpolationType(impl.da_fine, DMDA_Q0);

    // --- Build hierarchy using PETSc's DMCoarsen, which naturally handles odd sizes. ---

    // IMPORTANT for DMDA_Q0 (cell-centered) interpolation used by PCMG:
    //   DMCreateInterpolation(DMDA_Q0) requires that the fine-grid points be an
    //   integer multiple of the coarse-grid points along each axis. When halving
    //   an odd size (e.g., 129 -> 65), that multiple condition is violated
    //   (129 % 65 != 0) and PCSetUp/DMCreateInterpolation will error.
    //   We therefore refuse to coarsen along any axis whose size is odd.
    //
    //   Result: fewer MG levels on odd-sized axes, but robust setup without
    //   "Fine grid points must be multiple of coarse grid points" failures.
    auto can_coarsen_q0_once = [](PetscInt a) -> bool
    {
        // Q0 (cell-centered) interpolation needs fine%coarse==0 at each pair
        // AND PETSc requires the *coarse* grid to have at least 2 points.
        // With DMDA coarsening a2 = ceil(a/2):
        //  - disallow odd a (fine%coarse != 0)
        //  - disallow the last step 2 -> 1 (coarse<2)
        if (a % 2 != 0)
            return false;                 // must be even
        const PetscInt a2 = (a + 1) >> 1; // proposed coarse size
        return a2 >= 2;                   // keep coarse >= 2
    };
    auto next_size = [](PetscInt a) -> PetscInt
    {
        // PETSc DMDA coarsening halves with ceil
        return (a + 1) >> 1;
    };

    impl.da_lvl.clear();
    {
        // Build hierarchy but stop BEFORE a coarse level would violate proc-per-axis counts.
        // If a prospective (M2,N2,P2) would have M2 < px (or Y/Z analogues), stop coarsening.
        PetscInt M = 0, N = 0, P = 0, px = 0, py = 0, pz = 0;
        DMDAGetInfo(impl.da_fine, nullptr, &M, &N, &P, &px, &py, &pz, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr);

        auto can_coarsen_once = [&](PetscInt a, PetscInt pa) -> bool
        {
            // PETSc DA halves with ceil; require:
            //  - Q0 compatibility AND coarse >= 2 (handled by can_coarsen_q0_once)
            //  - next >= procs on that axis (layout validity)
            const PetscInt a2 = (a + 1) >> 1;
            return can_coarsen_q0_once(a) && (a2 >= pa);
        };

        DM cur = impl.da_fine;
        PetscObjectReference((PetscObject) cur);
        std::vector<DM> chain_f2c;
        chain_f2c.push_back(cur);
        while (true)
        {
            // Check prospective sizes *before* asking PETSc to coarsen
            PetscInt Mc = 0, Nc = 0, Pc = 0, px_c = 0, py_c = 0, pz_c = 0;
            DMDAGetInfo(cur, nullptr, &Mc, &Nc, &Pc, &px_c, &py_c, &pz_c, nullptr, nullptr, nullptr,
                        nullptr, nullptr, nullptr);
            // Additional DMDA_Q0 + coarse>=2 constraint per axis.
            const bool q0_ok_x = can_coarsen_q0_once(Mc);
            const bool q0_ok_y = can_coarsen_q0_once(Nc);
            const bool q0_ok_z = can_coarsen_q0_once(Pc);

            // Also maintain processor-layout validity on the next level
            const bool layout_ok_x = can_coarsen_once(Mc, px_c);
            const bool layout_ok_y = can_coarsen_once(Nc, py_c);
            const bool layout_ok_z = can_coarsen_once(Pc, pz_c);

            if (!(q0_ok_x && q0_ok_y && q0_ok_z && layout_ok_x && layout_ok_y && layout_ok_z))
            {
                break; // next level would be invalid for current proc layout
            }
            DM next = NULL;
            PetscErrorCode ierr = DMCoarsen(cur, MPI_COMM_NULL, &next);
            if (ierr || !next)
                break;
            chain_f2c.push_back(next);
            cur = next;
        }
        // PCMG expects 0=coarsest ... L-1=finest
        impl.da_lvl.assign(chain_f2c.rbegin(), chain_f2c.rend());
    }
    const PetscInt L = (PetscInt) impl.da_lvl.size();

    // ----- nullspace applicability (robust) -----
    //
    // A pressure Poisson with *no* Dirichlet on any *physical* boundary is singular:
    // the constant vector is a nullspace (gauge freedom). Periodic faces do not fix
    // the gauge. So we mark the operator as singular iff there is no Dirichlet face
    // on any axis that is NOT periodic. Otherwise (any Dirichlet present) it is SPD.
    //
    auto is_dir = [](const BcSpec* s) -> bool
    { return s && norm_p_type(s->type) == BcSpec::Type::dirichlet; };

    // Dirichlet may contribute only on non-periodic axes
    const bool hasDirX = (!perX) && (is_dir(pbc.W) || is_dir(pbc.E));
    const bool hasDirY = (!perY) && (is_dir(pbc.S) || is_dir(pbc.N));
    const bool hasDirZ = (!perZ) && (is_dir(pbc.B) || is_dir(pbc.T));
    const bool hasAnyDir = hasDirX || hasDirY || hasDirZ;

    // Singular iff no Dirichlet on physical boundary (includes fully periodic boxes).
    impl.all_neumann = !hasAnyDir;

    // ----- build MatShell operators and KSP/PCMG -----
    impl.ctx_lvl.clear();
    impl.ctx_lvl.resize(L); // only [L-1] will be populated

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
        KSPSetPCSide(impl.ksp, PC_LEFT);
        KSPSetNormType(impl.ksp, KSP_NORM_PRECONDITIONED); // Approx control of ||r||_2
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

    // Only create finest-level MatShell/ctx (Amat). Coarser levels are handled by PMAT+Galerkin.
    {
        const PetscInt lf = L - 1;
        PetscInt nxi_l = 0, nyi_l = 0, nzi_l = 0;
        DMDAGetInfo(impl.da_lvl[lf], nullptr, &nxi_l, &nyi_l, &nzi_l, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        // Create a temporary DMDA vector to read both GLOBAL and LOCAL sizes.
        Vec vtmp;
        DMCreateGlobalVector(impl.da_lvl[lf], &vtmp);
        PetscInt M = 0, N = 0, mloc = 0;
        VecGetSize(vtmp, &M);
        VecGetLocalSize(vtmp, &mloc);
        N = M;
        VecDestroy(&vtmp);

        auto ctx = std::make_unique<ShellCtx>();
        ctx->da = impl.da_lvl[lf];
        ctx->nxi = (int) nxi_l;
        ctx->nyi = (int) nyi_l;
        ctx->nzi = (int) nzi_l;
        ctx->ng = impl.ng;
        ctx->nxc_tot = impl.nxc_tot;
        ctx->nyc_tot = impl.nyc_tot;
        ctx->nzc_tot = impl.nzc_tot;
        ctx->dx = impl.dx;
        ctx->dy = impl.dy;
        ctx->dz = impl.dz;
        ctx->trans = nullptr; // per-level Tx/Ty/Tz dropped (computed on the fly/fine)
        ctx->beta = nullptr;  // see execute(): finest PMAT is (re)assembled as needed
        ctx->use_const_beta = !impl.varrho;
        ctx->beta_const = impl.varrho ? 1.0 : (1.0 / impl.rho_const);
        ctx->pbc = pbc;
        ctx->perX = perX;
        ctx->perY = perY;
        ctx->perZ = perZ;

        Mat A;
        // IMPORTANT: pass LOCAL sizes so the MatShell distribution exactly matches
        // the DMDA vector distribution used by KSP/MatMult.
        MatCreateShell(comm,
                       /*m(local rows)*/ mloc,
                       /*n(local cols)*/ mloc,
                       /*M(global rows)*/ M,
                       /*N(global cols)*/ N, (void*) ctx.get(), &A);
        MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) ShellMult);
        MatShellSetOperation(A, MATOP_GET_DIAGONAL, (void (*)(void)) ShellGetDiagonal);
        MatShellSetOperation(A, MATOP_MULT_TRANSPOSE, (void (*)(void)) ShellMult);
        MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
        if (!impl.all_neumann)
            MatSetOption(A, MAT_SPD, PETSC_TRUE);

        impl.ctx_lvl[lf] = std::move(ctx);
        impl.A_shell = A;
    }

    // MG cycle config
    // For all-Neumann (singular) prefer additive MG to keep the preconditioner benign/symmetric for
    // MINRES; otherwise multiplicative V is fine.
    if (impl.all_neumann)
    {
        PCMGSetType(pc, PC_MG_ADDITIVE);
    }
    else
    {
        PCMGSetType(pc, PC_MG_MULTIPLICATIVE);
    }
    PCMGSetNumberSmooth(pc, 2);

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

            // Use Chebyshev + Jacobi smoothing on all levels.
            PetscCallAbort(comm, KSPSetType(kspl, KSPCHEBYSHEV));
            PetscCallAbort(comm, KSPGetPC(kspl, &pcl));
            PetscCallAbort(comm, PCSetType(pcl, PCJACOBI));
            PCJacobiSetType(pcl, PC_JACOBI_DIAGONAL);
            // Auto-estimate spectral bounds (cheap power iterations), robust across levels.
            PetscCallAbort(comm, KSPChebyshevEstEigSet(kspl, PETSC_DEFAULT, PETSC_DEFAULT,
                                                       PETSC_DEFAULT, PETSC_DEFAULT));
            PetscCallAbort(comm, KSPSetConvergenceTest(kspl, KSPConvergedSkip, NULL, NULL));
            PetscCallAbort(comm, KSPSetNormType(kspl, KSP_NORM_NONE));
            PetscCallAbort(comm,
                           KSPSetTolerances(kspl, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 2));
        }
    }

    // ---- Coarse level ----
    // SPD: preonly + Cholesky (fast).  Pure-Neumann (singular): preonly + SVD.
    {
        KSP kspc = NULL;
        PC pcc = NULL;
        PetscCallAbort(comm, PCMGGetCoarseSolve(pc, &kspc));
        if (kspc)
        {
            PetscCallAbort(comm, KSPSetType(kspc, KSPPREONLY));
            PetscCallAbort(comm, KSPGetPC(kspc, &pcc));
            if (impl.all_neumann)
            {
                // Robust for rank-deficient coarse operators
                PetscCallAbort(comm, PCSetType(pcc, PCSVD));
            }
            else
            {
                PetscCallAbort(comm, PCSetType(pcc, PCCHOLESKY));
            }
        }
    }

    // fine vectors
    DMCreateGlobalVector(impl.da_fine, &impl.x);
    DMCreateGlobalVector(impl.da_fine, &impl.b);

    // Assemble fine-level AIJ (PMAT) from the fine ShellCtx and wire operators:
    //  - Outer KSP: Amat = fine MatShell; Pmat = fine AIJ
    //  - PCMG: set only the finest level operators; Galerkin builds coarse PMATs
    // Create on first build; later we will reuse the same matrix pattern.
    if (!impl.Afine_aij)
        impl.Afine_aij = nullptr;
    PetscCallAbort(comm, AssembleAIJFromShell(*impl.ctx_lvl[L - 1], &impl.Afine_aij));
    MatSetOption(impl.Afine_aij, MAT_SYMMETRIC, PETSC_TRUE);
    if (!impl.all_neumann)
    {
        MatSetOption(impl.Afine_aij, MAT_SPD, PETSC_TRUE);
    }

    // Outer KSP operators
    KSPSetOperators(impl.ksp, impl.A_shell, impl.Afine_aij);

    // Finest-level MG operators; coarser levels via Galerkin(PMAT)
    PC pc_for_ops = NULL;
    KSPGetPC(impl.ksp, &pc_for_ops);
    PCMGSetOperators(pc_for_ops, L - 1, impl.A_shell, impl.Afine_aij);

    retune_chebyshev_on_levels(impl.ksp);

    // ---------- Proper nullspace propagation (pure Neumann) ----------
    if (impl.all_neumann)
    {
        // One constant nullspace object reused everywhere
        if (!impl.ns_const)
        {
            PetscCallAbort(comm, MatNullSpaceCreate(comm, PETSC_TRUE, 0, NULL, &impl.ns_const));
        }
        // Attach to fine Amat and Pmat so KSP can handle the singular operator.
        PetscCallAbort(comm, MatSetNullSpace(impl.A_shell, impl.ns_const));
        PetscCallAbort(comm, MatSetNullSpace(impl.Afine_aij, impl.ns_const));
        // Help MG/AMG via near-nullspace on PMAT; coarse PMATs built by Galerkin will inherit.
        PetscCallAbort(comm, MatSetNearNullSpace(impl.Afine_aij, impl.ns_const));
    }

    KSPSetFromOptions(impl.ksp); // allow CLI/options to refine everything

    // ---- Lightweight KSP monitor into our logger (DEBUG only) ----
    // Logs: [poisson] it=… ||r||=…   only when SOLVER_LOG>=debug.
    {
        auto monitor = [](KSP ksp, PetscInt it, PetscReal rnorm, void*) -> PetscErrorCode
        {
            using core::master::logx::Level;
            if (core::master::logx::g_level.load() >= Level::Debug)
            {
                LOGD("[poisson] it=%d  ||r||=%.6e\n", (int) it, (double) rnorm);
            }
            return 0;
        };
        // Attach after KSP is configured so options may still add PETSc's own monitors.
        // (Multiple monitors can coexist.)
        PetscErrorCode ierr = KSPMonitorSet(impl.ksp, monitor, nullptr, nullptr);
        if (ierr)
        {
            LOGW("[poisson] unable to attach KSP monitor (ierr=%d)\n", (int) ierr);
        }
    }
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

    // Single source of truth for halo width
    const int ng = tile.mesh ? tile.mesh->ng : 0;
    // Optional sanity: extents must be consistent with mesh halo
    // (interior along x is nxc_tot - 2*ng)
    const int nx_interior = nxc_tot - 2 * ng;
    (void) nx_interior; // suppress unused warnings in release
    assert(nx_interior >= 1 && "Mesh halo inconsistent with pressure extent totals");

    // Translate DMDA global indices (xs,ys,zs) to this tile's local array offsets.
    // These are the global starting cell-indices (center-based) for this tile.
    const int i0 = (tile.mesh ? tile.mesh->global_lo[0] : 0);
    const int j0 = (tile.mesh ? tile.mesh->global_lo[1] : 0);
    const int k0 = (tile.mesh ? tile.mesh->global_lo[2] : 0);

    // Borrowed communicator: if null, fall back to self for single-rank use.
    // Robust unbox + validate. Prefer the app comm; fall back to PETSC_COMM_WORLD.
    MPI_Comm user_comm = PETSC_COMM_WORLD;
    if (mpi_comm_)
    {
        // Our public headers promise "boxed" MPI_Comm via mpi_box(); unbox it if so.
        MPI_Comm cand = mpi_unbox(mpi_comm_);
        user_comm = valid_comm_or(cand, PETSC_COMM_WORLD);
    }

    // rebuild hierarchy if geometry changed
    if (!impl_ || impl_->nxc_tot != nxc_tot || impl_->nyc_tot != nyc_tot ||
        impl_->nzc_tot != nzc_tot || impl_->ng != ng || // compare against mesh halo
        impl_->dx != dx_ || impl_->dy != dy_ || impl_->dz != dz_)
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
        impl_->nxi_glob = tile.mesh->global[0];
        impl_->nyi_glob = tile.mesh->global[1];
        impl_->nzi_glob = tile.mesh->global[2];
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

        // Read desired process grid from the mesh (if provided), else auto-compute.
        std::array<int, 3> proc_grid = {0, 0, 0};
        if (tile.mesh)
            proc_grid = tile.mesh->proc_grid;
        build_hierarchy(*impl_, user_comm, pbc, meshPerX, meshPerY, meshPerZ, proc_grid);
    }

    // Partition sanity check: with the reordered child communicator, the owned corners
    // should now agree with the application's Cartesian decomposition.
    {
        PetscInt xs_chk = 0, ys_chk = 0, zs_chk = 0, xm_chk = 0, ym_chk = 0, zm_chk = 0;
        DMDAGetCorners(impl_->da_fine, &xs_chk, &ys_chk, &zs_chk, &xm_chk, &ym_chk, &zm_chk);
        const int i0_mesh = (tile.mesh ? tile.mesh->global_lo[0] : 0);
        const int j0_mesh = (tile.mesh ? tile.mesh->global_lo[1] : 0);
        const int k0_mesh = (tile.mesh ? tile.mesh->global_lo[2] : 0);
        if (xs_chk != i0_mesh || ys_chk != j0_mesh || zs_chk != k0_mesh)
        {
            SETERRABORT(impl_->comm, PETSC_ERR_PLIB,
                        "DMDA corner (%d,%d,%d) != mesh.global_lo (%d,%d,%d) after reordering — "
                        "check proc_grid.",
                        (int) xs_chk, (int) ys_chk, (int) zs_chk, i0_mesh, j0_mesh, k0_mesh);
        }
    }

    // ---------------- RHS = - (1/dt) * div(u*) plus BC contributions ----------------

    {
        core::master::exchange_named_fields(fields, *tile.mesh, mpi_comm_, {"u", "v", "w"});
    }

    std::vector<double> div(std::size_t(nxc_tot) * nyc_tot * nzc_tot, 0.0);
    numerics::kernels::divergence(
        static_cast<const double*>(vu.host_ptr), static_cast<const double*>(vv.host_ptr),
        static_cast<const double*>(vw.host_ptr), vu.extents[0], vu.extents[1], vu.extents[2],
        vv.extents[0], vv.extents[1], vv.extents[2], vw.extents[0], vw.extents[1], vw.extents[2],
        nxc_tot, nyc_tot, nzc_tot, ng, dx_, dy_, dz_, div.data());

    // (New) Variable-density path: build only a fine-level β and attach to finest ctx.
    if (impl_->varrho)
    {
        // Size from finest DMDA
        PetscInt nxi_l = 0, nyi_l = 0, nzi_l = 0;
        DMDAGetInfo(impl_->da_lvl.back(), nullptr, &nxi_l, &nyi_l, &nzi_l, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        // persistent storage in impl_
        impl_->beta_fine.nxi = (int) nxi_l;
        impl_->beta_fine.nyi = (int) nyi_l;
        impl_->beta_fine.nzi = (int) nzi_l;
        impl_->beta_fine.data.assign(std::size_t(nxi_l) * nyi_l * nzi_l, 0.0);
        // Fill only this rank's owned DMDA box using global (xs,ys,zs), translating
        // to the tile-local storage via (i - i0) etc.
        PetscInt xs_b, ys_b, zs_b, xm_b, ym_b, zm_b;
        DMDAGetCorners(impl_->da_lvl.back(), &xs_b, &ys_b, &zs_b, &xm_b, &ym_b, &zm_b);
        for (int k = zs_b; k < zs_b + zm_b; ++k)
            for (int j = ys_b; j < ys_b + ym_b; ++j)
                for (int i = xs_b; i < xs_b + xm_b; ++i)
                {
                    const int ic = (i - i0) + impl_->ng;
                    const int jc = (j - j0) + impl_->ng;
                    const int kc = (k - k0) + impl_->ng;
                    // Guard in case a rank owns a DM slab that this tile does not cover
                    if (ic >= 0 && jc >= 0 && kc >= 0 && ic < nxc_tot && jc < nyc_tot &&
                        kc < nzc_tot)
                    {
                        const std::size_t c =
                            std::size_t(ic) +
                            std::size_t(nxc_tot) *
                                (std::size_t(jc) + std::size_t(nyc_tot) * std::size_t(kc));
                        const std::size_t g =
                            std::size_t(i) +
                            std::size_t(impl_->beta_fine.nxi) *
                                (std::size_t(j) +
                                 std::size_t(impl_->beta_fine.nyi) * std::size_t(k));
                        impl_->beta_fine.data[g] =
                            1.0 / static_cast<const double*>(vrho.host_ptr)[c];
                    }
                }
        // attach to finest MatShell ctx and mark diagonal dirty
        auto* fctx = impl_->ctx_lvl.back().get();
        fctx->beta = &impl_->beta_fine;
        fctx->use_const_beta = false;
        fctx->diag_built = false;

        // Reuse existing PMAT pattern; AssembleAIJFromShell() will MatZeroEntries()
        // and overwrite values in-place if impl_->Afine_aij already exists.
        PetscCallAbort(impl_->comm, AssembleAIJFromShell(*fctx, &impl_->Afine_aij));
        if (!impl_->all_neumann)
            MatSetOption(impl_->Afine_aij, MAT_SPD, PETSC_TRUE);

        // Re-attach nullspace (pure Neumann)
        if (impl_->all_neumann)
        {
            if (!impl_->ns_const)
                PetscCallAbort(impl_->comm, MatNullSpaceCreate(impl_->comm, PETSC_TRUE, 0, NULL,
                                                               &impl_->ns_const));
            PetscCallAbort(impl_->comm, MatSetNullSpace(impl_->A_shell, impl_->ns_const));
            PetscCallAbort(impl_->comm, MatSetNullSpace(impl_->Afine_aij, impl_->ns_const));
            PetscCallAbort(impl_->comm, MatSetNearNullSpace(impl_->Afine_aij, impl_->ns_const));
        }
        // refresh KSP/PCMG finest operators (coarse PMATs rebuilt by Galerkin)
        KSPSetOperators(impl_->ksp, impl_->A_shell, impl_->Afine_aij);
        PC pc_refresh = NULL;
        KSPGetPC(impl_->ksp, &pc_refresh);
        const PetscInt L = (PetscInt) impl_->da_lvl.size();
        PCMGSetOperators(pc_refresh, L - 1, impl_->A_shell, impl_->Afine_aij);

        retune_chebyshev_on_levels(impl_->ksp);
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

        // Local DMDA (global) corner and owned size (no ghosts)
        PetscInt xs, ys, zs, xm, ym, zm;
        DMDAGetCorners(impl_->da_fine, &xs, &ys, &zs, &xm, &ym, &zm);

        PetscScalar*** barr;
        DMDAVecGetArray(impl_->da_fine, impl_->b, &barr);
        const int nxc_tot = vp.extents[0], nyc_tot = vp.extents[1], nzc_tot = vp.extents[2];

        // Face transmissibilities on-the-fly (no stored Tx/Ty/Tz):
        const double invV = 1.0 / (dx_ * dy_ * dz_);
        const ShellCtx* fctx = impl_->ctx_lvl.back().get();
        const bool constB = !impl_->varrho;
        const double Tcx = (constB ? (1.0 / rho_) * (dy_ * dz_ / dx_) : 0.0);
        const double Tcy = (constB ? (1.0 / rho_) * (dx_ * dz_ / dy_) : 0.0);
        const double Tcz = (constB ? (1.0 / rho_) * (dx_ * dy_ / dz_) : 0.0);
        const int nxi = nxc_tot - 2 * ng, nyi = nyc_tot - 2 * ng, nzi = nzc_tot - 2 * ng;
        auto Tw = [&](int i, int j, int k)
        {
            if (constB)
                return Tcx;
            const double bC = fctx->beta->B(i, j, k);
            if (i > 0)
                return hmean(fctx->beta->B(i - 1, j, k), bC) * (dy_ * dz_ / dx_);
            if (fctx->perX)
                return hmean(fctx->beta->B(nxi - 1, j, k), bC) * (dy_ * dz_ / dx_);
            return bC * (dy_ * dz_ / dx_);
        };
        auto Te = [&](int i, int j, int k)
        {
            if (constB)
                return Tcx;
            const double bC = fctx->beta->B(i, j, k);
            if (i < nxi - 1)
                return hmean(fctx->beta->B(i + 1, j, k), bC) * (dy_ * dz_ / dx_);
            if (fctx->perX)
                return hmean(fctx->beta->B(0, j, k), bC) * (dy_ * dz_ / dx_);
            return bC * (dy_ * dz_ / dx_);
        };
        auto Ts = [&](int i, int j, int k)
        {
            if (constB)
                return Tcy;
            const double bC = fctx->beta->B(i, j, k);
            if (j > 0)
                return hmean(fctx->beta->B(i, j - 1, k), bC) * (dx_ * dz_ / dy_);
            if (fctx->perY)
                return hmean(fctx->beta->B(i, nyi - 1, k), bC) * (dx_ * dz_ / dy_);
            return bC * (dx_ * dz_ / dy_);
        };
        auto Tn = [&](int i, int j, int k)
        {
            if (constB)
                return Tcy;
            const double bC = fctx->beta->B(i, j, k);
            if (j < nyi - 1)
                return hmean(fctx->beta->B(i, j + 1, k), bC) * (dx_ * dz_ / dy_);
            if (fctx->perY)
                return hmean(fctx->beta->B(i, 0, k), bC) * (dx_ * dz_ / dy_);
            return bC * (dx_ * dz_ / dy_);
        };
        auto Tb = [&](int i, int j, int k)
        {
            if (constB)
                return Tcz;
            const double bC = fctx->beta->B(i, j, k);
            if (k > 0)
                return hmean(fctx->beta->B(i, j, k - 1), bC) * (dx_ * dy_ / dz_);
            if (fctx->perZ)
                return hmean(fctx->beta->B(i, j, nzi - 1), bC) * (dx_ * dy_ / dz_);
            return bC * (dx_ * dy_ / dz_);
        };
        auto Tt = [&](int i, int j, int k)
        {
            if (constB)
                return Tcz;
            const double bC = fctx->beta->B(i, j, k);
            if (k < nzi - 1)
                return hmean(fctx->beta->B(i, j, k + 1), bC) * (dx_ * dy_ / dz_);
            if (fctx->perZ)
                return hmean(fctx->beta->B(i, j, 0), bC) * (dx_ * dy_ / dz_);
            return bC * (dx_ * dy_ / dz_);
        };

        for (int k = zs; k < zs + zm; ++k)
            for (int j = ys; j < ys + ym; ++j)
                for (int i = xs; i < xs + xm; ++i)
                {
                    const int ii = (int) (i - i0) + ng;
                    const int jj = (int) (j - j0) + ng;
                    const int kk = (int) (k - k0) + ng;
                    const std::size_t c =
                        std::size_t(ii) +
                        std::size_t(nxc_tot) *
                            (std::size_t(jj) + std::size_t(nyc_tot) * std::size_t(kk));
                    // A is -div(β∇p) (SPD). Projection gives -div(β∇p) = -(1/dt) div(u*).
                    double rhs = -div[c] / dt;

                    // Use the finest-level periodic flags
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
                    auto add_dirichlet =
                        [&](const BcSpec* s, bool per_axis, double coef, double value)
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
                        add_dirichlet(pbc.W, perX, Tw(i, j, k), pbc.W ? pbc.W->value : 0.0);
                        add_neumann(pbc.W, Tw(i, j, k), dx_, -1);
                    }
                    if (!hasE && !perX)
                    {
                        add_dirichlet(pbc.E, perX, Te(i, j, k), pbc.E ? pbc.E->value : 0.0);
                        add_neumann(pbc.E, Te(i, j, k), dx_, +1);
                    }
                    // SOUTH/NORTH
                    if (!hasS && !perY)
                    {
                        add_dirichlet(pbc.S, perY, Ts(i, j, k), pbc.S ? pbc.S->value : 0.0);
                        add_neumann(pbc.S, Ts(i, j, k), dy_, -1);
                    }
                    if (!hasN && !perY)
                    {
                        add_dirichlet(pbc.N, perY, Tn(i, j, k), pbc.N ? pbc.N->value : 0.0);
                        add_neumann(pbc.N, Tn(i, j, k), dy_, +1);
                    }
                    // BOTTOM/TOP
                    if (!hasB && !perZ)
                    {
                        add_dirichlet(pbc.B, perZ, Tb(i, j, k), pbc.B ? pbc.B->value : 0.0);
                        add_neumann(pbc.B, Tb(i, j, k), dz_, -1);
                    }
                    if (!hasT && !perZ)
                    {
                        add_dirichlet(pbc.T, perZ, Tt(i, j, k), pbc.T ? pbc.T->value : 0.0);
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
                    const int ii = (int) (i - i0) + ng;
                    const int jj = (int) (j - j0) + ng;
                    const int kk = (int) (k - k0) + ng;
                    const std::size_t c =
                        std::size_t(ii) +
                        std::size_t(nxc_tot) *
                            (std::size_t(jj) + std::size_t(nyc_tot) * std::size_t(kk));
                    xarr[k][j][i] = impl_->p_host[c];
                }
        DMDAVecRestoreArray(impl_->da_fine, impl_->x, &xarr);
    }

    PetscCallAbort(impl_->comm, KSPSetInitialGuessNonzero(impl_->ksp, PETSC_TRUE));

    // Remove constant component from BOTH b and x when a nullspace is present (pure Neumann).
    // This keeps the singular system well-posed for Krylov and avoids stalls.
    {
        // Use the vectors' communicator (same as impl.comm)
        MPI_Comm vcomm = impl_->comm;
        PetscCallAbort(vcomm, PetscObjectGetComm((PetscObject) impl_->b, &vcomm));

        MatNullSpace nsp = NULL;
        PetscCallAbort(vcomm, MatGetNullSpace(impl_->A_shell, &nsp)); // may be NULL

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
        KSPSetTolerances(impl_->ksp, rtol, atol, 1e50, iters_);
    }
    // solve
    MPI_Comm pcomm;
    PetscObjectGetComm((PetscObject) impl_->ksp, &pcomm);
    PetscCallAbort(pcomm, KSPSolve(impl_->ksp, impl_->b, impl_->x));

    // ---- One-shot status line after solve ----
    {
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        PetscInt its = 0;
        PetscReal rnorm = 0.0;
        KSPGetConvergedReason(impl_->ksp, &reason);
        KSPGetIterationNumber(impl_->ksp, &its);
        KSPGetResidualNorm(impl_->ksp, &rnorm);
        if (reason < 0)
        {
            LOGE("[poisson] solve failed: reason=%d  iters=%d  ||r||=%.6e\n", (int) reason,
                 (int) its, (double) rnorm);
        }
        else
        {
            LOGD("[poisson] converged: reason=%d  iters=%d  ||r||=%.6e\n", (int) reason, (int) its,
                 (double) rnorm);
        }
    }

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
                    const int ii = (int) (i - i0) + ng;
                    const int jj = (int) (j - j0) + ng;
                    const int kk = (int) (k - k0) + ng;
                    const std::size_t c =
                        std::size_t(ii) +
                        std::size_t(nxc_tot) *
                            (std::size_t(jj) + std::size_t(nyc_tot) * std::size_t(kk));
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
        core::master::exchange_named_fields(fields, *tile.mesh, mpi_comm_, {"p"});
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
    const int iters = std::stoi(get("iters", "400"));
    const Params p = parse_params(kv);
    auto obj = std::make_shared<PressurePoisson>(rho, dx, dy, dz, iters, rc.mpi_comm, p.bcs);
    // propagate div_tol from user params if present (fallback 1e-7)
    obj->set_div_tol(std::stod(get("div_tol", "1e-7")));
    return obj;
}

} // namespace fluids
