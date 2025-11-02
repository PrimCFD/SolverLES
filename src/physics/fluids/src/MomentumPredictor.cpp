#include "MomentumPredictor.hpp"
#include "master/FieldCatalog.hpp"
#include "master/HaloOps.hpp"
#include <algorithm> // std::transform, std::copy_n
#include <array>
#include <cctype> // std::tolower
#include <cmath>  // std::sqrt
#include <cstdlib>
#include <cstring> // std::memcpy
#include <numeric>
#include <stdexcept>
#include <string>
#include "kernels_fluids.h"
#include <mpi.h>

#include "memory/MpiBox.hpp"

static std::string to_lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

using namespace core::master;

namespace fluids
{

static inline double to_d(const std::string& s, double dflt)
{
    char* e = nullptr;
    double v = std::strtod(s.c_str(), &e);
    return (e && *e == 0) ? v : dflt;
}

Predictor::Predictor(double rho, double nu, double dx, double dy, double dz, Mode mode)
    : rho_(rho), nu_(nu), dx_(dx), dy_(dy), dz_(dz), mode_(mode)
{
    info_.name = "momentum_predictor";
    info_.phases = plugin::Phase::Interior;
}

std::shared_ptr<plugin::IAction> make_predictor(const plugin::KV& kv,
                                                const core::master::RunContext& rc)
{
    auto get = [&](const char* k, const char* dflt) -> std::string
    {
        if (auto it = kv.find(k); it != kv.end())
            return it->second;
        return dflt;
    };

    const std::string scheme = to_lower(get("time_scheme", "fe"));
    Predictor::Mode mode = Predictor::Mode::FE;
    if (scheme == "be" || scheme == "backward_euler")
    {
        mode = Predictor::Mode::BE;
    }
    else if (scheme == "abm3" || scheme == "ab3" || scheme == "ab3am3")
    {
        mode = Predictor::Mode::ABM3;
    }

    auto P = std::make_shared<Predictor>(
        to_d(get("rho", "1.0"), 1.0), to_d(get("nu", "1e-3"), 1e-3), to_d(get("dx", "1.0"), 1.0),
        to_d(get("dy", "1.0"), 1.0), to_d(get("dz", "1.0"), 1.0), mode);

    // BE inner solver knobs (default: 50 iters, 1e-8 rtol)
    const int pred_imp_iters =
        (int) std::strtol(get("pred_imp_max_iters", "50").c_str(), nullptr, 10);
    const double pred_imp_rtol = to_d(get("pred_imp_rtol", "1e-8"), 1e-8);
    const std::string pred_imp_solver_str =
        to_lower(get("pred_imp_solver", "rbgs")); // "jacobi"|"rbgs"
    Predictor::IMPSolver pred_imp_solver =
        pred_imp_solver_str == "jacobi" ? Predictor::IMPSolver::Jacobi : Predictor::IMPSolver::RBGS;
    P->set_pred_imp_controls(pred_imp_iters, pred_imp_rtol, rc.mpi_comm,
                             pred_imp_solver); // Keep MPI communicator (opaque) for halo_ops
    const int imp_order =
        std::clamp((int) std::strtol(get("pred_imp_order", "3").c_str(), nullptr, 10), 1, 3);
    P->set_imp_order(imp_order);
    // --- KK advection knobs ---
    // mode: "fixed" uses adv_kkC (or adv_kk_C); "dynamic" uses gamma/Cmax/speed
    const std::string kk_mode = to_lower(get("adv_kk_mode", "fixed")); // "fixed"|"dynamic"
    // Accept both "adv_kkC" and "adv_kk_C" (prefer explicit adv_kkC if present)
    double adv_kkC_fixed = to_d(get("adv_kkC", "nan"), std::numeric_limits<double>::quiet_NaN());
    if (std::isnan(adv_kkC_fixed))
    {
        adv_kkC_fixed = to_d(get("adv_kk_C", "0.3333333333333333"), 1.0 / 3.0);
    }
    const double adv_kk_gamma = to_d(get("adv_kk_gamma", "5e-3"), 5e-3); // target γ*
    const double adv_kk_Cmax = to_d(get("adv_kk_Cmax", "1.0"), 1.0);
    const std::string kk_speed = to_lower(get("adv_kk_speed", "linf")); // "linf"|"rms"

    // Apply KK settings
    P->set_adv_kkC(adv_kkC_fixed);
    P->set_adv_kk_mode(kk_mode == "dynamic");
    P->set_adv_kk_gamma(adv_kk_gamma);
    P->set_adv_kk_Cmax(adv_kk_Cmax);
    P->set_adv_kk_speed(kk_speed);
    const std::string adv = to_lower(get("advect", "on"));
    const bool adv_on = !(adv == "0" || adv == "false" || adv == "off" || adv == "no");
    P->set_advect_enabled(adv_on);

    return P;
}

// ---- ABM history helper ----
void Predictor::push_adv_hist(const std::vector<double>& Nu, const std::vector<double>& Nv,
                              const std::vector<double>& Nw)
{
    // Keep newest at front; trim to 3 entries.
    auto push3 = [](std::deque<std::vector<double>>& Q, const std::vector<double>& a)
    {
        Q.push_front(a);
        if (Q.size() > 3)
            Q.pop_back();
    };
    push3(Nu_hist_, Nu);
    push3(Nv_hist_, Nv);
    push3(Nw_hist_, Nw);
    adv_hist_len_ = (int) std::min<std::size_t>(
        3, std::min({Nu_hist_.size(), Nv_hist_.size(), Nw_hist_.size()}));
}

void Predictor::push_diff_hist(const std::vector<double>& Du, const std::vector<double>& Dv,
                               const std::vector<double>& Dw)
{
    auto push3 = [](std::deque<std::vector<double>>& Q, const std::vector<double>& a)
    {
        Q.push_front(a);
        if (Q.size() > 3)
            Q.pop_back();
    };
    push3(Du_hist_, Du);
    push3(Dv_hist_, Dv);
    push3(Dw_hist_, Dw);
    diff_hist_len_ = (int) std::min<std::size_t>(
        3, std::min({Du_hist_.size(), Dv_hist_.size(), Dw_hist_.size()}));
}

void Predictor::execute(const MeshTileView& tile, FieldCatalog& fields, double dt)
{
    if (!fields.contains("u") || !fields.contains("v") || !fields.contains("w"))
        throw std::runtime_error("[fluids.predictor] u/v/w required.");

    auto vu = fields.view("u");
    auto vv = fields.view("v");
    auto vw = fields.view("w");
    auto vp = fields.view("p");

    if (!tile.mesh)
        throw std::runtime_error("[fluids.predictor] tile.mesh is null; Scheduler must set it.");
    const auto& mesh = *tile.mesh;
    const int ng = mesh.ng;
    const int nx_i = mesh.local[0], ny_i = mesh.local[1], nz_i = mesh.local[2];
    const int nx_tot = nx_i + 2 * ng;
    const int ny_tot = ny_i + 2 * ng;
    const int nz_tot = nz_i + 2 * ng;
    const int nxc_tot = nx_tot, nyc_tot = ny_tot, nzc_tot = nz_tot;

    // Face totals implied by MAC staggering (normal axis +1 interior cell)
    const int nxu_tot = (nx_i + 1) + 2 * ng, nyu_tot = ny_tot, nzu_tot = nz_tot;
    const int nxv_tot = nx_tot, nyv_tot = (ny_i + 1) + 2 * ng, nzv_tot = nz_tot;
    const int nxw_tot = nx_tot, nyw_tot = ny_tot, nzw_tot = (nz_i + 1) + 2 * ng;

    // Sanity: ensure allocated views match the mesh-derived totals
    auto check = [](const char* name, const std::array<int, 3>& e, int ex, int ey, int ez)
    {
        if (e[0] != ex || e[1] != ey || e[2] != ez)
            throw std::runtime_error(std::string("[fluids.predictor] view '") + name +
                                     "' extents do not match mesh totals.");
    };
    check("p", vp.extents, nx_tot, ny_tot, nz_tot);
    check("u", vu.extents, nxu_tot, nyu_tot, nzu_tot);
    check("v", vv.extents, nxv_tot, nyv_tot, nzv_tot);
    check("w", vw.extents, nxw_tot, nyw_tot, nzw_tot);

    const std::size_t Ncenter = (std::size_t) nx_tot * ny_tot * nz_tot;

    // KK3 advection accesses (i-2,j-2,k-2): require at least 2 ghost cells.
    if (advect_enabled_ && (ng < 2))
    {
        throw std::runtime_error(
            "[fluids.predictor] advect=on requires ng>=2 (got ng=" + std::to_string(ng) + ").");
    }

    const std::size_t Nu = (size_t) nxu_tot * nyu_tot * nzu_tot;
    const std::size_t Nv = (size_t) nxv_tot * nyv_tot * nzv_tot;
    const std::size_t Nw = (size_t) nxw_tot * nyw_tot * nzw_tot;
    // Freeze RHS = u^n, v^n, w^n for this timestep
    std::vector<double> urhs(Nu), vrhs(Nv), wrhs(Nw);
    std::copy_n((const double*) vu.host_ptr, Nu, urhs.data());
    std::copy_n((const double*) vv.host_ptr, Nv, vrhs.data());
    std::copy_n((const double*) vw.host_ptr, Nw, wrhs.data());

    // Effective viscosity array: nu + nu_t (if present)
    std::vector<double> nu_eff(Ncenter, nu_);
    if (fields.contains("nu_t"))
    {
        auto vt = fields.view("nu_t");
        const double* t = static_cast<const double*>(vt.host_ptr);
        for (std::size_t q = 0; q < Ncenter; ++q)
            nu_eff[q] = nu_ + t[q];
    }

    // Scratch storage for FE or BE next iterate
    if (us_.size() != Nu || vs_.size() != Nv || ws_.size() != Nw)
    {
        us_.assign(Nu, 0.0);
        vs_.assign(Nv, 0.0);
        ws_.assign(Nw, 0.0);
        reset_abm_history();
    }

    // Ensure ghosts are up-to-date before any wide-stencil reads
    core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

    // ---- Advection: KK3 tendency ----
    // Dynamic KK coefficient: recompute each call from current speeds & dt
    if (advect_enabled_ && adv_kk_dynamic_)
    {
        auto speed_ref = [&](const void* ptr, std::size_t N) -> double
        {
            const double* a = static_cast<const double*>(ptr);
            if (kk_speed_ref_ == KkSpeedRef::Linf)
            {
                double m = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                    m = std::max(m, std::abs(a[i]));
                return m;
            }
            else
            {
                long double s2 = 0.0L;
                for (std::size_t i = 0; i < N; ++i)
                {
                    long double v = a[i];
                    s2 += v * v;
                }
                return std::sqrt(static_cast<double>(s2 / std::max<std::size_t>(1, N)));
            }
        };

        const std::size_t Nu_tot = (std::size_t) nxu_tot * nyu_tot * nzu_tot;
        const std::size_t Nv_tot = (std::size_t) nxv_tot * nyv_tot * nzv_tot;
        const std::size_t Nw_tot = (std::size_t) nxw_tot * nyw_tot * nzw_tot;

        // Reference speeds on MAC faces (after halo exchange above)
        const double Urefx = speed_ref(vu.host_ptr, Nu_tot);
        const double Urefy = speed_ref(vv.host_ptr, Nv_tot);
        const double Urefz = speed_ref(vw.host_ptr, Nw_tot);

        // Nominal Courant numbers per axis
        const double eps = 1e-14;
        const double Cox = Urefx * dt / std::max(dx_, eps);
        const double Coy = Urefy * dt / std::max(dy_, eps);
        const double Coz = Urefz * dt / std::max(dz_, eps);
        const double Co_nom = std::max(Cox, std::max(Coy, Coz));

        if (Co_nom > eps)
        {
            // γ_KK ≈ (C/4)*Co  ⇒  C = 4 γ*/Co   (cap to Cmax)
            adv_kkC_ = std::min(4.0 * adv_kk_gamma_ / Co_nom, adv_kk_Cmax_);
        }
        else
        {
            adv_kkC_ = 0.0;
        }
    }
    // If advection is disabled, clear ONLY advection history;
    // keep diffusion history for AM2/AM3 implicit RHS.
    if (!advect_enabled_)
    {
        reset_adv_history();
    }
    std::vector<double> Nu_t, Nv_t, Nw_t;
    if (advect_enabled_)
    {
        const std::size_t NuN = Nu;
        const std::size_t NvN = Nv;
        const std::size_t NwN = Nw;
        Nu_t.assign(NuN, 0.0);
        Nv_t.assign(NvN, 0.0);
        Nw_t.assign(NwN, 0.0);
        advect_velocity_kk3_mac_c(static_cast<const double*>(vu.host_ptr), nxu_tot, nyu_tot,
                                  nzu_tot, static_cast<const double*>(vv.host_ptr), nxv_tot,
                                  nyv_tot, nzv_tot, static_cast<const double*>(vw.host_ptr),
                                  nxw_tot, nyw_tot, nzw_tot, ng, dx_, dy_, dz_, adv_kkC_,
                                  Nu_t.data(), Nv_t.data(), Nw_t.data());
    }

    // ---- Diffusion operator apply D(u)=L u at faces for AM2/AM3 RHS terms ----
    // We can reuse the FE kernel with dt_apply=1: out = u + 1*L(u) => L(u) = out - u
    std::vector<double> Du_t, Dv_t, Dw_t;
    {
        const double dt_apply = 1.0;
        // temp buffers
        std::vector<double> u0(Nu), v0(Nv), w0(Nw);
        std::memcpy(u0.data(), vu.host_ptr, Nu * sizeof(double));
        std::memcpy(v0.data(), vv.host_ptr, Nv * sizeof(double));
        std::memcpy(w0.data(), vw.host_ptr, Nw * sizeof(double));
        diffuse_velocity_fe_mac_c(
            /*u in  */ u0.data(), nxu_tot, nyu_tot, nzu_tot,
            /*v in  */ v0.data(), nxv_tot, nyv_tot, nzv_tot,
            /*w in  */ w0.data(), nxw_tot, nyw_tot, nzw_tot,
            /*nu eff*/ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
            /*geom  */ ng, dx_, dy_, dz_, dt_apply,
            /*out   */ us_.data(), vs_.data(), ws_.data());
        Du_t.assign(Nu, 0.0);
        Dv_t.assign(Nv, 0.0);
        Dw_t.assign(Nw, 0.0);
        for (std::size_t i = 0; i < Nu; ++i)
            Du_t[i] = us_[i] - u0[i];
        for (std::size_t i = 0; i < Nv; ++i)
            Dv_t[i] = vs_[i] - v0[i];
        for (std::size_t i = 0; i < Nw; ++i)
            Dw_t[i] = ws_[i] - w0[i];
    }

    // ---------------- ABM3 advection comb (AB1/AB2/AB3 start-up automatically) ----------------
    auto combine_AB = [&](const std::deque<std::vector<double>>& H, std::vector<double>& out)
    {
        // H[0]=N^n, H[1]=N^{n-1}, H[2]=N^{n-2}
        // AB1:  1.0
        // AB2:  3/2, -1/2
        // AB3:  23/12, -16/12, 5/12
        const int h = adv_hist_len_;
        const double a0 = (h >= 1) ? (h == 1 ? 1.0 : (h == 2 ? 1.5 : 23.0 / 12.0)) : 0.0;
        const double a1 = (h >= 2) ? (h == 2 ? -0.5 : -16.0 / 12.0) : 0.0;
        const double a2 = (h >= 3) ? (5.0 / 12.0) : 0.0;
        const std::size_t n = out.size();
        const double* p0 = H.size() >= 1 ? H[0].data() : nullptr;
        const double* p1 = H.size() >= 2 ? H[1].data() : nullptr;
        const double* p2 = H.size() >= 3 ? H[2].data() : nullptr;
        for (std::size_t i = 0; i < n; ++i)
        {
            double acc = 0.0;
            if (p0)
                acc += a0 * p0[i];
            if (p1)
                acc += a1 * p1[i];
            if (p2)
                acc += a2 * p2[i];
            out[i] = acc;
        }
    };

    // Always push current tendencies if advection is enabled
    if (advect_enabled_)
    {
        push_adv_hist(Nu_t, Nv_t, Nw_t);
    }
    push_diff_hist(Du_t, Dv_t, Dw_t);

    // Build explicit-advection combination per mode
    std::vector<double> Nu_ab(Nu, 0.0), Nv_ab(Nv, 0.0), Nw_ab(Nw, 0.0);
    if (advect_enabled_)
    {
        if (mode_ == Mode::ABM3)
        {
            // AB1/AB2/AB3 depending on history length
            combine_AB(Nu_hist_, Nu_ab);
            combine_AB(Nv_hist_, Nv_ab);
            combine_AB(Nw_hist_, Nw_ab);
        }
        else
        {
            // Legacy paths use Euler (AB1)
            std::memcpy(Nu_ab.data(), Nu_t.data(), Nu * sizeof(double));
            std::memcpy(Nv_ab.data(), Nv_t.data(), Nv * sizeof(double));
            std::memcpy(Nw_ab.data(), Nw_t.data(), Nw * sizeof(double));
        }
    }

    double* u_ptr = static_cast<double*>(vu.host_ptr);
    double* v_ptr = static_cast<double*>(vv.host_ptr);
    double* w_ptr = static_cast<double*>(vw.host_ptr);

    // Pure FE path
    if (mode_ == Mode::FE)
    {
        // ---- FE: single explicit update ----
        diffuse_velocity_fe_mac_c(
            /*u in  */ u_ptr, nxu_tot, nyu_tot, nzu_tot,
            /*v in  */ v_ptr, nxv_tot, nyv_tot, nzv_tot,
            /*w in  */ w_ptr, nxw_tot, nyw_tot, nzw_tot,
            /*nu eff*/ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
            /*geom  */ ng, dx_, dy_, dz_, dt,
            /*out   */ us_.data(), vs_.data(), ws_.data());

        // Add explicit advection: AB1/AB2/AB3 combo already in Nu_ab/Nv_ab/Nw_ab
        if (advect_enabled_)
        {
            for (std::size_t q = 0; q < Nu; ++q)
                us_[q] += dt * Nu_ab[q];
            for (std::size_t q = 0; q < Nv; ++q)
                vs_[q] += dt * Nv_ab[q];
            for (std::size_t q = 0; q < Nw; ++q)
                ws_[q] += dt * Nw_ab[q];
        }

        std::memcpy(u_ptr, us_.data(), Nu * sizeof(double));
        std::memcpy(v_ptr, vs_.data(), Nv * sizeof(double));
        std::memcpy(w_ptr, ws_.data(), Nw * sizeof(double));
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});
    }

    else
    {

        // ---- BE: host-controlled solve (Jacobi or RBGS), centralized halos ----

        // Initial guess u^(0) := u^n; ensure halos are consistent before first sweep
        core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});

        // --- IMEX RHS build ---
        // Start from u^n
        // Add explicit convection combo (AB1/2/3 depending on history) if enabled
        if (advect_enabled_)
        {
            for (std::size_t q = 0; q < Nu; ++q)
                urhs[q] += dt * Nu_ab[q];
            for (std::size_t q = 0; q < Nv; ++q)
                vrhs[q] += dt * Nv_ab[q];
            for (std::size_t q = 0; q < Nw; ++q)
                wrhs[q] += dt * Nw_ab[q];
        }
        // Implicit diffusion diagonal and RHS (user-selectable order):
        // Order 1: BE      -> theta=1, RHS += 0
        // Order 2: CN/AM2  -> theta=1/2, RHS += (1/2) dt D^n
        // Order 3: AM3     -> theta=5/12, RHS += dt*(8/12 D^n - 1/12 D^{n-1})
        double theta = 1.0; // default (BE)
        if (mode_ == Mode::ABM3)
        {
            if (imp_order_ == 1)
            {
                theta = 1.0; // BE
            }
            else if (imp_order_ == 2)
            {
                theta = 0.5; // CN
                if (diff_hist_len_ >= 1)
                {
                    const auto& Du_n = Du_hist_[0];
                    const auto& Dv_n = Dv_hist_[0];
                    const auto& Dw_n = Dw_hist_[0];
                    for (std::size_t q = 0; q < Nu; ++q)
                        urhs[q] += dt * (0.5 * Du_n[q]);
                    for (std::size_t q = 0; q < Nv; ++q)
                        vrhs[q] += dt * (0.5 * Dv_n[q]);
                    for (std::size_t q = 0; q < Nw; ++q)
                        wrhs[q] += dt * (0.5 * Dw_n[q]);
                }
            }
            else /* imp_order_ == 3 */
            {
                if (diff_hist_len_ >= 2)
                {
                    theta = 5.0 / 12.0;
                    const auto& Du_n = Du_hist_[0];
                    const auto& Du_nm = Du_hist_[1];
                    const auto& Dv_n = Dv_hist_[0];
                    const auto& Dv_nm = Dv_hist_[1];
                    const auto& Dw_n = Dw_hist_[0];
                    const auto& Dw_nm = Dw_hist_[1];
                    for (std::size_t q = 0; q < Nu; ++q)
                        urhs[q] += dt * ((8.0 / 12.0) * Du_n[q] + (-1.0 / 12.0) * Du_nm[q]);
                    for (std::size_t q = 0; q < Nv; ++q)
                        vrhs[q] += dt * ((8.0 / 12.0) * Dv_n[q] + (-1.0 / 12.0) * Dv_nm[q]);
                    for (std::size_t q = 0; q < Nw; ++q)
                        wrhs[q] += dt * ((8.0 / 12.0) * Dw_n[q] + (-1.0 / 12.0) * Dw_nm[q]);
                }
                else if (diff_hist_len_ == 1)
                {
                    // graceful startup: fall back to CN for the very first step
                    theta = 0.5;
                    const auto& Du_n = Du_hist_[0];
                    const auto& Dv_n = Dv_hist_[0];
                    const auto& Dw_n = Dw_hist_[0];
                    for (std::size_t q = 0; q < Nu; ++q)
                        urhs[q] += dt * (0.5 * Du_n[q]);
                    for (std::size_t q = 0; q < Nv; ++q)
                        vrhs[q] += dt * (0.5 * Dv_n[q]);
                    for (std::size_t q = 0; q < Nw; ++q)
                        wrhs[q] += dt * (0.5 * Dw_n[q]);
                }
                else
                {
                    // no history yet -> BE startup
                    theta = 1.0;
                }
            }
        }
        const double dt_imp = theta * dt; // scale the implicit block to match AMk diagonal

        const int kmax = std::max(1, pred_imp_max_iters_);
        const double rtol = (pred_imp_rtol_ > 0.0 ? pred_imp_rtol_ : 1e-8);

        for (int k = 0; k < kmax; ++k)
        {
            if (pred_imp_solver_ == IMPSolver::Jacobi)
            {
                // ---------- Jacobi (as before) ----------
                diffuse_velocity_be_sweep_mac_c(
                    /*rhs   */ urhs.data(), vrhs.data(), wrhs.data(),
                    /*iter  */ u_ptr, v_ptr, w_ptr,
                    /*nu    */ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
                    /*u ext */ nxu_tot, nyu_tot, nzu_tot,
                    /*v ext */ nxv_tot, nyv_tot, nzv_tot,
                    /*w ext */ nxw_tot, nyw_tot, nzw_tot,
                    /*geom  */ ng, dx_, dy_, dz_, dt_imp,
                    /*next  */ us_.data(), vs_.data(), ws_.data());

                // swap into the field buffers for halo exchange
                std::memcpy(u_ptr, us_.data(), Nu * sizeof(double));
                std::memcpy(v_ptr, vs_.data(), Nv * sizeof(double));
                std::memcpy(w_ptr, ws_.data(), Nw * sizeof(double));
                core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});
            }
            else
            {
                // ---------- RBGS (in-place color sweeps) ----------
                // Red
                diffuse_velocity_be_gs_color_mac_c(
                    /*inout*/ u_ptr, v_ptr, w_ptr,
                    /*rhs  */ urhs.data(), vrhs.data(), wrhs.data(),
                    /*nu   */ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
                    /*u ext*/ nxu_tot, nyu_tot, nzu_tot,
                    /*v ext*/ nxv_tot, nyv_tot, nzv_tot,
                    /*w ext*/ nxw_tot, nyw_tot, nzw_tot,
                    /*geom */ ng, dx_, dy_, dz_, dt_imp, /*color*/ 0);
                core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});
                // Black
                diffuse_velocity_be_gs_color_mac_c(
                    /*inout*/ u_ptr, v_ptr, w_ptr,
                    /*rhs  */ urhs.data(), vrhs.data(), wrhs.data(),
                    /*nu   */ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
                    /*u ext*/ nxu_tot, nyu_tot, nzu_tot,
                    /*v ext*/ nxv_tot, nyv_tot, nzv_tot,
                    /*w ext*/ nxw_tot, nyw_tot, nzw_tot,
                    /*geom */ ng, dx_, dy_, dz_, dt_imp, /*color*/ 1);
                core::master::exchange_named_fields(fields, mesh, mpi_comm_, {"u", "v", "w"});
            }

            // Residual (needs valid halos; we just exchanged)
            double res2 = 0.0, rhs2 = 0.0;
            diffuse_velocity_be_residual_mac_c(
                /*rhs  */ urhs.data(), vrhs.data(), wrhs.data(),
                /*next */ u_ptr, v_ptr, w_ptr,
                /*nu   */ nu_eff.data(), nxc_tot, nyc_tot, nzc_tot,
                /*u ext*/ nxu_tot, nyu_tot, nzu_tot,
                /*v ext*/ nxv_tot, nyv_tot, nzv_tot,
                /*w ext*/ nxw_tot, nyw_tot, nzw_tot,
                /*geom */ ng, dx_, dy_, dz_, dt_imp, &res2, &rhs2);
            if (mpi_comm_)
            {
                MPI_Comm comm = mpi_unbox(mpi_comm_);
                double pair[2] = {res2, rhs2};
                MPI_Allreduce(MPI_IN_PLACE, pair, 2, MPI_DOUBLE, MPI_SUM, comm);
                res2 = pair[0];
                rhs2 = pair[1];
            }

            if (rhs2 > 0.0 && std::sqrt(res2 / rhs2) < rtol)
                break;
        }
    }
}

} // namespace fluids
