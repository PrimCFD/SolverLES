#include "MacOps.hpp"
#include "kernels.h"

namespace numerics::kernels {

void sgs_smagorinsky(const double* u, const double* v, const double* w,
                     int nxu_tot, int nyu_tot, int nzu_tot,
                     int nxv_tot, int nyv_tot, int nzv_tot,
                     int nxw_tot, int nyw_tot, int nzw_tot,
                     int nxc_tot, int nyc_tot, int nzc_tot,
                     int ng, double dx, double dy, double dz,
                     double Cs, double* nu_t_c)
{
    ::sgs_smagorinsky(u, v, w,
                          nxu_tot, nyu_tot, nzu_tot,
                          nxv_tot, nyv_tot, nzv_tot,
                          nxw_tot, nyw_tot, nzw_tot,
                          nxc_tot, nyc_tot, nzc_tot,
                          ng, dx, dy, dz, Cs, nu_t_c);
}

void divergence(const double* u, const double* v, const double* w,
                int nxu_tot, int nyu_tot, int nzu_tot,
                int nxv_tot, int nyv_tot, int nzv_tot,
                int nxw_tot, int nyw_tot, int nzw_tot,
                int nxc_tot, int nyc_tot, int nzc_tot,
                int ng, double dx, double dy, double dz,
                double* div_out)
{
    ::divergence(u, v, w,
                     nxu_tot, nyu_tot, nzu_tot,
                     nxv_tot, nyv_tot, nzv_tot,
                     nxw_tot, nyw_tot, nzw_tot,
                     nxc_tot, nyc_tot, nzc_tot,
                     ng, dx, dy, dz, div_out);
}

void grad_p_faces(const double* p,
                  int nxc_tot, int nyc_tot, int nzc_tot,
                  int ng, double dx, double dy, double dz,
                  double* dpx_u, int nxu_tot, int nyu_tot, int nzu_tot,
                  double* dpy_v, int nxv_tot, int nyv_tot, int nzv_tot,
                  double* dpz_w, int nxw_tot, int nyw_tot, int nzw_tot)
{
    ::grad_p_faces(p,
                  nxc_tot, nyc_tot, nzc_tot,
                  ng, dx, dy, dz,
                  dpx_u, nxu_tot, nyu_tot, nzu_tot,
                  dpy_v, nxv_tot, nyv_tot, nzv_tot,
                  dpz_w, nxw_tot, nyw_tot, nzw_tot);
}

void correct_velocity_const_rho(double* u, double* v, double* w,
                                const double* dpx_u,
                                const double* dpy_v,
                                const double* dpz_w,
                                int nxu_tot, int nyu_tot, int nzu_tot,
                                int nxv_tot, int nyv_tot, int nzv_tot,
                                int nxw_tot, int nyw_tot, int nzw_tot,
                                int nxc_tot, int nyc_tot, int nzc_tot,
                                int ng, double rho, double dt)
{
    ::correct_velocity_const_rho(u, v, w,
                           dpx_u, dpy_v, dpz_w,
                           nxu_tot, nyu_tot, nzu_tot,
                           nxv_tot, nyv_tot, nzv_tot,
                           nxw_tot, nyw_tot, nzw_tot,
                           nxc_tot, nyc_tot, nzc_tot,
                           ng, rho, dt);
}

void correct_velocity_varrho(double* u, double* v, double* w,
                             const double* dpx_u,
                             const double* dpy_v,
                             const double* dpz_w,
                             int nxu_tot, int nyu_tot, int nzu_tot,
                             int nxv_tot, int nyv_tot, int nzv_tot,
                             int nxw_tot, int nyw_tot, int nzw_tot,
                             int nxc_tot, int nyc_tot, int nzc_tot,
                             int ng, const double* rho_c, double dt)
{
    ::correct_velocity_varrho(u, v, w,
                                  dpx_u, dpy_v, dpz_w,
                                  nxu_tot, nyu_tot, nzu_tot,
                                  nxv_tot, nyv_tot, nzv_tot,
                                  nxw_tot, nyw_tot, nzw_tot,
                                  nxc_tot, nyc_tot, nzc_tot,
                                  ng, rho_c, dt);
}

void diffuse_fe(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                const double* w, int nxw_tot, int nyw_tot, int nzw_tot,
                const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
                int ng, double dx, double dy, double dz, double dt,
                double* uo, double* vo, double* wo)
{
    ::diffuse_fe(
        u, nxu_tot, nyu_tot, nzu_tot,
        v, nxv_tot, nyv_tot, nzv_tot,
        w, nxw_tot, nyw_tot, nzw_tot,
        nu_eff, nxc_tot, nyc_tot, nzc_tot,
        ng, dx, dy, dz, dt,
        uo, vo, wo);
}

void diffuse_be_jacobi_sweep(
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* u_iter, const double* v_iter, const double* w_iter,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    double* u_next, double* v_next, double* w_next)
{
    ::diffuse_be_jacobi_sweep(
        u_rhs, v_rhs, w_rhs,
        u_iter, v_iter, w_iter,
        nu_eff, nxc_tot, nyc_tot, nzc_tot,
        nxu_tot, nyu_tot, nzu_tot,
        nxv_tot, nyv_tot, nzv_tot,
        nxw_tot, nyw_tot, nzw_tot,
        ng, dx, dy, dz, dt,
        u_next, v_next, w_next);
}

void diffuse_be_rbgs_color(
    double* u, double* v, double* w,
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    int color)
{
    ::diffuse_be_rbgs_color(
        u, v, w,
        u_rhs, v_rhs, w_rhs,
        nu_eff, nxc_tot, nyc_tot, nzc_tot,
        nxu_tot, nyu_tot, nzu_tot,
        nxv_tot, nyv_tot, nzv_tot,
        nxw_tot, nyw_tot, nzw_tot,
        ng, dx, dy, dz, dt, color);
}

void diffuse_be_residual(
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* u_next, const double* v_next, const double* w_next,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    double& res2, double& rhs2)
{
    double res2_local = 0.0;
    double rhs2_local = 0.0;

    ::diffuse_be_residual(
        u_rhs, v_rhs, w_rhs,
        u_next, v_next, w_next,
        nu_eff, nxc_tot, nyc_tot, nzc_tot,
        nxu_tot, nyu_tot, nzu_tot,
        nxv_tot, nyv_tot, nzv_tot,
        nxw_tot, nyw_tot, nzw_tot,
        ng, dx, dy, dz, dt,
        &res2_local, &rhs2_local);

    res2 = res2_local;
    rhs2 = rhs2_local;
}

void advect_kk3(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                const double* w, int nxw_tot, int nyw_tot, int nzw_tot,
                int ng, double dx, double dy, double dz,
                double* Nu, double* Nv, double* Nw)
{
    ::advect_kk3(
        u, nxu_tot, nyu_tot, nzu_tot,
        v, nxv_tot, nyv_tot, nzv_tot,
        w, nxw_tot, nyw_tot, nzw_tot,
        ng, dx, dy, dz,
        Nu, Nv, Nw);
}

} // namespace numerics::kernels
