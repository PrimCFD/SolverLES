#pragma once

#include <cstddef>

namespace numerics::kernels {

// =================
// Low-level kernels
// =================

// Smagorinsky SGS viscosity on a MAC grid
void sgs_smagorinsky(const double* u, const double* v, const double* w,
                     int nxu_tot, int nyu_tot, int nzu_tot,
                     int nxv_tot, int nyv_tot, int nzv_tot,
                     int nxw_tot, int nyw_tot, int nzw_tot,
                     int nxc_tot, int nyc_tot, int nzc_tot,
                     int ng, double dx, double dy, double dz,
                     double Cs, double* nu_t_c);

// Divergence of MAC velocity to cell centres
void divergence(const double* u, const double* v, const double* w,
                int nxu_tot, int nyu_tot, int nzu_tot,
                int nxv_tot, int nyv_tot, int nzv_tot,
                int nxw_tot, int nyw_tot, int nzw_tot,
                int nxc_tot, int nyc_tot, int nzc_tot,
                int ng, double dx, double dy, double dz,
                double* div_out);

// Gradient of cell-centred scalar p to faces
void grad_p_faces(const double* p,
                  int nxc_tot, int nyc_tot, int nzc_tot,
                  int ng, double dx, double dy, double dz,
                  double* dpx_u, int nxu_tot, int nyu_tot, int nzu_tot,
                  double* dpy_v, int nxv_tot, int nyv_tot, int nzv_tot,
                  double* dpz_w, int nxw_tot, int nyw_tot, int nzw_tot);

// Velocity correction for constant density
void correct_velocity_const_rho(double* u, double* v, double* w,
                                const double* dpx_u,
                                const double* dpy_v,
                                const double* dpz_w,
                                int nxu_tot, int nyu_tot, int nzu_tot,
                                int nxv_tot, int nyv_tot, int nzv_tot,
                                int nxw_tot, int nyw_tot, int nzw_tot,
                                int nxc_tot, int nyc_tot, int nzc_tot,
                                int ng, double rho, double dt);

// Velocity correction for variable density (rho at centres)
void correct_velocity_varrho(double* u, double* v, double* w,
                             const double* dpx_u,
                             const double* dpy_v,
                             const double* dpz_w,
                             int nxu_tot, int nyu_tot, int nzu_tot,
                             int nxv_tot, int nyv_tot, int nzv_tot,
                             int nxw_tot, int nyw_tot, int nzw_tot,
                             int nxc_tot, int nyc_tot, int nzc_tot,
                             int ng, const double* rho_c, double dt);

// Explicit FE diffusion of MAC velocity with cell-centred nu_eff
void diffuse_fe(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                const double* w, int nxw_tot, int nyw_tot, int nzw_tot,
                const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
                int ng, double dx, double dy, double dz, double dt,
                double* uo, double* vo, double* wo);

// One Jacobi BE diffusion sweep on MAC velocity (INTERIOR only)
void diffuse_be_jacobi_sweep(
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* u_iter, const double* v_iter, const double* w_iter,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    double* u_next, double* v_next, double* w_next);

// One Gauss–Seidel colour BE diffusion sweep (in-place, colour = 0 or 1)
void diffuse_be_rbgs_color(
    double* u, double* v, double* w,
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    int color);

// BE diffusion residual on MAC velocity (INTERIOR)
void diffuse_be_residual(
    const double* u_rhs, const double* v_rhs, const double* w_rhs,
    const double* u_next, const double* v_next, const double* w_next,
    const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
    int nxu_tot, int nyu_tot, int nzu_tot,
    int nxv_tot, int nyv_tot, int nzv_tot,
    int nxw_tot, int nyw_tot, int nzw_tot,
    int ng, double dx, double dy, double dz, double dt,
    double& res2, double& rhs2);

// KK3 advection of MAC velocity (produces Nu,Nv,Nw at faces)
void advect_kk3(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                const double* w, int nxw_tot, int nyw_tot, int nzw_tot,
                int ng, double dx, double dy, double dz,
                double* Nu, double* Nv, double* Nw);


} // namespace numerics::kernels
