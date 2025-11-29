#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

    void sgs_smagorinsky(const double* u, const double* v, const double* w, int nxu_tot,
                               int nyu_tot, int nzu_tot, int nxv_tot, int nyv_tot, int nzv_tot,
                               int nxw_tot, int nyw_tot, int nzw_tot, int nxc_tot, int nyc_tot,
                               int nzc_tot, int ng, double dx, double dy, double dz, double Cs,
                               double* nu_t_c);

    void divergence(const double* u, const double* v, const double* w, int nxu_tot,
                          int nyu_tot, int nzu_tot, int nxv_tot, int nyv_tot, int nzv_tot,
                          int nxw_tot, int nyw_tot, int nzw_tot, int nxc_tot, int nyc_tot,
                          int nzc_tot, int ng, double dx, double dy, double dz, double* div);

    void grad_p_faces(const double* p, int nxc_tot, int nyc_tot, int nzc_tot, int ng, double dx,
                       double dy, double dz, double* dpx_u, int nxu_tot, int nyu_tot, int nzu_tot,
                       double* dpy_v, int nxv_tot, int nyv_tot, int nzv_tot, double* dpz_w,
                       int nxw_tot, int nyw_tot, int nzw_tot);

    void poisson_jacobi_varcoef_c(const double* rhs, const double* beta, int nx_tot, int ny_tot,
                                  int nz_tot, int ng, double dx, double dy, double dz, int iters,
                                  double* p);

    void correct_velocity_const_rho(double* u, double* v, double* w, const double* dpx_u,
                                const double* dpy_v, const double* dpz_w, int nxu_tot, int nyu_tot,
                                int nzu_tot, int nxv_tot, int nyv_tot, int nzv_tot, int nxw_tot,
                                int nyw_tot, int nzw_tot, int nxc_tot, int nyc_tot, int nzc_tot,
                                int ng, double rho, double dt);

    void correct_velocity_varrho(double* u, double* v, double* w, const double* dpx_u,
                                       const double* dpy_v, const double* dpz_w, int nxu_tot,
                                       int nyu_tot, int nzu_tot, int nxv_tot, int nyv_tot,
                                       int nzv_tot, int nxw_tot, int nyw_tot, int nzw_tot,
                                       int nxc_tot, int nyc_tot, int nzc_tot, int ng,
                                       const double* rho_c, double dt);

    // FE (one shot)
    void diffuse_fe(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                                   const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                                   const double* w, int nxw_tot, int nyw_tot, int nzw_tot,
                                   const double* nu_eff, int nxc_tot, int nyc_tot, int nzc_tot,
                                   int ng, double dx, double dy, double dz, double dt, double* uo,
                                   double* vo, double* wo);

    // BE: one Jacobi sweep (interior only)
    void diffuse_be_jacobi_sweep(
        const double* u_rhs, const double* v_rhs, const double* w_rhs, const double* u_iter,
        const double* v_iter, const double* w_iter, const double* nu_eff, int nxc_tot, int nyc_tot,
        int nzc_tot, int nxu_tot, int nyu_tot, int nzu_tot, int nxv_tot, int nyv_tot, int nzv_tot,
        int nxw_tot, int nyw_tot, int nzw_tot, int ng, double dx, double dy, double dz, double dt,
        double* u_next, double* v_next, double* w_next);

    // BE: one red–black Gauss/Streidel sweep (interior only)
    void diffuse_be_rbgs_color(double* u, double* v, double* w, const double* u_rhs,
                                            const double* v_rhs, const double* w_rhs,
                                            const double* nu_eff, int nxc_tot, int nyc_tot,
                                            int nzc_tot, int nxu_tot, int nyu_tot, int nzu_tot,
                                            int nxv_tot, int nyv_tot, int nzv_tot, int nxw_tot,
                                            int nyw_tot, int nzw_tot, int ng, double dx, double dy,
                                            double dz, double dt, int color /* 0=red, 1=black */);

    // BE: residual (needs valid halos on u_next, v_next, w_next)
    void diffuse_be_residual(const double* u_rhs, const double* v_rhs,
                                            const double* w_rhs, const double* u_next,
                                            const double* v_next, const double* w_next,
                                            const double* nu_eff, int nxc_tot, int nyc_tot,
                                            int nzc_tot, int nxu_tot, int nyu_tot, int nzu_tot,
                                            int nxv_tot, int nyv_tot, int nzv_tot, int nxw_tot,
                                            int nyw_tot, int nzw_tot, int ng, double dx, double dy,
                                            double dz, double dt, double* res2, double* rhs2);

    /* KK MAC advection: central 4th + (C/4)|U| * 4th-diff */
    void advect_kk3(const double* u, int nxu_tot, int nyu_tot, int nzu_tot,
                                   const double* v, int nxv_tot, int nyv_tot, int nzv_tot,
                                   const double* w, int nxw_tot, int nyw_tot, int nzw_tot, int ng,
                                   double dx, double dy, double dz, double* Nu, double* Nv,
                                   double* Nw);

#ifdef __cplusplus
}
#endif
