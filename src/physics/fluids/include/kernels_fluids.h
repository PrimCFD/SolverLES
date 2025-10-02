#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

    void fluids_kernels_free_scratch();

    void sgs_smagorinsky_c(const double* u, const double* v, const double* w, int nx_tot,
                           int ny_tot, int nz_tot, int ng, double dx, double dy, double dz,
                           double Cs, double* nu_t_out);

    void divergence_c(const double* u, const double* v, const double* w, int nx_tot, int ny_tot,
                      int nz_tot, int ng, double dx, double dy, double dz, double* div_out);

    void divergence_rhie_chow_c(const double* u, const double* v, const double* w, const double* p,
                                const double* rho, int nx_tot, int ny_tot, int nz_tot, int ng,
                                double dx, double dy, double dz, double dt, double* div);

    void gradp_c(const double* p, int nx_tot, int ny_tot, int nz_tot, int ng, double dx, double dy,
                 double dz, double* dpx, double* dpy, double* dpz);

    void poisson_jacobi_c(const double* rhs, int nx_tot, int ny_tot, int nz_tot, int ng, double dx,
                          double dy, double dz, int iters, double* p_io);

    void poisson_jacobi_varcoef_c(const double* rhs, const double* beta, int nx_tot, int ny_tot,
                                  int nz_tot, int ng, double dx, double dy, double dz, int iters,
                                  double* p);

    void correct_velocity_c(double* u, double* v, double* w, const double* dpx, const double* dpy,
                            const double* dpz, int nx_tot, int ny_tot, int nz_tot, int ng,
                            double rho, double dt);

    void correct_velocity_varrho_c(double* u, double* v, double* w, const double* dpx,
                                   const double* dpy, const double* dpz, int nx_tot, int ny_tot,
                                   int nz_tot, int ng, const double* rho, double dt);

    // FE (one shot)
    void diffuse_velocity_fe_c(const double* u, const double* v, const double* w,
                               const double* nu_eff, int nx_tot, int ny_tot, int nz_tot, int ng,
                               double dx, double dy, double dz, double dt, double* uo, double* vo,
                               double* wo);

    // BE: one Jacobi sweep (interior only)
    void diffuse_velocity_be_sweep_c(const double* u_rhs, const double* v_rhs, const double* w_rhs,
                                     const double* u_iter, const double* v_iter,
                                     const double* w_iter, const double* nu_eff, int nx_tot,
                                     int ny_tot, int nz_tot, int ng, double dx, double dy,
                                     double dz, double dt, double* u_next, double* v_next,
                                     double* w_next);

    // BE: one red–black Gauss/Streidel sweep (interior only)
    void diffuse_velocity_be_gs_color_c(double* u, double* v, double* w, const double* u_rhs,
                                        const double* v_rhs, const double* w_rhs,
                                        const double* nu_eff, int nx_tot, int ny_tot, int nz_tot,
                                        int ng, double dx, double dy, double dz, double dt,
                                        int color /* 0=red, 1=black */);

    // BE: residual (needs valid halos on u_next, v_next, w_next)
    void diffuse_velocity_be_residual_c(const double* u_rhs, const double* v_rhs,
                                        const double* w_rhs, const double* u_next,
                                        const double* v_next, const double* w_next,
                                        const double* nu_eff, int nx_tot, int ny_tot, int nz_tot,
                                        int ng, double dx, double dy, double dz, double dt,
                                        double* res2, double* rhs2);

#ifdef __cplusplus
}
#endif
