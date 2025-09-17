module fluids_kernels
  use, intrinsic :: iso_c_binding
  implicit none
contains

  pure integer(c_size_t) function lin(i, j, k, nx, ny, nz) result(L)
    integer(c_int), value :: i, j, k, nx, ny, nz
    L = int(i, c_size_t)+int(nx, c_size_t)*(int(j, c_size_t)+int(ny, c_size_t)*int(k, c_size_t))
  end function lin

  subroutine sgs_smagorinsky_c(u, v, w, nx_tot, ny_tot, nz_tot, ng, dx, dy, &
                               dz, Cs, nu_t) bind(C, name="sgs_smagorinsky_c")
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, Cs
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: nu_t(*)
    integer :: i, j, k, nx, ny, nz
    real(c_double) :: fx, fy, fz
    real(c_double) :: Sxx, Syy, Szz, Sxy, Sxz, Syz, Smag, Cs2, Delta
    real(c_double) :: dudx, dvdy, dwdz, dudy, dvdx, dudz, dwdx, dvdz, dwdy
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    Cs2 = Cs*Cs
    Delta = (dx*dy*dz)**(1.0d0/3.0d0)
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          dudx = (u(ip)-u(im))*0.5d0*fx
          dvdy = (v(jp)-v(jm))*0.5d0*fy
          dwdz = (w(kp)-w(km))*0.5d0*fz

          dudy = (u(jp)-u(jm))*0.5d0*fy
          dvdx = (v(ip)-v(im))*0.5d0*fx

          dudz = (u(kp)-u(km))*0.5d0*fz
          dwdx = (w(ip)-w(im))*0.5d0*fx

          dvdz = (v(kp)-v(km))*0.5d0*fz
          dwdy = (w(jp)-w(jm))*0.5d0*fy

          Sxx = dudx; Syy = dvdy; Szz = dwdz
          Sxy = 0.5d0*(dudy+dvdx)
          Sxz = 0.5d0*(dudz+dwdx)
          Syz = 0.5d0*(dvdz+dwdy)

          Smag = sqrt(2.d0*(Sxx*Sxx+Syy*Syy+Szz*Szz) &
                      +4.d0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz))

          nu_t(c) = (Cs2*Delta*Delta)*Smag
        end do
      end do
    end do
  end subroutine sgs_smagorinsky_c

  subroutine divergence_c(u, v, w, nx_tot, ny_tot, nz_tot, ng, dx, dy, &
                          dz, div) bind(C, name="divergence_c")
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: div(*)
    integer :: i, j, k, nx, ny, nz
    real(c_double) :: fx, fy, fz
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          div(c) = (u(ip)-u(im))*0.5d0*fx+(v(jp)-v(jm))*0.5d0*fy+(w(kp)-w(km))*0.5d0*fz
        end do
      end do
    end do
  end subroutine divergence_c

  subroutine poisson_jacobi_c(rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, iters, p) &
    bind(C, name="poisson_jacobi_c")
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng, iters
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)    :: rhs(*)
    real(c_double), intent(inout) :: p(*)
    integer :: i, j, k, nx, ny, nz, iter
    real(c_double), allocatable :: pn(:)
    real(c_double) :: ax, ay, az, ap, invap
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz
    integer(c_size_t) :: Ntot

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    Ntot = int(nx_tot, c_size_t)*int(ny_tot, c_size_t)*int(nz_tot, c_size_t)
    allocate (pn(Ntot))
    pn = 0.0d0

    ax = 1.0d0/(dx*dx); ay = 1.0d0/(dy*dy); az = 1.0d0/(dz*dz)
    ap = 2.0d0*(ax+ay+az); invap = 1.0d0/ap

    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do iter = 1, iters
      do k = 0, nz-1
        do j = 0, ny-1
          do i = 0, nx-1
            c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
            ip = c+sx; im = c-sx
            jp = c+sy; jm = c-sy
            kp = c+sz; km = c-sz
            pn(c) = (ax*(p(ip)+p(im))+ay*(p(jp)+p(jm))+az*(p(kp)+p(km))-rhs(c))*invap
          end do
        end do
      end do
      ! copy back **interior only**; keep halos (they will be exchanged by the caller)
      do k = 0, nz-1
        do j = 0, ny-1
          do i = 0, nx-1
            c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
            p(c) = pn(c)
          end do
        end do
      end do
    end do

    deallocate (pn)
  end subroutine poisson_jacobi_c

  subroutine gradp_c(p, nx_tot, ny_tot, nz_tot, ng, dx, dy, &
                     dz, dpx, dpy, dpz) bind(C, name="gradp_c")
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)  :: p(*)
    real(c_double), intent(out) :: dpx(*), dpy(*), dpz(*)
    integer :: i, j, k, nx, ny, nz
    real(c_double) :: fx, fy, fz
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          dpx(c) = (p(ip)-p(im))*0.5d0*fx
          dpy(c) = (p(jp)-p(jm))*0.5d0*fy
          dpz(c) = (p(kp)-p(km))*0.5d0*fz
        end do
      end do
    end do
  end subroutine gradp_c

  subroutine correct_velocity_c(u, v, w, dpx, dpy, dpz, nx_tot, ny_tot, nz_tot, ng, rho, &
                                dt) bind(C, name="correct_velocity_c")
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: rho, dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: dpx(*), dpy(*), dpz(*)
    integer(c_size_t) :: N, q
    N = int(nx_tot, c_size_t)*int(ny_tot, c_size_t)*int(nz_tot, c_size_t)
    do q = 1, N
      u(q) = u(q)-dt/rho*dpx(q)
      v(q) = v(q)-dt/rho*dpy(q)
      w(q) = w(q)-dt/rho*dpz(q)
    end do
  end subroutine correct_velocity_c

  ! =========================
  ! FE: single explicit step
  ! =========================
  subroutine diffuse_velocity_fe_c(u, v, w, nu_eff, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dt, &
                                   uo, vo, wo) bind(C, name="diffuse_velocity_fe_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u(*), v(*), w(*), nu_eff(*)
    real(c_double), intent(out) :: uo(*), vo(*), wo(*)
    integer :: i, j, k, nx, ny, nz
    real(c_double) :: ax, ay, az, nu_loc
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          nu_loc = nu_eff(c)
          uo(c) = u(c)+nu_loc*(ax*(u(ip)-2d0*u(c)+u(im))+ay*(u(jp)-2d0*u(c)+u(jm)) &
                               +az*(u(kp)-2d0*u(c)+u(km)))
          vo(c) = v(c)+nu_loc*(ax*(v(ip)-2d0*v(c)+v(im))+ay*(v(jp)-2d0*v(c)+v(jm)) &
                               +az*(v(kp)-2d0*v(c)+v(km)))
          wo(c) = w(c)+nu_loc*(ax*(w(ip)-2d0*w(c)+w(im))+ay*(w(jp)-2d0*w(c)+w(jm)) &
                               +az*(w(kp)-2d0*w(c)+w(km)))
        end do
      end do
    end do
  end subroutine diffuse_velocity_fe_c

  ! =========================================
  ! BE: one Jacobi sweep (INTERIOR ONLY)
  !    u_rhs = u^n, u_iter = current iterate
  !    writes u_next on interior; halos untouched
  ! =========================================
  subroutine diffuse_velocity_be_sweep_c(u_rhs, v_rhs, w_rhs, u_iter, v_iter, w_iter, nu_eff, &
                                         nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dt, &
                                         u_next, v_next, w_next) &
    bind(C, name="diffuse_velocity_be_sweep_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u_rhs(*), v_rhs(*), w_rhs(*), &
                                   u_iter(*), v_iter(*), w_iter(*), nu_eff(*)
    real(c_double), intent(out) :: u_next(*), v_next(*), w_next(*)

    integer :: i, j, k, nx, ny, nz
    real(c_double) :: ax, ay, az, ap, invap, nu_loc
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          nu_loc = nu_eff(c)
          ap = 1d0+2d0*nu_loc*(ax+ay+az)
          invap = 1d0/ap

          u_next(c) = (u_rhs(c)+nu_loc* &
                       (ax*(u_iter(ip)+u_iter(im))+ay*(u_iter(jp)+u_iter(jm)) &
                        +az*(u_iter(kp)+u_iter(km))))*invap
          v_next(c) = (v_rhs(c)+nu_loc* &
                       (ax*(v_iter(ip)+v_iter(im))+ay*(v_iter(jp)+v_iter(jm)) &
                        +az*(v_iter(kp)+v_iter(km))))*invap
          w_next(c) = (w_rhs(c)+nu_loc* &
                       (ax*(w_iter(ip)+w_iter(im))+ay*(w_iter(jp)+w_iter(jm)) &
                        +az*(w_iter(kp)+w_iter(km))))*invap
        end do
      end do
    end do
  end subroutine diffuse_velocity_be_sweep_c

  subroutine diffuse_velocity_be_gs_color_c(u, v, w, u_rhs, v_rhs, w_rhs, nu_eff, &
                                            nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dt, color) &
    bind(C, name="diffuse_velocity_be_gs_color_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng, color
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: u_rhs(*), v_rhs(*), w_rhs(*), nu_eff(*)

    integer :: i, j, k, nx, ny, nz, parity
    real(c_double) :: ax, ay, az, ap, invap, nu_loc
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          parity = i+j+k
          if (mod(parity, 2) /= color) cycle  ! 0 = red, 1 = black

          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          nu_loc = nu_eff(c)
          ap = 1d0+2d0*nu_loc*(ax+ay+az)
          invap = 1d0/ap

          u(c) = (u_rhs(c)+nu_loc*(ax*(u(ip)+u(im))+ay*(u(jp)+u(jm))+az*(u(kp)+u(km))))*invap
          v(c) = (v_rhs(c)+nu_loc*(ax*(v(ip)+v(im))+ay*(v(jp)+v(jm))+az*(v(kp)+v(km))))*invap
          w(c) = (w_rhs(c)+nu_loc*(ax*(w(ip)+w(im))+ay*(w(jp)+w(jm))+az*(w(kp)+w(km))))*invap
        end do
      end do
    end do
  end subroutine diffuse_velocity_be_gs_color_c

  ! =========================================
  ! BE residual on INTERIOR (needs valid halos
  ! on u_next/v_next/w_next before calling)
  ! =========================================
  subroutine diffuse_velocity_be_residual_c(u_rhs, v_rhs, w_rhs, u_next, v_next, w_next, nu_eff, &
                                            nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, dt, &
                                            res2, rhs2) &
    bind(C, name="diffuse_velocity_be_residual_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u_rhs(*), v_rhs(*), w_rhs(*), &
                                   u_next(*), v_next(*), w_next(*), nu_eff(*)
    real(c_double), intent(out) :: res2, rhs2

    integer :: i, j, k, nx, ny, nz
    real(c_double) :: ax, ay, az, ap, nu_loc, r_u, r_v, r_w
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    res2 = 0d0; rhs2 = 0d0

    do k = 0, nz-1
      do j = 0, ny-1
        do i = 0, nx-1
          c = 1+lin(i+ng, j+ng, k+ng, nx_tot, ny_tot, nz_tot)
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          nu_loc = nu_eff(c)
          ap = 1d0+2d0*nu_loc*(ax+ay+az)

          r_u = u_rhs(c)-(ap*u_next(c)-nu_loc* &
                          (ax*(u_next(ip)+u_next(im))+ay*(u_next(jp)+u_next(jm)) &
                           +az*(u_next(kp)+u_next(km))))
          r_v = v_rhs(c)-(ap*v_next(c)-nu_loc* &
                          (ax*(v_next(ip)+v_next(im))+ay*(v_next(jp)+v_next(jm)) &
                           +az*(v_next(kp)+v_next(km))))
          r_w = w_rhs(c)-(ap*w_next(c)-nu_loc* &
                          (ax*(w_next(ip)+w_next(im))+ay*(w_next(jp)+w_next(jm)) &
                           +az*(w_next(kp)+w_next(km))))

          res2 = res2+r_u*r_u+r_v*r_v+r_w*r_w
          rhs2 = rhs2+u_rhs(c)*u_rhs(c)+v_rhs(c)*v_rhs(c)+w_rhs(c)*w_rhs(c)
        end do
      end do
    end do
  end subroutine diffuse_velocity_be_residual_c

end module fluids_kernels
