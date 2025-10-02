module fluids_kernels
  use, intrinsic :: iso_c_binding
  implicit none
  ! Reusable scratch for poisson_jacobi_c
  real(c_double), allocatable, save :: pn_scratch(:)
contains

  subroutine fluids_kernels_free_scratch() bind(C, name="fluids_kernels_free_scratch")
    use, intrinsic :: iso_c_binding
    implicit none
    if (allocated(pn_scratch)) deallocate (pn_scratch)
  end subroutine fluids_kernels_free_scratch

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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    Cs2 = Cs*Cs
    Delta = (dx*dy*dz)**(1.0d0/3.0d0)
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km, &
    !$omp   dudx,dvdy,dwdz,dudy,dvdx,dudz,dwdx,dvdz,dwdy,Sxx,Syy,Szz,Sxy,Sxz,Syz,Smag)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
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
          c = c+1_c_size_t
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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          div(c) = (u(ip)-u(im))*0.5d0*fx+(v(jp)-v(jm))*0.5d0*fy+(w(kp)-w(km))*0.5d0*fz
          c = c+1_c_size_t
        end do
      end do
    end do
  end subroutine divergence_c

  subroutine divergence_rhie_chow_c(u, v, w, p, rho, nx_tot, ny_tot, nz_tot, &
                                    ng, dx, dy, dz, dt, div) &
    bind(C, name="divergence_rhie_chow_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u(*), v(*), w(*), p(*), rho(*)
    real(c_double), intent(out) :: div(*)

    integer :: i, j, k, nx, ny, nz
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz
    integer(c_size_t) :: c0
    real(c_double) :: fx, fy, fz
    real(c_double) :: u_e, u_w, v_n, v_s, w_t, w_b
    real(c_double) :: rho_e, rho_w, rho_n, rho_s, rho_t, rho_b
    real(c_double) :: D_e, D_w, D_n, D_s, D_t, D_b

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz

    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,ip,im,jp,jm,kp,km,&
    !$omp   u_e,u_w,v_n,v_s,w_t,w_b,rho_e,rho_w,rho_n,rho_s,rho_t,rho_b,D_e,D_w,D_n,D_s,D_t,D_b)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          ! Face densities (simple arithmetic average)
          rho_e = 0.5d0*(rho(c)+rho(ip)); rho_w = 0.5d0*(rho(c)+rho(im))
          rho_n = 0.5d0*(rho(c)+rho(jp)); rho_s = 0.5d0*(rho(c)+rho(jm))
          rho_t = 0.5d0*(rho(c)+rho(kp)); rho_b = 0.5d0*(rho(c)+rho(km))

          D_e = dt/max(rho_e, 1.0d-300); D_w = dt/max(rho_w, 1.0d-300)
          D_n = dt/max(rho_n, 1.0d-300); D_s = dt/max(rho_s, 1.0d-300)
          D_t = dt/max(rho_t, 1.0d-300); D_b = dt/max(rho_b, 1.0d-300)

          ! RC face velocities (no area scaling; consistent with divergence_c)
          u_e = 0.5d0*(u(c)+u(ip))-D_e*(p(ip)-p(c))*fx
          u_w = 0.5d0*(u(im)+u(c))-D_w*(p(c)-p(im))*fx

          v_n = 0.5d0*(v(c)+v(jp))-D_n*(p(jp)-p(c))*fy
          v_s = 0.5d0*(v(jm)+v(c))-D_s*(p(c)-p(jm))*fy

          w_t = 0.5d0*(w(c)+w(kp))-D_t*(p(kp)-p(c))*fz
          w_b = 0.5d0*(w(km)+w(c))-D_b*(p(c)-p(km))*fz

          div(c) = (u_e-u_w)*fx+(v_n-v_s)*fy+(w_t-w_b)*fz
          c = c+1_c_size_t
        end do
      end do
    end do
  end subroutine divergence_rhie_chow_c

  subroutine poisson_jacobi_c(rhs, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, iters, p) &
    bind(C, name="poisson_jacobi_c")
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng, iters
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)    :: rhs(*)
    real(c_double), intent(inout) :: p(*)
    integer :: i, j, k, nx, ny, nz, iter
    real(c_double) :: ax, ay, az, ap, invap
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz
    integer(c_size_t) :: Ntot, c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    Ntot = int(nx_tot, c_size_t)*int(ny_tot, c_size_t)*int(nz_tot, c_size_t)

    ! Allocate or resize the module-static scratch once; reuse across calls.
    if (.not. allocated(pn_scratch)) then
      allocate (pn_scratch(Ntot))
    else if (size(pn_scratch) /= int(Ntot)) then
      deallocate (pn_scratch)
      allocate (pn_scratch(Ntot))
    end if

    ax = 1.0d0/(dx*dx); ay = 1.0d0/(dy*dy); az = 1.0d0/(dz*dz)
    ap = 2.0d0*(ax+ay+az); invap = 1.0d0/ap

    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel default(shared) private(iter,i,j,k,c,c0,ip,im,jp,jm,kp,km)
    do iter = 1, iters
      !$omp do collapse(2) schedule(static)
      do k = 0, nz-1
        do j = 0, ny-1
          c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
          c = c0
          !$omp simd linear(c:1)
          do i = 0, nx-1
            ip = c+sx; im = c-sx
            jp = c+sy; jm = c-sy
            kp = c+sz; km = c-sz
            pn_scratch(c) = (ax*(p(ip)+p(im))+ay*(p(jp)+p(jm))+az*(p(kp)+p(km))-rhs(c))*invap
            c = c+1_c_size_t
          end do
        end do
      end do

      !$omp do collapse(2) schedule(static)
      do k = 0, nz-1
        do j = 0, ny-1
          c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
          c = c0
          !$omp simd linear(c:1)
          do i = 0, nx-1
            p(c) = pn_scratch(c)
            c = c+1_c_size_t
          end do
        end do
      end do
    end do
    !$omp end parallel

  end subroutine poisson_jacobi_c

  subroutine poisson_jacobi_varcoef_c(rhs, beta, nx_tot, ny_tot, nz_tot, &
                                      ng, dx, dy, dz, iters, p) &
    bind(C, name="poisson_jacobi_varcoef_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng, iters
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)    :: rhs(*), beta(*)
    real(c_double), intent(inout) :: p(*)

    integer :: i, j, k, nx, ny, nz, iter
    real(c_double) :: ax, ay, az, denom, be, bw, bn, bs, bt, bb
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz, c0, Ntot

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    Ntot = int(nx_tot, c_size_t)*int(ny_tot, c_size_t)*int(nz_tot, c_size_t)

    if (.not. allocated(pn_scratch)) then
      allocate (pn_scratch(Ntot))
    else if (size(pn_scratch) /= int(Ntot)) then
      deallocate (pn_scratch); allocate (pn_scratch(Ntot))
    end if

    ax = 1.0d0/(dx*dx); ay = 1.0d0/(dy*dy); az = 1.0d0/(dz*dz)
    sx = 1_c_size_t; sy = int(nx_tot, c_size_t); sz = sy*int(ny_tot, c_size_t)

    !$omp parallel default(shared) private(iter,i,j,k,c,c0,ip,im,jp,jm,kp,km,be,bw,bn,bs,bt,bb,denom)
    do iter = 1, iters
      !$omp do collapse(2) schedule(static)
      do k = 0, nz-1
        do j = 0, ny-1
          c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
          c = c0
          !$omp simd linear(c:1)
          do i = 0, nx-1
            ip = c+sx; im = c-sx; jp = c+sy; jm = c-sy; kp = c+sz; km = c-sz
            be = 0.5d0*(beta(c)+beta(ip)); bw = 0.5d0*(beta(c)+beta(im))
            bn = 0.5d0*(beta(c)+beta(jp)); bs = 0.5d0*(beta(c)+beta(jm))
            bt = 0.5d0*(beta(c)+beta(kp)); bb = 0.5d0*(beta(c)+beta(km))

            denom = ax*(be+bw)+ay*(bn+bs)+az*(bt+bb)
            pn_scratch(c) = (ax*(be*p(ip)+bw*p(im))+ay*(bn*p(jp)+bs*p(jm)) &
                             +az*(bt*p(kp)+bb*p(km))-rhs(c))/denom
            c = c+1_c_size_t
          end do
        end do
      end do

      !$omp do collapse(2) schedule(static)
      do k = 0, nz-1
        do j = 0, ny-1
          c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
          c = c0
          !$omp simd linear(c:1)
          do i = 0, nx-1
            p(c) = pn_scratch(c)
            c = c+1_c_size_t
          end do
        end do
      end do
    end do
    !$omp end parallel
  end subroutine poisson_jacobi_varcoef_c

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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          dpx(c) = (p(ip)-p(im))*0.5d0*fx
          dpy(c) = (p(jp)-p(jm))*0.5d0*fy
          dpz(c) = (p(kp)-p(km))*0.5d0*fz
          c = c+1_c_size_t
        end do
      end do
    end do
  end subroutine gradp_c

  subroutine correct_velocity_c(u, v, w, dpx, dpy, dpz, nx_tot, &
                                ny_tot, nz_tot, ng, rho, dt) bind(C, name="correct_velocity_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value  :: rho, dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: dpx(*), dpy(*), dpz(*)
    integer :: i, j, k, nx, ny, nz
    integer(c_size_t) :: sx, sy, sz, c0, c
    real(c_double) :: fac

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)
    fac = dt/rho

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          u(c) = u(c)-fac*dpx(c)
          v(c) = v(c)-fac*dpy(c)
          w(c) = w(c)-fac*dpz(c)
          c = c+sx
        end do
      end do
    end do
  end subroutine correct_velocity_c

  subroutine correct_velocity_varrho_c(u, v, w, dpx, dpy, dpz, nx_tot, ny_tot, nz_tot, &
                                       ng, rho, dt) &
    bind(C, name="correct_velocity_varrho_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: dpx(*), dpy(*), dpz(*), rho(*)
    integer :: i, j, k, nx, ny, nz
    integer(c_size_t) :: sx, sy, sz, c0, c
    real(c_double) :: fac

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    sx = 1_c_size_t; sy = int(nx_tot, c_size_t); sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,fac)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          fac = dt/max(rho(c), 1.0d-300)
          u(c) = u(c)-fac*dpx(c)
          v(c) = v(c)-fac*dpy(c)
          w(c) = w(c)-fac*dpz(c)
          c = c+sx
        end do
      end do
    end do
  end subroutine correct_velocity_varrho_c

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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km,nu_loc)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
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
          c = c+1_c_size_t
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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km,nu_loc,ap,invap)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
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
          c = c+1_c_size_t
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

    integer :: i, j, k, nx, ny, nz
    real(c_double) :: ax, ay, az, ap, invap, nu_loc
    integer(c_size_t) :: c, ip, im, jp, jm, kp, km, sx, sy, sz
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km,nu_loc,ap,invap)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0

        !$omp simd linear(c:1)
        do i = 0, nx-1
          if (iand(i+j+k, 1) == int(color)) then
            ip = c+sx; im = c-sx
            jp = c+sy; jm = c-sy
            kp = c+sz; km = c-sz

            nu_loc = nu_eff(c)
            ap = 1d0+2d0*nu_loc*(ax+ay+az)
            invap = 1d0/ap

            u(c) = (u_rhs(c)+nu_loc*(ax*(u(ip)+u(im))+ay*(u(jp)+u(jm))+az*(u(kp)+u(km))))*invap
            v(c) = (v_rhs(c)+nu_loc*(ax*(v(ip)+v(im))+ay*(v(jp)+v(jm))+az*(v(kp)+v(km))))*invap
            w(c) = (w_rhs(c)+nu_loc*(ax*(w(ip)+w(im))+ay*(w(jp)+w(jm))+az*(w(kp)+w(km))))*invap
          end if
          c = c+1_c_size_t
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
    integer(c_size_t) :: c0

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    res2 = 0d0; rhs2 = 0d0

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c,c0,ip,im,jp,jm,kp,km,nu_loc,ap,r_u,r_v,r_w) &
    !$omp reduction(+:res2,rhs2)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1) reduction(+:res2,rhs2)
        do i = 0, nx-1
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
          c = c+1_c_size_t
        end do
      end do
    end do
  end subroutine diffuse_velocity_be_residual_c

end module fluids_kernels
