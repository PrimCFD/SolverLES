module fluids_kernels
  use, intrinsic :: iso_c_binding
  implicit none
contains

  subroutine fluids_kernels_free_scratch() bind(C, name="fluids_kernels_free_scratch")
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine fluids_kernels_free_scratch

  subroutine sgs_smagorinsky_mac_c(u, v, w, &
                                   nxu_tot, nyu_tot, nzu_tot, &
                                   nxv_tot, nyv_tot, nzv_tot, &
                                   nxw_tot, nyw_tot, nzw_tot, &
                                   nxc_tot, nyc_tot, nzc_tot, &
                                   ng, dx, dy, dz, Cs, nu_t_c) bind(C, name="sgs_smagorinsky_mac_c")
    use, intrinsic :: iso_c_binding
    implicit none
    ! Face-centered velocities
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot
    integer(c_int), value :: nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot
    ! Cell-centered output grid
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: dx, dy, dz, Cs
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: nu_t_c(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cC0, cC, cu0, cv0, cw0
    integer(c_size_t) :: cu_w, cu_e, cv_s, cv_n, cw_b, cw_t
    integer(c_size_t) :: t0
    real(c_double) :: fx, fy, fz, Cs2, Delta
    real(c_double) :: dudx, dvdy, dwdz, dudy, dvdx, dudz, dwdx, dvdz, dwdy
    real(c_double) :: Sxx, Syy, Szz, Sxy, Sxz, Syz, Smag

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    Cs2 = Cs*Cs
    Delta = (dx*dy*dz)**(1.0d0/3.0d0)

    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)

    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cC0,cC,cu0,cv0,cw0, &
    !$omp   cu_w,cu_e,cv_s,cv_n,cw_b,cw_t, dudx,dvdy,dwdz,dudy,dvdx,dudz,dwdx,dvdz,dwdy, &
    !$omp   Sxx,Syy,Szz,Sxy,Sxz,Syz,Smag, t0)
    do k = 0, nzc-1
      do j = 0, nyc-1
        ! base pointers for this (j,k)
        cC0 = 1_c_size_t+int(ng, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        cu0 = 1_c_size_t+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        cv0 = 1_c_size_t+int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        cw0 = 1_c_size_t+int(ng, c_size_t)+syw*int(j+ng, c_size_t)

        cC = cC0
        !$omp simd linear(cC:1)
        do i = 0, nxc-1
          ! Faces straddling the cell center (west/east, south/north, bottom/top)
          cu_w = cu0+int(ng+i, c_size_t)   ! u at i-1/2
          cu_e = cu_w+sxu   ! u at i+1/2

          cv_s = (1_c_size_t+int(ng+i, c_size_t)+syv* &
                  int(ng+j, c_size_t)+szv*int(k+ng, c_size_t))  ! v at j-1/2
          cv_n = cv_s+syv   ! v at j+1/2

          cw_b = (1_c_size_t+int(ng+i, c_size_t)+syw* &
                  int(j+ng, c_size_t)+szw*int(ng+k, c_size_t))  ! w at k-1/2
          cw_t = cw_b+szw    ! w at k+1/2

          ! Normal derivatives (centered from face differences)
          dudx = (u(cu_e)-u(cu_w))*fx
          dvdy = (v(cv_n)-v(cv_s))*fy
          dwdz = (w(cw_t)-w(cw_b))*fz

          ! Cross derivatives (average the face-wise central diffs to center)
          ! dudy: average over the two x-faces:
          dudy = 0.25d0*fy*((u(cu_w+syu)-u(cu_w-syu))+(u(cu_e+syu)-u(cu_e-syu)))
          ! dvdx: average over the two y-faces:
          dvdx = 0.25d0*fx*((v(cv_s+sxv)-v(cv_s-sxv))+(v(cv_n+sxv)-v(cv_n-sxv)))
          ! dudz: average over the two x-faces:
          dudz = 0.25d0*fz*((u(cu_w+szu)-u(cu_w-szu))+(u(cu_e+szu)-u(cu_e-szu)))
          ! dwdx: average over the two z-faces:
          dwdx = 0.25d0*fx*((w(cw_b+sxw)-w(cw_b-sxw))+(w(cw_t+sxw)-w(cw_t-sxw)))
          ! dvdz: average over the two y-faces:
          dvdz = 0.25d0*fz*((v(cv_s+szv)-v(cv_s-szv))+(v(cv_n+szv)-v(cv_n-szv)))
          ! dwdy: average over the two z-faces:
          dwdy = 0.25d0*fy*((w(cw_b+syw)-w(cw_b-syw))+(w(cw_t+syw)-w(cw_t-syw)))

          ! Strain tensor at the cell center
          Sxx = dudx; Syy = dvdy; Szz = dwdz
          Sxy = 0.5d0*(dudy+dvdx)
          Sxz = 0.5d0*(dudz+dwdx)
          Syz = 0.5d0*(dvdz+dwdy)

          Smag = sqrt(2.d0*(Sxx*Sxx+Syy*Syy+Szz*Szz) &
                      +4.d0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz))

          nu_t_c(cC) = (Cs2*Delta*Delta)*Smag

          cC = cC+sxc
        end do
      end do
    end do
  end subroutine sgs_smagorinsky_mac_c

  subroutine divergence_mac_c(u, v, w, &
                              nxu_tot, nyu_tot, nzu_tot, &
                              nxv_tot, nyv_tot, nzv_tot, &
                              nxw_tot, nyw_tot, nzw_tot, &
                              nxc_tot, nyc_tot, nzc_tot, &
                              ng, dx, dy, dz, div) bind(C, name="divergence_mac_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: div(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cC0, cC, cu0, cv0, cw0, cu_w, cu_e, cv_s, cv_n, cw_b, cw_t
    real(c_double) :: fx, fy, fz

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz

    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cC0,cC,cu0,cv0,cw0,cu_w,cu_e,cv_s,cv_n,cw_b,cw_t)
    do k = 0, nzc-1
      do j = 0, nyc-1
        ! base indices at (i=0) for this j,k
        cC0 = 1_c_size_t+int(ng, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        cu0 = 1_c_size_t+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        cv0 = 1_c_size_t+int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        cw0 = 1_c_size_t+int(ng, c_size_t)+syw*int(j+ng, c_size_t)

        cC = cC0
        !$omp simd linear(cC:1)
        do i = 0, nxc-1
          ! u-faces: west at i -> (ng+i), east -> +1
          cu_w = cu0+int(ng+i, c_size_t)
          cu_e = cu_w+sxu

          ! v-faces: south at j -> (ng+j), north -> +syv
          cv_s = (1_c_size_t+int(ng+i, c_size_t)+syv*int(ng+j, c_size_t)+szv*int(k+ng, c_size_t))
          cv_n = cv_s+syv

          ! w-faces: bottom at k -> (ng+k), top -> +szw
          cw_b = (1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng+k, c_size_t))
          cw_t = cw_b+szw

          div(cC) = (u(cu_e)-u(cu_w))*fx+(v(cv_n)-v(cv_s))*fy+(w(cw_t)-w(cw_b))*fz
          cC = cC+sxc
        end do
      end do
    end do
  end subroutine divergence_mac_c

  subroutine gradp_faces_c(p, nxc_tot, nyc_tot, nzc_tot, ng, dx, dy, dz, &
                           dpx_u, nxu_tot, nyu_tot, nzu_tot, &
                           dpy_v, nxv_tot, nyv_tot, nzv_tot, &
                           dpz_w, nxw_tot, nyw_tot, nzw_tot) bind(C, name="gradp_faces_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot, nxu_tot, nyu_tot, nzu_tot
    integer(c_int), value :: nxv_tot, nyv_tot, nzv_tot, nxw_tot, nyw_tot, nzw_tot, ng
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)  :: p(*)
    real(c_double), intent(out) :: dpx_u(*), dpy_v(*), dpz_w(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxc, syc, szc, sxu, syu, szu, sxv, syv, szv, sxw, syw, szw
    integer(c_size_t) :: cL, cR, cu, cv, cw, baseC, baseU, baseV, baseW

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)

    ! u-faces: interior faces i=1..nxc-1
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseU,cL,cR,cu)
    do k = 0, nzc-1
      do j = 0, nyc-1
        baseC = 1_c_size_t+int(ng, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        baseU = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 1, nxc-1
          cL = baseC+int(i-1, c_size_t)
          cR = cL+sxc
          cu = baseU+int(i, c_size_t)     ! write at i = 1..nxc-1 -> faces ng+1..ng+nxc-1
          dpx_u(cu) = (p(cR)-p(cL))/dx
        end do
      end do
    end do

    ! v-faces: interior faces j=1..nyc-1
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseV,cL,cR,cv)
    do k = 0, nzc-1
      do i = 0, nxc-1
        baseC = 1_c_size_t+int(ng+i, c_size_t)+szc*int(k+ng, c_size_t)
        baseV = 1_c_size_t+int(ng+i, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 1, nyc-1
          cL = baseC+syc*int(j-1+ng, c_size_t)
          cR = cL+syc
          cv = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng+j, c_size_t)+szv*int(k+ng, c_size_t)
          dpy_v(cv) = (p(cR)-p(cL))/dy
        end do
      end do
    end do

    ! w-faces: interior faces k=1..nzc-1
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseW,cL,cR,cw)
    do j = 0, nyc-1
      do i = 0, nxc-1
        baseC = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)
        baseW = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)
        do k = 1, nzc-1
          cL = baseC+szc*int(k-1+ng, c_size_t)
          cR = cL+szc
          cw = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng+k, c_size_t)
          dpz_w(cw) = (p(cR)-p(cL))/dz
        end do
      end do
    end do
  end subroutine gradp_faces_c

  subroutine correct_velocity_mac_c(u, v, w, dpx_u, dpy_v, dpz_w, &
                                    nxu_tot, nyu_tot, nzu_tot, &
                                    nxv_tot, nyv_tot, nzv_tot, &
                                    nxw_tot, nyw_tot, nzw_tot, &
                                    nxc_tot, nyc_tot, nzc_tot, ng, rho, dt) &
    bind(C, name="correct_velocity_mac_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: rho, dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: dpx_u(*), dpy_v(*), dpz_w(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw
    integer(c_size_t) :: cu
    real(c_double) :: fac

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    fac = dt/rho

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu)
    do k = 0, nzc-1
      do j = 0, nyc-1
        do i = 1, nxc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
          u(cu) = u(cu)-fac*dpx_u(cu)
        end do
      end do
    end do

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu)
    do k = 0, nzc-1
      do i = 0, nxc-1
        do j = 1, nyc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)
          v(cu) = v(cu)-fac*dpy_v(cu)
        end do
      end do
    end do

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu)
    do j = 0, nyc-1
      do i = 0, nxc-1
        do k = 1, nzc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)
          w(cu) = w(cu)-fac*dpz_w(cu)
        end do
      end do
    end do
  end subroutine correct_velocity_mac_c

  subroutine correct_velocity_varrho_mac_c(u, v, w, dpx_u, dpy_v, dpz_w, &
                                           nxu_tot, nyu_tot, nzu_tot, &
                                           nxv_tot, nyv_tot, nzv_tot, &
                                           nxw_tot, nyw_tot, nzw_tot, &
                                           nxc_tot, nyc_tot, nzc_tot, ng, rho_c, dt) &
    bind(C, name="correct_velocity_varrho_mac_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: dpx_u(*), dpy_v(*), dpz_w(*), rho_c(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxc, syc, szc, sxu, syu, szu, sxv, syv, szv, sxw, syw, szw
    integer(c_size_t) :: cu, cL, cR
    real(c_double), parameter :: eps_rho = 1.0d-300
    real(c_double) :: betaL, betaR, beta_f, fac

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)

    ! -------- u faces: (i+1/2, j, k) => left/right cells (i, j, k) and (i+1, j, k)
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu,cL,cR,betaL,betaR,beta_f,fac)
    do k = 0, nzc-1
      do j = 0, nyc-1
        do i = 1, nxc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i-1, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          cR = cL+sxc

          ! beta = 1/rho at adjacent centers (guard against tiny rho)
          betaL = 1.0d0/max(rho_c(cL), eps_rho)
          betaR = 1.0d0/max(rho_c(cR), eps_rho)

          ! face beta via harmonic mean (consistent with Poisson operator):
          ! beta_f = 2 * betaL * betaR / (betaL + betaR)
          beta_f = (2.0d0*betaL*betaR)/(betaL+betaR)

          fac = dt*beta_f
          u(cu) = u(cu)-fac*dpx_u(cu)
        end do
      end do
    end do

    ! -------- v faces: (i, j+1/2, k) => below/above cells (i, j, k) and (i, j+1, k)
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu,cL,cR,betaL,betaR,beta_f,fac)
    do k = 0, nzc-1
      do i = 0, nxc-1
        do j = 1, nyc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j-1+ng, c_size_t)+szc*int(k+ng, c_size_t)
          cR = cL+syc

          betaL = 1.0d0/max(rho_c(cL), eps_rho)
          betaR = 1.0d0/max(rho_c(cR), eps_rho)
          beta_f = (2.0d0*betaL*betaR)/(betaL+betaR)

          fac = dt*beta_f
          v(cu) = v(cu)-fac*dpy_v(cu)
        end do
      end do
    end do

    ! -------- w faces: (i, j, k+1/2) => below/above cells (i, j, k) and (i, j, k+1)
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu,cL,cR,betaL,betaR,beta_f,fac)
    do j = 0, nyc-1
      do i = 0, nxc-1
        do k = 1, nzc-1
          cu = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k-1+ng, c_size_t)
          cR = cL+szc

          betaL = 1.0d0/max(rho_c(cL), eps_rho)
          betaR = 1.0d0/max(rho_c(cR), eps_rho)
          beta_f = (2.0d0*betaL*betaR)/(betaL+betaR)

          fac = dt*beta_f
          w(cu) = w(cu)-fac*dpz_w(cu)
        end do
      end do
    end do
  end subroutine correct_velocity_varrho_mac_c

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

  pure function kk3_deriv(qm2, qm1, q0, qp1, s, h, blend) result(dq)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(in) :: qm2, qm1, q0, qp1, s, h, blend
    real(c_double) :: dq, up3, cen4, sigma

    ! 3rd-order upwind part (mirrored by sign of s)
    if (s >= 0.d0) then
      ! wind > 0 : {i-2,i-1,i,i+1}
      up3 = (2.0d0*qm1+3.0d0*q0-6.0d0*qp1+qm2)/(6.0d0*h)
    else
      ! wind < 0 : mirrored {i+2,i+1,i,i-1}
      up3 = (-2.0d0*qp1-3.0d0*q0+6.0d0*qm1-qm2)/(6.0d0*h)
    end if

    ! small central correction (KK-style blend)
    sigma = merge(1.0d0, -1.0d0, s >= 0.d0)
    cen4 = (qp1-qm1)/(2.0d0*h)-sigma*(qp1-2.0d0*q0+qm1)/(6.0d0*h)

    dq = (1.0d0-blend)*up3+blend*cen4
  end function kk3_deriv

  subroutine advect_velocity_kk3_c(u, v, w, nx_tot, ny_tot, nz_tot, ng, dx, dy, dz, blend, &
                                   nu_u, nu_v, nu_w) bind(C, name="advect_velocity_kk3_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nx_tot, ny_tot, nz_tot, ng
    real(c_double), value :: dx, dy, dz, blend
    ! input face-velocity components (same layout as existing kernels)
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    ! output tendencies N(u) = -(u·∇)u, same shapes as inputs
    real(c_double), intent(out) :: nu_u(*), nu_v(*), nu_w(*)

    integer :: i, j, k, nx, ny, nz
    integer(c_size_t) :: sx, sy, sz, c, c0, ip, im, jp, jm, kp, km
    real(c_double) :: fx, fy, fz
    real(c_double) :: um, up, vm, vp, wm, wp
    real(c_double) :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
    real(c_double) :: uc, vc, wc

    nx = nx_tot-2*ng; ny = ny_tot-2*ng; nz = nz_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz
    sx = 1_c_size_t
    sy = int(nx_tot, c_size_t)
    sz = sy*int(ny_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,ip,im,jp,jm,kp,km,  &
    !$omp   uc,vc,wc, um,up,vm,vp,wm,wp, du_dx,du_dy,du_dz, dv_dx,dv_dy,dv_dz, dw_dx,dw_dy,dw_dz)
    do k = 0, nz-1
      do j = 0, ny-1
        c0 = 1_c_size_t+int(ng, c_size_t)+sy*int(j+ng, c_size_t)+sz*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nx-1
          ip = c+sx; im = c-sx
          jp = c+sy; jm = c-sy
          kp = c+sz; km = c-sz

          uc = u(c); vc = v(c); wc = w(c)

          ! ---- x-derivatives (3rd-upwind with KK blend) ----
          du_dx = kk3_deriv(u(c-2*sx), u(im), u(c), u(ip), uc, dx, blend)
          dv_dx = kk3_deriv(v(c-2*sx), v(im), v(c), v(ip), uc, dx, blend)
          dw_dx = kk3_deriv(w(c-2*sx), w(im), w(c), w(ip), uc, dx, blend)

          ! ---- y-derivatives ----
          du_dy = kk3_deriv(u(c-2*sy), u(c-sy), u(c), u(c+sy), vc, dy, blend)
          dv_dy = kk3_deriv(v(c-2*sy), v(c-sy), v(c), v(c+sy), vc, dy, blend)
          dw_dy = kk3_deriv(w(c-2*sy), w(c-sy), w(c), w(c+sy), vc, dy, blend)

          ! ---- z-derivatives ----
          du_dz = kk3_deriv(u(c-2*sz), u(c-sz), u(c), u(c+sz), wc, dz, blend)
          dv_dz = kk3_deriv(v(c-2*sz), v(c-sz), v(c), v(c+sz), wc, dz, blend)
          dw_dz = kk3_deriv(w(c-2*sz), w(c-sz), w(c), w(c+sz), wc, dz, blend)

          ! N(u) = -(u·∇)u (component-wise)
          nu_u(c) = -(uc*du_dx+vc*du_dy+wc*du_dz)
          nu_v(c) = -(uc*dv_dx+vc*dv_dy+wc*dv_dz)
          nu_w(c) = -(uc*dw_dx+vc*dw_dy+wc*dw_dz)

          c = c+1_c_size_t
        end do
      end do
    end do

  end subroutine advect_velocity_kk3_c

  subroutine advect_velocity_kk3_mac_c( &
    u, nxu_tot, nyu_tot, nzu_tot, &
    v, nxv_tot, nyv_tot, nzv_tot, &
    w, nxw_tot, nyw_tot, nzw_tot, &
    ng, dx, dy, dz, blend, &
    Nu, Nv, Nw) bind(C, name="advect_velocity_kk3_mac_c")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot
    integer(c_int), value :: nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot
    integer(c_int), value :: ng
    real(c_double), value :: dx, dy, dz, blend
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: Nu(*), Nv(*), Nw(*)

    integer :: i, j, k, nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw
    integer(c_size_t) :: c, c0, ip, im, jp, jm, kp, km
    real(c_double)    :: uc, vc, wc
    real(c_double)    :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz

    ! sizes (interior counts)
    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng

    ! strides
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)

    !-------------------------------
    ! (1) N_u on u-faces
    !-------------------------------
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,ip,im,jp,jm, &
    !$omp& kp, km, uc, vc, wc, du_dx, du_dy, du_dz)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0 = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:1)
        do i = 0, nxu-1
          ip = c+sxu; im = c-sxu
          jp = c+syu; jm = c-syu
          kp = c+szu; km = c-szu

          ! u-component located at u-face:
          uc = u(c)
          ! at u-face (i,j,k): collocate v(i..i+1, j+1/2, k) and w(i..i+1, j, k+1/2) across x
          vc = 0.5d0*(v(1_c_size_t+int(ng+i, c_size_t) &
                        +syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)) &  ! v(i,   j+1/2, k)
                      +v(1_c_size_t+int(ng+i+1, c_size_t) &
                         +syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)))  ! v(i+1, j+1/2, k)

          wc = 0.5d0*(w(1_c_size_t+int(ng+i, c_size_t) &
                        +syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)) &  ! w(i,   j, k+1/2)
                      +w(1_c_size_t+int(ng+i+1, c_size_t) &
                         +syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)))  ! w(i+1, j, k+1/2)

          du_dx = kk3_deriv(u(c-2*sxu), u(im), u(c), u(ip), uc, dx, blend)
          du_dy = kk3_deriv(u(c-2*syu), u(jm), u(c), u(jp), vc, dy, blend)
          du_dz = kk3_deriv(u(c-2*szu), u(km), u(c), u(kp), wc, dz, blend)
          Nu(c) = -(uc*du_dx+vc*du_dy+wc*du_dz)
          c = c+sxu
        end do
      end do
    end do

    !-------------------------------
    ! (2) N_v on v-faces
    !-------------------------------
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,ip,im,jp,jm, &
    !$omp& kp, km, uc, vc, wc, dv_dx, dv_dy, dv_dz)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        c = c0
        !$omp simd linear(c:syv)
        do j = 0, nyv-1
          ip = c+sxv; im = c-sxv
          jp = c+syv; jm = c-syv
          kp = c+szv; km = c-szv

          vc = v(c)
          ! at v-face (i,j,k): collocate u(i+1/2, j..j+1, k) and w(i, j..j+1, k+1/2) across y
          uc = 0.5d0*(u(1_c_size_t+int(ng+i, c_size_t) &
                        +syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)) &  ! u(i+1/2, j,   k)
                      +u(1_c_size_t+int(ng+i, c_size_t) &
                         +syu*int(j+ng+1, c_size_t)+szu*int(k+ng, c_size_t)))  ! u(i+1/2, j+1, k)

          wc = 0.5d0*(w(1_c_size_t+int(ng+i, c_size_t) &
                        +syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)) &  ! w(i, j,   k+1/2)
                      +w(1_c_size_t+int(ng+i, c_size_t) &
                         +syw*int(j+ng+1, c_size_t)+szw*int(k+ng, c_size_t)))  ! w(i, j+1, k+1/2)

          dv_dx = kk3_deriv(v(c-2*sxv), v(im), v(c), v(ip), uc, dx, blend)
          dv_dy = kk3_deriv(v(c-2*syv), v(jm), v(c), v(jp), vc, dy, blend)
          dv_dz = kk3_deriv(v(c-2*szv), v(km), v(c), v(kp), wc, dz, blend)
          Nv(c) = -(uc*dv_dx+vc*dv_dy+wc*dv_dz)
          c = c+syv
        end do
      end do
    end do

    !-------------------------------
    ! (3) N_w on w-faces
    !-------------------------------
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,c0,c,ip,im,jp,jm, &
    !$omp& kp, km, uc, vc, wc, dw_dx, dw_dy, dw_dz)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        c = c0
        !$omp simd linear(c:szw)
        do k = 0, nzw-1
          ip = c+sxw; im = c-sxw
          jp = c+syw; jm = c-syw
          kp = c+szw; km = c-szw

          wc = w(c)
          ! at w-face (i,j,k): collocate u(i+1/2, j, k..k+1) and v(i, j+1/2, k..k+1) across z
          uc = 0.5d0*(u(1_c_size_t+int(ng+i, c_size_t) &
                        +syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)) &  ! u(i+1/2, j, k)
                      +u(1_c_size_t+int(ng+i, c_size_t) &
                         +syu*int(j+ng, c_size_t)+szu*int(k+ng+1, c_size_t)))  ! u(i+1/2, j, k+1)

          vc = 0.5d0*(v(1_c_size_t+int(ng+i, c_size_t) &
                        +syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)) &  ! v(i, j+1/2, k)
                      +v(1_c_size_t+int(ng+i, c_size_t) &
                         +syv*int(j+ng, c_size_t)+szv*int(k+ng+1, c_size_t)))  ! v(i, j+1/2, k+1)

          dw_dx = kk3_deriv(w(c-2*sxw), w(im), w(c), w(ip), uc, dx, blend)
          dw_dy = kk3_deriv(w(c-2*syw), w(jm), w(c), w(jp), vc, dy, blend)
          dw_dz = kk3_deriv(w(c-2*szw), w(km), w(c), w(kp), wc, dz, blend)
          Nw(c) = -(uc*dw_dx+vc*dw_dy+wc*dw_dz)
          c = c+szw
        end do
      end do
    end do
  end subroutine advect_velocity_kk3_mac_c

end module fluids_kernels
