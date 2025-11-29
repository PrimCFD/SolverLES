module kernels
  use, intrinsic :: iso_c_binding
  implicit none
contains

  subroutine sgs_smagorinsky(u, v, w, &
                             nxu_tot, nyu_tot, nzu_tot, &
                             nxv_tot, nyv_tot, nzv_tot, &
                             nxw_tot, nyw_tot, nzw_tot, &
                             nxc_tot, nyc_tot, nzc_tot, &
                             ng, dx, dy, dz, Cs, nu_t_c) bind(C, name="sgs_smagorinsky")
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
    integer(c_size_t) :: cC0, cC, cu0
    integer(c_size_t) :: cu_w, cu_e, cv_s, cv_n, cw_b, cw_t
    integer(c_size_t) :: baseV_jk, baseW_jk
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

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cC0,cC,cu0, &
    !$omp   cu_w,cu_e,cv_s,cv_n,cw_b,cw_t, dudx,dvdy,dwdz,dudy,dvdx,dudz,dwdx,dvdz,dwdy, &
    !$omp   Sxx,Syy,Szz,Sxy,Sxz,Syz,Smag, t0)
    do k = 0, nzc-1
      do j = 0, nyc-1
        ! base pointers for this (j,k)
        cC0 = 1_c_size_t+int(ng, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        cu0 = 1_c_size_t+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)

        cC = cC0
        ! bases for v- and w-faces that only depend on (j,k) **and include x halo**
        baseV_jk = 1_c_size_t+int(ng, c_size_t) &
                   +syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)
        baseW_jk = 1_c_size_t+int(ng, c_size_t) &
                   +syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)
        !$omp simd linear(cC:sxc)
        do i = 0, nxc-1
          ! Faces straddling the cell center (west/east, south/north, bottom/top)
          cu_w = cu0+int(ng+i, c_size_t)   ! u at i-1/2
          cu_e = cu_w+sxu   ! u at i+1/2

          ! v at j-1/2 and j+1/2 using a consistent base + x offset
          cv_s = baseV_jk+int(i, c_size_t)
          cv_n = cv_s+syv   ! v at j+1/2

          ! w at k-1/2 and k+1/2 using a consistent base + x offset
          cw_b = baseW_jk+int(i, c_size_t)
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
  end subroutine sgs_smagorinsky

  subroutine divergence(u, v, w, &
                        nxu_tot, nyu_tot, nzu_tot, &
                        nxv_tot, nyv_tot, nzv_tot, &
                        nxw_tot, nyw_tot, nzw_tot, &
                        nxc_tot, nyc_tot, nzc_tot, &
                        ng, dx, dy, dz, div) bind(C, name="divergence")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: dx, dy, dz
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: div(*)

    integer :: i, j, k, nxc, nyc, nzc
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cC0, cC, cu0, cu_w, cu_e, cv_s, cv_n, cw_b, cw_t
    integer(c_size_t) :: baseV_jk, baseW_jk
    real(c_double) :: fx, fy, fz

    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng
    fx = 1.0d0/dx; fy = 1.0d0/dy; fz = 1.0d0/dz

    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cC0,cC,cu0,cu_w,cu_e,cv_s,cv_n,cw_b,cw_t)
    do k = 0, nzc-1
      do j = 0, nyc-1
        ! base indices at (i=0) for this j,k
        cC0 = 1_c_size_t+int(ng, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        cu0 = 1_c_size_t+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)

        cC = cC0
        ! bases for v- and w-faces; fold x-halo into base and add +i only
        baseV_jk = 1_c_size_t+int(ng, c_size_t) &
                   +syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)
        baseW_jk = 1_c_size_t+int(ng, c_size_t) &
                   +syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)
        !$omp simd linear(cC:sxc)
        do i = 0, nxc-1
          ! u-faces: west at i -> (ng+i), east -> +1
          cu_w = cu0+int(ng+i, c_size_t)
          cu_e = cu_w+sxu

          ! v-faces: use base + x offset for south/north
          cv_s = baseV_jk+int(i, c_size_t)
          cv_n = cv_s+syv

          ! w-faces: use base + x offset for bottom/top
          cw_b = baseW_jk+int(i, c_size_t)
          cw_t = cw_b+szw

          div(cC) = (u(cu_e)-u(cu_w))*fx+(v(cv_n)-v(cv_s))*fy+(w(cw_t)-w(cw_b))*fz
          cC = cC+sxc
        end do
      end do
    end do
  end subroutine divergence

  subroutine grad_p_faces(p, nxc_tot, nyc_tot, nzc_tot, ng, dx, dy, dz, &
                          dpx_u, nxu_tot, nyu_tot, nzu_tot, &
                          dpy_v, nxv_tot, nyv_tot, nzv_tot, &
                          dpz_w, nxw_tot, nyw_tot, nzw_tot) bind(C, name="grad_p_faces")
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

    ! u-faces: interior faces i=0..nxc (faces ng .. ng+nxc)
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseU,cL,cR,cu)
    do k = 0, nzc-1
      do j = 0, nyc-1
        baseC = 1_c_size_t+int(ng, c_size_t) &
                +syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
        baseU = 1_c_size_t+int(ng, c_size_t) &
                +syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 0, nxc
          cL = baseC+sxc*int(i, c_size_t)-sxc
          cR = cL+sxc
          cu = baseU+int(i, c_size_t)
          dpx_u(cu) = (p(cR)-p(cL))/dx
        end do
      end do
    end do

    ! v-faces: interior faces j=0..nyc (faces ng .. ng+nyc)
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseV,cL,cR,cv)
    do k = 0, nzc-1
      do i = 0, nxc-1
        baseC = 1_c_size_t+int(ng+i, c_size_t) &
                +szc*int(k+ng, c_size_t)
        baseV = 1_c_size_t+int(ng+i, c_size_t) &
                +syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 0, nyc
          cL = baseC+syc*int(j+ng, c_size_t)-syc
          cR = cL+syc
          cv = baseV+syv*int(j, c_size_t)
          dpy_v(cv) = (p(cR)-p(cL))/dy
        end do
      end do
    end do

    ! w-faces: interior faces k=0..nzc (faces ng .. ng+nzc)
    !$omp parallel do collapse(2) schedule(static) private(j,k,i,baseC,baseW,cL,cR,cw)
    do j = 0, nyc-1
      do i = 0, nxc-1
        baseC = 1_c_size_t+int(ng+i, c_size_t) &
                +syc*int(j+ng, c_size_t)
        baseW = 1_c_size_t+int(ng+i, c_size_t) &
                +syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        do k = 0, nzc
          cL = baseC+szc*int(k+ng, c_size_t)-szc
          cR = cL+szc
          cw = baseW+szw*int(k, c_size_t)
          dpz_w(cw) = (p(cR)-p(cL))/dz
        end do
      end do
    end do
  end subroutine grad_p_faces

  subroutine correct_velocity_const_rho(u, v, w, dpx_u, dpy_v, dpz_w, &
                                        nxu_tot, nyu_tot, nzu_tot, &
                                        nxv_tot, nyv_tot, nzv_tot, &
                                        nxw_tot, nyw_tot, nzw_tot, &
                                        nxc_tot, nyc_tot, nzc_tot, ng, rho, dt) &
    bind(C, name="correct_velocity_const_rho")
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
        do i = 0, nxc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
          u(cu) = u(cu)-fac*dpx_u(cu)
        end do
      end do
    end do

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu)
    do k = 0, nzc-1
      do i = 0, nxc-1
        do j = 0, nyc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)
          v(cu) = v(cu)-fac*dpy_v(cu)
        end do
      end do
    end do

    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu)
    do j = 0, nyc-1
      do i = 0, nxc-1
        do k = 0, nzc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)
          w(cu) = w(cu)-fac*dpz_w(cu)
        end do
      end do
    end do
  end subroutine correct_velocity_const_rho

  subroutine correct_velocity_varrho(u, v, w, dpx_u, dpy_v, dpz_w, &
                                     nxu_tot, nyu_tot, nzu_tot, &
                                     nxv_tot, nyv_tot, nzv_tot, &
                                     nxw_tot, nyw_tot, nzw_tot, &
                                     nxc_tot, nyc_tot, nzc_tot, ng, rho_c, dt) &
    bind(C, name="correct_velocity_varrho")
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
        do i = 0, nxc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i-1, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          cR = cL+sxc

          ! If the two adjacent centers have exactly the same density, match the
          ! constant-ρ path *bit-for-bit*: fac = dt / rho. Otherwise use the
          ! harmonic mean of beta = 1/rho (consistent with the Poisson operator).
          if (rho_c(cL) == rho_c(cR)) then
            fac = dt/rho_c(cL)
          else
            betaL = 1.0d0/max(rho_c(cL), eps_rho)
            betaR = 1.0d0/max(rho_c(cR), eps_rho)
            beta_f = (2.0d0*betaL*betaR)/max(betaL+betaR, eps_rho)
            fac = dt*beta_f
          end if

          u(cu) = u(cu)-fac*dpx_u(cu)
        end do
      end do
    end do

    ! -------- v faces: (i, j+1/2, k) => below/above cells (i, j, k) and (i, j+1, k)
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu,cL,cR,betaL,betaR,beta_f,fac)
    do k = 0, nzc-1
      do i = 0, nxc-1
        do j = 0, nyc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j-1+ng, c_size_t)+szc*int(k+ng, c_size_t)
          cR = cL+syc

          if (rho_c(cL) == rho_c(cR)) then
            fac = dt/rho_c(cL)
          else
            betaL = 1.0d0/max(rho_c(cL), eps_rho)
            betaR = 1.0d0/max(rho_c(cR), eps_rho)
            beta_f = (2.0d0*betaL*betaR)/(betaL+betaR)
            fac = dt*beta_f
          end if

          v(cu) = v(cu)-fac*dpy_v(cu)
        end do
      end do
    end do

    ! -------- w faces: (i, j, k+1/2) => below/above cells (i, j, k) and (i, j, k+1)
    !$omp parallel do collapse(2) schedule(static) private(i,j,k,cu,cL,cR,betaL,betaR,beta_f,fac)
    do j = 0, nyc-1
      do i = 0, nxc-1
        do k = 0, nzc
          cu = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)

          cL = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k-1+ng, c_size_t)
          cR = cL+szc

          if (rho_c(cL) == rho_c(cR)) then
            fac = dt/rho_c(cL)
          else
            betaL = 1.0d0/max(rho_c(cL), eps_rho)
            betaR = 1.0d0/max(rho_c(cR), eps_rho)
            beta_f = (2.0d0*betaL*betaR)/(betaL+betaR)
            fac = dt*beta_f
          end if

          w(cu) = w(cu)-fac*dpz_w(cu)
        end do
      end do
    end do
  end subroutine correct_velocity_varrho

  ! =========================
  ! FE: single explicit step
  ! =========================
  subroutine diffuse_fe( &
    u, nxu_tot, nyu_tot, nzu_tot, &
    v, nxv_tot, nyv_tot, nzv_tot, &
    w, nxw_tot, nyw_tot, nzw_tot, &
    nu_eff, nxc_tot, nyc_tot, nzc_tot, &
    ng, dx, dy, dz, dt, &
    uo, vo, wo) bind(C, name="diffuse_fe")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot
    integer(c_int), value :: nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u(*), v(*), w(*), nu_eff(*)
    real(c_double), intent(out) :: uo(*), vo(*), wo(*)

    integer :: i, j, k, nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw
    integer :: nxc, nyc, nzc
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cu, cv, cw, c0u, c0v, c0w
    integer(c_size_t) :: c_l, c_r, c_s, c_n, c_b, c_t, c_s2, c_n2, c_b2, c_t2
    real(c_double) :: ax, ay, az
    real(c_double) :: nu_x, nu_y, nu_z
    real(c_double) :: ue, uw, un, us, ut, ub

    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng
    nxc = nxc_tot-2*ng; nyc = nyc_tot-2*ng; nzc = nzc_tot-2*ng

    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)

    ! ---- U faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(i,j,k, c0u, cu, c_l, c_r, c_s, c_n, c_b, c_t, c_s2, c_n2, c_b2, c_t2, &
    !$omp&         nu_x, nu_y, nu_z, ue, uw, un, us, ut, ub)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0u = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 0, nxu-1
          cu = c0u+int(i, c_size_t)

          ue = u(cu+sxu); uw = u(cu-sxu)
          un = u(cu+syu); us = u(cu-syu)
          ut = u(cu+szu); ub = u(cu-szu)

          c_l = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_r = c_l+sxc

          c_s = c_l; c_n = c_s+syc
          c_s2 = c_r; c_n2 = c_s2+syc
          c_b = c_l; c_t = c_b+szc
          c_b2 = c_r; c_t2 = c_b2+szc

          nu_x = 0.5d0*(nu_eff(c_l)+nu_eff(c_r))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s2)+nu_eff(c_n2))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b2)+nu_eff(c_t2))

          uo(cu) = u(cu)+(ax*nu_x*(ue-2d0*u(cu)+uw) &
                          +ay*nu_y*(un-2d0*u(cu)+us) &
                          +az*nu_z*(ut-2d0*u(cu)+ub))
        end do
      end do
    end do

    ! ---- V faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(i,j,k, c0v, cv, c_s, c_n, c_b, c_t, c_l, c_r, &
    !$omp&         nu_x, nu_y, nu_z, ue, uw, un, us, ut, ub)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0v = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 0, nyv-1
          cv = c0v+syv*int(j, c_size_t)

          ue = v(cv+sxv); uw = v(cv-sxv)
          un = v(cv+syv); us = v(cv-syv)
          ut = v(cv+szv); ub = v(cv-szv)

          c_s = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_n = c_s+syc

          c_l = c_s; c_r = c_l+sxc
          c_b = c_s; c_t = c_b+szc

          nu_y = 0.5d0*(nu_eff(c_s)+nu_eff(c_n))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+syc)+nu_eff(c_r+syc))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b+syc)+nu_eff(c_t+syc))

          vo(cv) = v(cv)+(ax*nu_x*(ue-2d0*v(cv)+uw) &
                          +ay*nu_y*(un-2d0*v(cv)+us) &
                          +az*nu_z*(ut-2d0*v(cv)+ub))
        end do
      end do
    end do

    ! ---- W faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(i,j,k, c0w, cw, c_b, c_t, c_l, c_r, c_s, c_n, &
    !$omp&         nu_x, nu_y, nu_z, ue, uw, un, us, ut, ub)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0w = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        do k = 0, nzw-1
          cw = c0w+szw*int(k, c_size_t)

          ue = w(cw+sxw); uw = w(cw-sxw)
          un = w(cw+syw); us = w(cw-syw)
          ut = w(cw+szw); ub = w(cw-szw)

          c_b = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_t = c_b+szc

          c_l = c_b; c_r = c_l+sxc
          c_s = c_b; c_n = c_s+syc

          nu_z = 0.5d0*(nu_eff(c_b)+nu_eff(c_t))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+szc)+nu_eff(c_r+szc))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+szc)+nu_eff(c_n+szc))

          wo(cw) = w(cw)+(ax*nu_x*(ue-2d0*w(cw)+uw) &
                          +ay*nu_y*(un-2d0*w(cw)+us) &
                          +az*nu_z*(ut-2d0*w(cw)+ub))
        end do
      end do
    end do
  end subroutine diffuse_fe

  ! =========================================
  ! BE: one Jacobi sweep (INTERIOR ONLY)
  !    u_rhs = u^n, u_iter = current iterate
  !    writes u_next on interior; halos untouched
  ! =========================================
  subroutine diffuse_be_jacobi_sweep( &
    u_rhs, v_rhs, w_rhs, &
    u_iter, v_iter, w_iter, &
    nu_eff, nxc_tot, nyc_tot, nzc_tot, &
    nxu_tot, nyu_tot, nzu_tot, &
    nxv_tot, nyv_tot, nzv_tot, &
    nxw_tot, nyw_tot, nzw_tot, &
    ng, dx, dy, dz, dt, &
    u_next, v_next, w_next) bind(C, name="diffuse_be_jacobi_sweep")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in)  :: u_rhs(*), v_rhs(*), w_rhs(*)
    real(c_double), intent(in)  :: u_iter(*), v_iter(*), w_iter(*), nu_eff(*)
    real(c_double), intent(out) :: u_next(*), v_next(*), w_next(*)

    integer :: nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, i, j, k
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cu, cv, cw, c0
    integer(c_size_t) :: c_l, c_r, c_s, c_n, c_b, c_t, c_s2, c_n2, c_b2, c_t2
    real(c_double) :: ax, ay, az, ap, nu_x, nu_y, nu_z

    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)

    ! ---- U faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cu, c_l, c_r, c_s, c_n, c_b, c_t, c_s2, c_n2, nu_x, nu_y, nu_z, ap)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0 = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 0, nxu-1
          cu = c0+int(i, c_size_t)

          c_l = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_r = c_l+sxc

          c_s = c_l; c_n = c_s+syc
          c_s2 = c_r; c_n2 = c_s2+syc
          c_b = c_l; c_t = c_b+szc
          c_b2 = c_r; c_t2 = c_b2+szc

          nu_x = 0.5d0*(nu_eff(c_l)+nu_eff(c_r))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s2)+nu_eff(c_n2))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b2)+nu_eff(c_t2))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          u_next(cu) = (u_rhs(cu) &
                        +ax*nu_x*(u_iter(cu+sxu)+u_iter(cu-sxu)) &
                        +ay*nu_y*(u_iter(cu+syu)+u_iter(cu-syu)) &
                        +az*nu_z*(u_iter(cu+szu)+u_iter(cu-szu)))/ap
        end do
      end do
    end do

    ! ---- V faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cv, c_s, c_n, c_l, c_r, c_b, c_t, nu_x, nu_y, nu_z, ap)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 0, nyv-1
          cv = c0+syv*int(j, c_size_t)

          c_s = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_n = c_s+syc
          c_l = c_s; c_r = c_l+sxc
          c_b = c_s; c_t = c_b+szc

          nu_y = 0.5d0*(nu_eff(c_s)+nu_eff(c_n))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+syc)+nu_eff(c_r+syc))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b+syc)+nu_eff(c_t+syc))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          v_next(cv) = (v_rhs(cv) &
                        +ax*nu_x*(v_iter(cv+sxv)+v_iter(cv-sxv)) &
                        +ay*nu_y*(v_iter(cv+syv)+v_iter(cv-syv)) &
                        +az*nu_z*(v_iter(cv+szv)+v_iter(cv-szv)))/ap
        end do
      end do
    end do

    ! ---- W faces ----
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cw, c_b, c_t, c_l, c_r, c_s, c_n, nu_x, nu_y, nu_z, ap)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        do k = 0, nzw-1
          cw = c0+szw*int(k, c_size_t)

          c_b = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_t = c_b+szc
          c_l = c_b; c_r = c_l+sxc
          c_s = c_b; c_n = c_s+syc

          nu_z = 0.5d0*(nu_eff(c_b)+nu_eff(c_t))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+szc)+nu_eff(c_r+szc))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+szc)+nu_eff(c_n+szc))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          w_next(cw) = (w_rhs(cw) &
                        +ax*nu_x*(w_iter(cw+sxw)+w_iter(cw-sxw)) &
                        +ay*nu_y*(w_iter(cw+syw)+w_iter(cw-syw)) &
                        +az*nu_z*(w_iter(cw+szw)+w_iter(cw-szw)))/ap
        end do
      end do
    end do
  end subroutine diffuse_be_jacobi_sweep

  subroutine diffuse_be_rbgs_color( &
    u, v, w, u_rhs, v_rhs, w_rhs, &
    nu_eff, nxc_tot, nyc_tot, nzc_tot, &
    nxu_tot, nyu_tot, nzu_tot, &
    nxv_tot, nyv_tot, nzv_tot, &
    nxw_tot, nyw_tot, nzw_tot, &
    ng, dx, dy, dz, dt, color) bind(C, name="diffuse_be_rbgs_color")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, ng, color
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(inout) :: u(*), v(*), w(*)
    real(c_double), intent(in)    :: u_rhs(*), v_rhs(*), w_rhs(*), nu_eff(*)

    integer :: nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, i, j, k
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: c, c0, c_l, c_r, c_s, c_n, c_b, c_t, c_b2, c_t2
    real(c_double) :: ax, ay, az, ap, nu_x, nu_y, nu_z

    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)

    ! U faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, c, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0 = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 0, nxu-1
          c = c0+int(i, c_size_t)
          if (iand(i+j+k, 1) == int(color)) then
            c_l = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
            c_r = c_l+sxc
            c_s = c_l; c_n = c_s+syc
            c_b = c_l; c_t = c_b+szc
            c_b2 = c_r; c_t2 = c_b2+szc

            nu_x = 0.5d0*(nu_eff(c_l)+nu_eff(c_r))
            nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+sxc)+nu_eff(c_n+sxc))
            nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b2)+nu_eff(c_t2))

            ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

            u(c) = (u_rhs(c) &
                    +ax*nu_x*(u(c+sxu)+u(c-sxu)) &
                    +ay*nu_y*(u(c+syu)+u(c-syu)) &
                    +az*nu_z*(u(c+szu)+u(c-szu)))/ap
          end if
        end do
      end do
    end do

    ! V faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, c, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 0, nyv-1
          c = c0+syv*int(j, c_size_t)
          if (iand(i+j+k, 1) == int(color)) then
            c_s = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
            c_n = c_s+syc
            c_l = c_s; c_r = c_l+sxc
            c_b = c_s; c_t = c_b+szc

            nu_y = 0.5d0*(nu_eff(c_s)+nu_eff(c_n))
            nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+syc)+nu_eff(c_r+syc))
            nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b+syc)+nu_eff(c_t+syc))

            ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

            v(c) = (v_rhs(c) &
                    +ax*nu_x*(v(c+sxv)+v(c-sxv)) &
                    +ay*nu_y*(v(c+syv)+v(c-syv)) &
                    +az*nu_z*(v(c+szv)+v(c-szv)))/ap
          end if
        end do
      end do
    end do

    ! W faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, c, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        do k = 0, nzw-1
          c = c0+szw*int(k, c_size_t)
          if (iand(i+j+k, 1) == int(color)) then
            c_b = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
            c_t = c_b+szc
            c_l = c_b; c_r = c_l+sxc
            c_s = c_b; c_n = c_s+syc

            nu_z = 0.5d0*(nu_eff(c_b)+nu_eff(c_t))
            nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+szc)+nu_eff(c_r+szc))
            nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+szc)+nu_eff(c_n+szc))

            ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

            w(c) = (w_rhs(c) &
                    +ax*nu_x*(w(c+sxw)+w(c-sxw)) &
                    +ay*nu_y*(w(c+syw)+w(c-syw)) &
                    +az*nu_z*(w(c+szw)+w(c-szw)))/ap
          end if
        end do
      end do
    end do
  end subroutine diffuse_be_rbgs_color

  ! =========================================
  ! BE residual on INTERIOR (needs valid halos
  ! on u_next/v_next/w_next before calling)
  ! =========================================
  subroutine diffuse_be_residual( &
    u_rhs, v_rhs, w_rhs, u_next, v_next, w_next, &
    nu_eff, nxc_tot, nyc_tot, nzc_tot, &
    nxu_tot, nyu_tot, nzu_tot, &
    nxv_tot, nyv_tot, nzv_tot, &
    nxw_tot, nyw_tot, nzw_tot, &
    ng, dx, dy, dz, dt, res2, rhs2) bind(C, name="diffuse_be_residual")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nxc_tot, nyc_tot, nzc_tot
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot, nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot, ng
    real(c_double), value :: dx, dy, dz, dt
    real(c_double), intent(in) :: u_rhs(*), v_rhs(*), w_rhs(*)
    real(c_double), intent(in) :: u_next(*), v_next(*), w_next(*), nu_eff(*)
    real(c_double), intent(out) :: res2, rhs2

    integer :: nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw, i, j, k
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw, sxc, syc, szc
    integer(c_size_t) :: cu, cv, cw, c0
    integer(c_size_t) :: c_l, c_r, c_s, c_n, c_b, c_t, c_b2, c_t2
    real(c_double) :: ax, ay, az, ap, nu_x, nu_y, nu_z, r

    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng
    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)
    sxc = 1_c_size_t; syc = int(nxc_tot, c_size_t); szc = syc*int(nyc_tot, c_size_t)

    ax = dt/(dx*dx); ay = dt/(dy*dy); az = dt/(dz*dz)
    res2 = 0d0; rhs2 = 0d0

    ! U faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cu, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap, r) &
    !$omp& reduction(+:res2,rhs2)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0 = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)
        do i = 0, nxu-1
          cu = c0+int(i, c_size_t)

          c_l = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_r = c_l+sxc
          c_s = c_l; c_n = c_s+syc
          c_b = c_l; c_t = c_b+szc
          c_b2 = c_r; c_t2 = c_b2+szc

          nu_x = 0.5d0*(nu_eff(c_l)+nu_eff(c_r))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+sxc)+nu_eff(c_n+sxc))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b2)+nu_eff(c_t2))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          r = u_rhs(cu)-(ap*u_next(cu) &
                         -ax*nu_x*(u_next(cu+sxu)+u_next(cu-sxu)) &
                         -ay*nu_y*(u_next(cu+syu)+u_next(cu-syu)) &
                         -az*nu_z*(u_next(cu+szu)+u_next(cu-szu)))
          res2 = res2+r*r
          rhs2 = rhs2+u_rhs(cu)*u_rhs(cu)
        end do
      end do
    end do

    ! V faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cv, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap, r) &
    !$omp& reduction(+:res2,rhs2)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)
        do j = 0, nyv-1
          cv = c0+syv*int(j, c_size_t)

          c_s = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_n = c_s+syc
          c_l = c_s; c_r = c_l+sxc
          c_b = c_s; c_t = c_b+szc

          nu_y = 0.5d0*(nu_eff(c_s)+nu_eff(c_n))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+syc)+nu_eff(c_r+syc))
          nu_z = 0.25d0*(nu_eff(c_b)+nu_eff(c_t)+nu_eff(c_b+syc)+nu_eff(c_t+syc))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          r = v_rhs(cv)-(ap*v_next(cv) &
                         -ax*nu_x*(v_next(cv+sxv)+v_next(cv-sxv)) &
                         -ay*nu_y*(v_next(cv+syv)+v_next(cv-syv)) &
                         -az*nu_z*(v_next(cv+szv)+v_next(cv-szv)))
          res2 = res2+r*r
          rhs2 = rhs2+v_rhs(cv)*v_rhs(cv)
        end do
      end do
    end do

    ! W faces
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(c0, cw, c_l, c_r, c_s, c_n, c_b, c_t, nu_x, nu_y, nu_z, ap, r) &
    !$omp& reduction(+:res2,rhs2)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)
        do k = 0, nzw-1
          cw = c0+szw*int(k, c_size_t)

          c_b = 1_c_size_t+int(ng+i, c_size_t)+syc*int(j+ng, c_size_t)+szc*int(k+ng, c_size_t)
          c_t = c_b+szc
          c_l = c_b; c_r = c_l+sxc
          c_s = c_b; c_n = c_s+syc

          nu_z = 0.5d0*(nu_eff(c_b)+nu_eff(c_t))
          nu_x = 0.25d0*(nu_eff(c_l)+nu_eff(c_r)+nu_eff(c_l+szc)+nu_eff(c_r+szc))
          nu_y = 0.25d0*(nu_eff(c_s)+nu_eff(c_n)+nu_eff(c_s+szc)+nu_eff(c_n+szc))

          ap = 1d0+2d0*(ax*nu_x+ay*nu_y+az*nu_z)

          r = w_rhs(cw)-(ap*w_next(cw) &
                         -ax*nu_x*(w_next(cw+sxw)+w_next(cw-sxw)) &
                         -ay*nu_y*(w_next(cw+syw)+w_next(cw-syw)) &
                         -az*nu_z*(w_next(cw+szw)+w_next(cw-szw)))
          res2 = res2+r*r
          rhs2 = rhs2+w_rhs(cw)*w_rhs(cw)
        end do
      end do
    end do
  end subroutine diffuse_be_residual

  pure function kk_alpha_from_phi(qm2, qm1, q0, qp1, qp2, h) result(alpha)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(in) :: qm2, qm1, q0, qp1, qp2, h
    real(c_double) :: alpha
    real(c_double) :: phi_p1, phi_p2, phi_m1, phi_m2, phi0
    real(c_double) :: dxx_p, dxx_m, lp, lm, num, den
    real(c_double), parameter :: eps = 1.0d-12
    phi0 = 0.0d0
    phi_p1 = (qp1-q0)/h
    phi_p2 = phi_p1+(qp2-qp1)/h
    phi_m1 = -(q0-qm1)/h
    phi_m2 = phi_m1-(qm1-qm2)/h
    dxx_p = (phi_p2-2.0d0*phi_p1+phi0)/2.0d0
    dxx_m = (phi0-2.0d0*phi_m1+phi_m2)/2.0d0
    lp = sqrt(1.0d0+phi_p1*phi_p1)
    lm = sqrt(1.0d0+phi_m1*phi_m1)
    num = abs(dxx_p+dxx_m)
    den = 4.0d0*(lp+lm)+eps
    alpha = min(1.0d0, num/den)
  end function kk_alpha_from_phi

  pure function kk_flux_hybrid_kk3(qm2, qm1, q0, qp1, qp2, U, h) result(flux)
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), intent(in) :: qm2, qm1, q0, qp1, qp2, U, h
    real(c_double)             :: flux
    real(c_double)             :: cen4, filt4, lap2, alpha, aU

    ! 4th-order central derivative (dq/dx)
    cen4 = (-qp2+8.0d0*qp1-8.0d0*qm1+qm2)/(12.0d0*h)

    ! second-difference and 4th-order filter in compact form
    lap2 = (qp1-2.0d0*q0+qm1)/(h*h)     ! -> (q_{k+1}-2q_k+q_{k-1})/h^2
    filt4 = (qp2-4.0d0*qp1+6.0d0*q0-4.0d0*qm1+qm2)/h

    alpha = kk_alpha_from_phi(qm2, qm1, q0, qp1, qp2, h)
    aU = abs(U)

    ! Adv3 + alpha*Dif1 + (1-alpha)*Dif3 (Eq. 10)
    ! Dif1 = -|U|/(2h) * (q_{k+1}-2q_k+q_{k-1}) = -|U| * 0.5 * (h*lap2)
    ! Dif3 =  |U|/(4h) * (q_{k+2}-4q_{k+1}+6q_k-4q_{k-1}+q_{k-2}) = |U| * 0.25 * filt4
    flux = U*cen4 &   ! Adv3
           +aU*(-0.5d0*alpha*(h*lap2) &   ! alpha * Dif1
                +0.25d0*(1.0d0-alpha)*filt4)    ! (1-alpha) * Dif3
  end function kk_flux_hybrid_kk3

  subroutine advect_kk3( &
    u, nxu_tot, nyu_tot, nzu_tot, &
    v, nxv_tot, nyv_tot, nzv_tot, &
    w, nxw_tot, nyw_tot, nzw_tot, &
    ng, dx, dy, dz, &
    Nu, Nv, Nw) bind(C, name="advect_kk3")
    use, intrinsic :: iso_c_binding
    implicit none
    ! sizes
    integer(c_int), value :: nxu_tot, nyu_tot, nzu_tot
    integer(c_int), value :: nxv_tot, nyv_tot, nzv_tot
    integer(c_int), value :: nxw_tot, nyw_tot, nzw_tot
    integer(c_int), value :: ng
    real(c_double), value :: dx, dy, dz
    ! fields
    real(c_double), intent(in)  :: u(*), v(*), w(*)
    real(c_double), intent(out) :: Nu(*), Nv(*), Nw(*)

    ! local sizes
    integer :: i, j, k, nxu, nyu, nzu, nxv, nyv, nzv, nxw, nyw, nzw
    ! strides
    integer(c_size_t) :: sxu, syu, szu, sxv, syv, szv, sxw, syw, szw
    ! rolling pointers / neighbors
    integer(c_size_t) :: c, c0, ip, im, jp, jm, kp, km, ipp, imm, jpp, jmm, kpp, kmm
    ! per-loop bases (predeclared to keep SIMD clean)
    integer(c_size_t) :: baseV_jk_u, baseW_jk_u
    integer(c_size_t) :: baseU_ik_v, baseW_ik_v
    integer(c_size_t) :: baseU_ij_w, baseV_ij_w
    ! face indices gathered by bases
    integer(c_size_t) :: ivL, ivR, iwB, iwT, iuB, iuT, iwB2, iwT2, iuL, iuU, ivL2, ivU2
    ! convective velocities at faces
    real(c_double) :: uc, vc, wc
    ! KK3 Fluxes
    real(c_double) :: Fx, Fy, Fz

    ! --------------------
    ! basic geometry
    ! --------------------
    nxu = nxu_tot-2*ng; nyu = nyu_tot-2*ng; nzu = nzu_tot-2*ng
    nxv = nxv_tot-2*ng; nyv = nyv_tot-2*ng; nzv = nzv_tot-2*ng
    nxw = nxw_tot-2*ng; nyw = nyw_tot-2*ng; nzw = nzw_tot-2*ng

    sxu = 1_c_size_t; syu = int(nxu_tot, c_size_t); szu = syu*int(nyu_tot, c_size_t)
    sxv = 1_c_size_t; syv = int(nxv_tot, c_size_t); szv = syv*int(nyv_tot, c_size_t)
    sxw = 1_c_size_t; syw = int(nxw_tot, c_size_t); szw = syw*int(nyw_tot, c_size_t)

    ! KK3 needs ±2 interior neighbors
    if (ng < 2_c_int) return

    !=======================================================
    ! (1) N_u on u-faces
    !=======================================================
    !$omp parallel do collapse(2) schedule(static) private( &
    !$omp   i,j,k,c0,c,ip,im,ipp,imm,jp,jm,jpp,jmm,kp,km,kpp,kmm, &
    !$omp   baseV_jk_u, baseW_jk_u, ivL,ivR,iwB,iwT, uc,vc,wc, &
    !$omp   Fx, Fy, Fz)
    do k = 0, nzu-1
      do j = 0, nyu-1
        c0 = 1_c_size_t+int(ng, c_size_t)+syu*int(j+ng, c_size_t)+szu*int(k+ng, c_size_t)

        ! bases that depend only on (j,k); fold x-halo
        baseV_jk_u = 1_c_size_t+int(ng, c_size_t)+syv*int(j+ng, c_size_t)+szv*int(k+ng, c_size_t)
        baseW_jk_u = 1_c_size_t+int(ng, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(k+ng, c_size_t)

        c = c0
        !$omp simd linear(c:sxu)
        do i = 0, nxu-1
          ip = c+sxu; im = c-sxu; ipp = c+2*sxu; imm = c-2*sxu
          jp = c+syu; jm = c-syu; jpp = c+2*syu; jmm = c-2*syu
          kp = c+szu; km = c-szu; kpp = c+2*szu; kmm = c-2*szu

          uc = u(c)

          ivL = baseV_jk_u+int(i, c_size_t)
          ivR = baseV_jk_u+int(i+1, c_size_t)
          iwB = baseW_jk_u+int(i, c_size_t)
          iwT = baseW_jk_u+int(i+1, c_size_t)

          vc = 0.5d0*(v(ivL)+v(ivR))
          wc = 0.5d0*(w(iwB)+w(iwT))

          Fx = kk_flux_hybrid_kk3(u(imm), u(im), u(c), u(ip), u(ipp), uc, dx)
          Fy = kk_flux_hybrid_kk3(u(jmm), u(jm), u(c), u(jp), u(jpp), vc, dy)
          Fz = kk_flux_hybrid_kk3(u(kmm), u(km), u(c), u(kp), u(kpp), wc, dz)

          Nu(c) = -(Fx+Fy+Fz)
          c = c+sxu
        end do
      end do
    end do

    !=======================================================
    ! (2) N_v on v-faces
    !=======================================================
    !$omp parallel do collapse(2) schedule(static) private( &
    !$omp   i,j,k,c0,c,ip,im,ipp,imm,jp,jm,jpp,jmm,kp,km,kpp,kmm, &
    !$omp   baseU_ik_v, baseW_ik_v, iuB,iuT,iwB2,iwT2, uc,vc,wc, &
    !$omp   Fx, Fy, Fz)
    do k = 0, nzv-1
      do i = 0, nxv-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syv*int(ng, c_size_t)+szv*int(k+ng, c_size_t)

        ! bases that depend only on (i,k); fold y-halo
        baseU_ik_v = 1_c_size_t+int(ng+i, c_size_t)+szu*int(k+ng, c_size_t)
        baseW_ik_v = 1_c_size_t+int(ng+i, c_size_t)+szw*int(k+ng, c_size_t)

        c = c0
        !$omp simd linear(c:syv)
        do j = 0, nyv-1
          ip = c+sxv; im = c-sxv; ipp = c+2*sxv; imm = c-2*sxv
          jp = c+syv; jm = c-syv; jpp = c+2*syv; jmm = c-2*syv
          kp = c+szv; km = c-szv; kpp = c+2*szv; kmm = c-2*szv

          vc = v(c)

          iuB = baseU_ik_v+syu*int(j+ng, c_size_t)
          iuT = baseU_ik_v+syu*int(j+ng+1, c_size_t)
          iwB2 = baseW_ik_v+syw*int(j+ng, c_size_t)
          iwT2 = baseW_ik_v+syw*int(j+ng+1, c_size_t)

          uc = 0.5d0*(u(iuB)+u(iuT))
          wc = 0.5d0*(w(iwB2)+w(iwT2))

          Fx = kk_flux_hybrid_kk3(v(imm), v(im), v(c), v(ip), v(ipp), uc, dx)
          Fy = kk_flux_hybrid_kk3(v(jmm), v(jm), v(c), v(jp), v(jpp), vc, dy)
          Fz = kk_flux_hybrid_kk3(v(kmm), v(km), v(c), v(kp), v(kpp), wc, dz)

          Nv(c) = -(Fx+Fy+Fz)
          c = c+syv
        end do
      end do
    end do

    !=======================================================
    ! (3) N_w on w-faces
    !=======================================================
    !$omp parallel do collapse(2) schedule(static) private( &
    !$omp   i,j,k,c0,c,ip,im,ipp,imm,jp,jm,jpp,jmm,kp,km,kpp,kmm, &
    !$omp   baseU_ij_w, baseV_ij_w, iuL,iuU,ivL2,ivU2, uc,vc,wc, &
    !$omp   Fx, Fy, Fz)
    do j = 0, nyw-1
      do i = 0, nxw-1
        c0 = 1_c_size_t+int(ng+i, c_size_t)+syw*int(j+ng, c_size_t)+szw*int(ng, c_size_t)

        ! bases that depend only on (i,j); fold z-halo
        baseU_ij_w = 1_c_size_t+int(ng+i, c_size_t)+syu*int(j+ng, c_size_t)
        baseV_ij_w = 1_c_size_t+int(ng+i, c_size_t)+syv*int(j+ng, c_size_t)

        c = c0
        !$omp simd linear(c:szw)
        do k = 0, nzw-1
          ip = c+sxw; im = c-sxw; ipp = c+2*sxw; imm = c-2*sxw
          jp = c+syw; jm = c-syw; jpp = c+2*syw; jmm = c-2*syw
          kp = c+szw; km = c-szw; kpp = c+2*szw; kmm = c-2*szw

          wc = w(c)

          iuL = baseU_ij_w+szu*int(k+ng, c_size_t)
          iuU = baseU_ij_w+szu*int(k+ng+1, c_size_t)
          ivL2 = baseV_ij_w+szv*int(k+ng, c_size_t)
          ivU2 = baseV_ij_w+szv*int(k+ng+1, c_size_t)

          uc = 0.5d0*(u(iuL)+u(iuU))
          vc = 0.5d0*(v(ivL2)+v(ivU2))

          Fx = kk_flux_hybrid_kk3(w(imm), w(im), w(c), w(ip), w(ipp), uc, dx)
          Fy = kk_flux_hybrid_kk3(w(jmm), w(jm), w(c), w(jp), w(jpp), vc, dy)
          Fz = kk_flux_hybrid_kk3(w(kmm), w(km), w(c), w(kp), w(kpp), wc, dz)

          Nw(c) = -(Fx+Fy+Fz)
          c = c+szw
        end do
      end do
    end do
  end subroutine advect_kk3

end module kernels
