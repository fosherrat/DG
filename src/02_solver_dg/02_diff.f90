subroutine rhs_diff()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: i, k, q, vl, vr
    double precision :: xi, phi, grad_phi, normal, eps_q, eps_face
    double precision :: uq(nvar), grad_uq(nvar), uhat(nvar,nsur_tot)
    double precision :: ul(nvar), ur(nvar), ql(nvar), qr(nvar)
    double precision :: jump_u(nvar), lift_l(nvar), lift_r(nvar)
    double precision :: fvl(nvar), fvr(nvar), fvhat(nvar)

    ! One-dimensional BR1/BR2 discretization for the physical and artificial
    ! viscous fluxes. Persson artificial viscosity uses the same lifted gradient.

    if (dim.ne.1) then
        write(*,*) 'ERROR: rhs_diff currently supports only dim = 1'
        stop
    endif
    if (shock_capture.eq.1 .and. .not.allocated(av_face)) then
        write(*,*) 'ERROR: artificial viscosity was not initialized'
        stop
    endif

    if (allocated(diff)) deallocate(diff)
    allocate(diff(nvar,ndof,nvol)); diff = 0.d0

    if (allocated(grad_cons)) deallocate(grad_cons)
    allocate(grad_cons(dim,nvar,ndof,nvol)); grad_cons = 0.d0

    !------------------------------------------------------------------------------------------
    ! Step 1: centered trace of the conservative variables.
    ! Boundary faces currently use the interior state, which corresponds to a natural
    ! zero-jump treatment. Physical viscous boundary conditions can replace this state later.
    uhat = 0.d0
    do i=1,nsur_tot
        vl = con_sur_vol(1,i)
        vr = con_sur_vol(2,i)

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        ul = 0.d0
        do k=1,ndof
            phi = poly_value(legendre(:,k),ndof,xi)
            ul(:) = ul(:) + cons(:,k,vl)*phi
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            ur = 0.d0
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                ur(:) = ur(:) + cons(:,k,vr)*phi
            enddo
        else
            ur(:) = ul(:)
        endif

        uhat(:,i) = 0.5d0*(ul(:)+ur(:))
    enddo

    !------------------------------------------------------------------------------------------
    ! Step 2a: volume part of
    ! int_E phi_k S dx = int_dE phi_k uhat*n ds - int_E d(phi_k)/dx U dx.
    do i=1,nvol
        do q=1,nquad
            xi = gauss_x(q)
            uq = 0.d0

            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                uq(:) = uq(:) + cons(:,k,i)*phi
            enddo

            do k=1,ndof
                grad_phi = poly_value(grad_legendre(:,k),ndof-1,xi)
                grad_cons(1,:,k,i) = grad_cons(1,:,k,i) &
                    - gauss_w(q)*uq(:)*grad_phi
            enddo
        enddo
    enddo

    ! Step 2b: interface part. sur_vec is the outward normal of the left element.
    do i=1,nsur_tot
        vl = con_sur_vol(1,i)
        vr = con_sur_vol(2,i)
        normal = sur_vec(1,i)

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        do k=1,ndof
            phi = poly_value(legendre(:,k),ndof,xi)
            grad_cons(1,:,k,vl) = grad_cons(1,:,k,vl) &
                + sur(i)*phi*uhat(:,i)*normal
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                grad_cons(1,:,k,vr) = grad_cons(1,:,k,vr) &
                    - sur(i)*phi*uhat(:,i)*normal
            enddo
        endif
    enddo

    do i=1,nvol
        do k=1,ndof
            grad_cons(1,:,k,i) = grad_cons(1,:,k,i)/(det_jcb(i)*mass(k,k))
        enddo
    enddo

    !------------------------------------------------------------------------------------------
    ! Step 3: viscous volume contribution, -int_E d(phi_k)/dx Fv dx.
    do i=1,nvol
        do q=1,nquad
            xi = gauss_x(q)
            uq = 0.d0
            grad_uq = 0.d0

            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                uq(:) = uq(:) + cons(:,k,i)*phi
                grad_uq(:) = grad_uq(:) + grad_cons(1,:,k,i)*phi
            enddo

            eps_q = 0.d0
            if (shock_capture.eq.1) then
                eps_q = 0.5d0*((1.d0-xi)*av_end(1,i) + (1.d0+xi)*av_end(2,i))
            endif
            call viscous_flux_1d(uq,grad_uq,eps_q,fvl)

            do k=1,ndof
                grad_phi = poly_value(grad_legendre(:,k),ndof-1,xi)
                diff(:,k,i) = diff(:,k,i) - gauss_w(q)*fvl(:)*grad_phi
            enddo
        enddo
    enddo

    !------------------------------------------------------------------------------------------
    ! Step 4: viscous interface flux.
    ! BR1:  Fvhat = average(Fv(U, grad(U) + global lift)).
    ! BR2:  Fvhat = average(Fv(U, grad(U) + lift local to this face)).
    do i=1,nsur_tot
        vl = con_sur_vol(1,i)
        vr = con_sur_vol(2,i)
        normal = sur_vec(1,i)

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        ul = 0.d0
        ql = 0.d0
        do k=1,ndof
            phi = poly_value(legendre(:,k),ndof,xi)
            ul(:) = ul(:) + cons(:,k,vl)*phi

            select case(diffusion)
            case(0)
                ql(:) = ql(:) + grad_cons(1,:,k,vl)*phi
            case(1)
                grad_phi = poly_value(grad_legendre(:,k),ndof-1,xi) &
                    *inv_jcb(1,1,vl)
                ql(:) = ql(:) + cons(:,k,vl)*grad_phi
            end select
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            ur = 0.d0
            qr = 0.d0
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                ur(:) = ur(:) + cons(:,k,vr)*phi

                select case(diffusion)
                case(0)
                    qr(:) = qr(:) + grad_cons(1,:,k,vr)*phi
                case(1)
                    grad_phi = poly_value(grad_legendre(:,k),ndof-1,xi) &
                        *inv_jcb(1,1,vr)
                    qr(:) = qr(:) + cons(:,k,vr)*grad_phi
                end select
            enddo
        else
            ur(:) = ul(:)
            qr(:) = ql(:)
        endif

        if (diffusion.eq.1 .and. vr.gt.0) then
            jump_u(:) = (ul(:)-ur(:))*normal
            lift_l = 0.d0
            lift_r = 0.d0

            xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                lift_l(:) = lift_l(:) - 0.5d0*sur(i)*phi*phi*jump_u(:) &
                    /(det_jcb(vl)*mass(k,k))
            enddo

            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                lift_r(:) = lift_r(:) - 0.5d0*sur(i)*phi*phi*jump_u(:) &
                    /(det_jcb(vr)*mass(k,k))
            enddo

            ql(:) = ql(:) + lift_l(:)
            qr(:) = qr(:) + lift_r(:)
        endif

        eps_face = 0.d0
        if (shock_capture.eq.1) eps_face = av_face(i)
        call viscous_flux_1d(ul,ql,eps_face,fvl)
        call viscous_flux_1d(ur,qr,eps_face,fvr)
        fvhat(:) = 0.5d0*(fvl(:)+fvr(:))

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        do k=1,ndof
            phi = poly_value(legendre(:,k),ndof,xi)
            diff(:,k,vl) = diff(:,k,vl) + sur(i)*fvhat(:)*normal*phi
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            do k=1,ndof
                phi = poly_value(legendre(:,k),ndof,xi)
                diff(:,k,vr) = diff(:,k,vr) - sur(i)*fvhat(:)*normal*phi
            enddo
        endif
    enddo

contains

    subroutine viscous_flux_1d(u,grad_u,eps_av,fv)
        double precision, intent(in) :: u(nvar), grad_u(nvar), eps_av
        double precision, intent(out) :: fv(nvar)

        double precision :: rho, vel, etot, grad_vel, grad_eint, tau

        ! Persson artificial viscosity applies a scalar Laplacian to every
        ! conservative equation: F_av = eps_av*grad(U).
        fv(:) = eps_av*grad_u(:)
        if (eq.ne.21) return

        rho = u(1)
        if (rho.le.0.d0) then
            write(*,*) 'ERROR: non-positive density in viscous_flux_1d'
            stop
        endif

        vel = u(2)/rho
        etot = u(3)/rho
        grad_vel = (grad_u(2)-vel*grad_u(1))/rho
        grad_eint = (grad_u(3)-etot*grad_u(1))/rho - vel*grad_vel

        ! Stokes hypothesis in a one-dimensional flow: tau_xx = 4/3*mu*du/dx.
        tau = (4.d0/3.d0)*viscosity*grad_vel

        fv(2) = fv(2) + tau
        fv(3) = fv(3) + vel*tau + viscosity*gamma*grad_eint/prandtl
    end subroutine viscous_flux_1d

    double precision function poly_value(coef,ncoef,x)
        integer, intent(in) :: ncoef
        double precision, intent(in) :: coef(ncoef), x

        integer :: ii

        poly_value = 0.d0
        do ii=ncoef,1,-1
            poly_value = poly_value*x + coef(ii)
        enddo
    end function poly_value

end subroutine rhs_diff
