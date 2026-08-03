subroutine rhs_conv()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow
    
    implicit none

    integer :: i, j, k, q, vl, vr
    double precision :: xq, xi, phi, grad_phi, normal
    double precision :: uq(nvar), fq(nvar), ul(nvar), ur(nvar), fhat(nvar)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! this is only for dim 1
    if (allocated(conv)) deallocate(conv)
    allocate(conv(nvar,ndof,nvol)); conv = 0.d0
    
    ! Gaussian quadrature integral approximation
    do i=1,nvol
        do q=1,nquad
            xq = gauss_x(q)
            uq = 0.d0

            do k=1,ndof
                phi = 0.d0
                do j=1,ndof
                    phi = phi + legendre(j,k)*xq**dble(j-1)
                enddo
                uq(:) = uq(:) + cons(:,k,i)*phi
            enddo

            call flux_euler(1,uq,fq)

            do k=1,ndof
                grad_phi = 0.d0
                do j=1,ndof-1
                    grad_phi = grad_phi + grad_legendre(j,k)*xq**dble(j-1)
                enddo
                conv(:,k,i) = conv(:,k,i) &
                    + gauss_w(q)*fq(:)*grad_phi*det_jcb(i)*inv_jcb(1,1,i)
            enddo
        enddo
    enddo

    do i=1,nsur_tot
        vl = con_sur_vol(1,i)
        vr = con_sur_vol(2,i)
        normal = sur_vec(1,i)

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        ul(:) = 0.d0
        do k=1,ndof
            phi = 0.d0
            do j=1,ndof
                phi = phi + legendre(j,k)*xi**dble(j-1)
            enddo
            ul(:) = ul(:) + cons(:,k,vl)*phi
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            ur(:) = 0.d0
            do k=1,ndof
                phi = 0.d0
                do j=1,ndof
                    phi = phi + legendre(j,k)*xi**dble(j-1)
                enddo
                ur(:) = ur(:) + cons(:,k,vr)*phi
            enddo
        else
            ur(:) = ul(:)
        endif

        select case(flux)
        case(0)
            call flux_llf(ul,ur,normal,fhat)
        case default
            write(*,*) 'ERROR: Unsupported flux = ', flux
            stop
        end select

        xi = (sur_cen(1,i)-vol_cen(1,vl))/jcb(1,1,vl)
        do k=1,ndof
            phi = 0.d0
            do j=1,ndof
                phi = phi + legendre(j,k)*xi**dble(j-1)
            enddo
            conv(:,k,vl) = conv(:,k,vl) - sur(i)*fhat(:)*phi
        enddo

        if (vr.gt.0) then
            xi = (sur_cen(1,i)-vol_cen(1,vr))/jcb(1,1,vr)
            do k=1,ndof
                phi = 0.d0
                do j=1,ndof
                    phi = phi + legendre(j,k)*xi**dble(j-1)
                enddo
                conv(:,k,vr) = conv(:,k,vr) + sur(i)*fhat(:)*phi
            enddo
        endif
    enddo

end subroutine rhs_conv





subroutine flux_euler(d,u,f)
    use m_parameter
    use m_input

    implicit none

    integer, intent(in) :: d
    double precision, intent(in) :: u(nvar)
    double precision, intent(out) :: f(nvar)

    double precision :: rho, vel, p

    ! Current implementation is one-dimensional.
    if (d.ne.1 .or. nvar.ne.3) then
        write(*,*) 'ERROR: flux_euler currently supports only the 1D Euler equations'
        stop
    endif

    rho = u(1)
    vel = u(2)/rho
    p = (gamma-1.d0)*(u(3)-0.5d0*rho*vel*vel)

    f(1) = u(2)
    f(2) = u(2)*vel + p
    f(3) = vel*(u(3)+p)

end subroutine flux_euler





subroutine flux_llf(ul,ur,normal,fhat)
    use m_parameter
    use m_input

    implicit none

    double precision, intent(in) :: ul(nvar), ur(nvar), normal
    double precision, intent(out) :: fhat(nvar)

    double precision :: pl, pr, alpha
    double precision :: fl(nvar), fr(nvar)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pl = (gamma-1.d0)*(ul(3)-0.5d0*ul(2)**2/ul(1))
    fl(:) = 0.d0
    fl(1) = ul(2)
    fl(2) = ul(2)**2/ul(1) + pl
    fl(3) = ul(2)*(ul(3)+pl)/ul(1)

    pr = (gamma-1.d0)*(ur(3)-0.5d0*ur(2)**2/ur(1))
    fr(:) = 0.d0
    fr(1) = ur(2)
    fr(2) = ur(2)**2/ur(1) + pr
    fr(3) = ur(2)*(ur(3)+pr)/ur(1)

    alpha = max(abs(ul(2)/ul(1))+sqrt(gamma*pl/ul(1)), &
                abs(ur(2)/ur(1))+sqrt(gamma*pr/ur(1)))
    fhat(:) = 0.5d0*(fl(:)+fr(:))*normal - 0.5d0*alpha*(ur(:)-ul(:))
end subroutine flux_llf
