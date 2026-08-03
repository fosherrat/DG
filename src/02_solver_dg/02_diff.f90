subroutine rhs_diff()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: i
    double precision :: uhat(nvar, nsur_tot), u_p(nvar), u_m(nvar)
    double precision :: phi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BR1(implementation)
    ! step1: calculate 0.5*(u+ + u-) for every interface
    uhat = 0.d0
    do i=1,nsur_tot
        if(con_sur_vol(2,i).eq.-1) then
            cycle !temporary boundary
        endif
        u_p = 0.d0; u_m = 0.d0
        do j=1,nvar
            do k=1,ndof
                call poly_1d(legendre(:,k),sur_cen(1,i),phi)
                u_p(j) = u_p(j) + cons(j,k, 1)*phi
                u_m(j) = u_m(j) + cons(j,k,-1)*phi
            enddo
        enddo
        uhat(:,i) = 0.5d0*(u_p(:) + u_m(:))   
    enddo

    ! step2: calculate gradient for all element
    allocate(grad_cons(dim,nvar,ndof,nvol)); grad_cons = 0.d0
    do i=1,nvol
        do j=1,dim
            do k=1,nvar
                do m=1,ndof
                    
                enddo
            enddo
        enddo
    enddo

    ! step3: calculate volume viscous flux with information from step2
    do i=1,nvol
    enddo

    ! step4: calculate interface viscous flux
    do i=1,nsur_tot
    enddo
end subroutine rhs_diff