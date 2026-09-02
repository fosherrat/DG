subroutine build_basis()
    use m_input
    use m_flow

    implicit none
    
    integer :: i,j,n
    double precision :: temp1(ndof), temp2(ndof-1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! (ndof-1) order solution in each element (ndof, ndof): number of terms, number of basis
    ! 1D Legendre basis: P_0 = 1, P_1 = x, P_i = [(2*i-1)*x*P_(i-1) - (i-1)*P_(i-2)]/(i)
    if (allocated(legendre)) deallocate(legendre)
    allocate(legendre(ndof,ndof)); legendre = 0.d0

    do i=1,ndof
        if(i.eq.1) legendre(1,i) = 1.d0
        if(i.eq.2) legendre(2,i) = 1.d0
        if(i.ge.3) then
            n = i-1
            legendre(:,i) = -(dble(n-1)/dble(n))*legendre(:,i-2)
            temp1 = 0.d0
            do j=2,ndof
                temp1(j) = (dble(2*n-1)/dble(n))*legendre(j-1,i-1)
            enddo
            legendre(:,i) = legendre(:,i) + temp1
        endif
    enddo

    ! (ndof-1, ndof): number of terms(1 order less than normal legendre), number of basis
    if (allocated(grad_legendre)) deallocate(grad_legendre)
    allocate(grad_legendre(ndof-1,ndof)); grad_legendre = 0.d0

    do i=1,ndof
        if(i.eq.1) cycle
        temp2 = 0.d0 
        do j=2,i
            temp2(j-1) = legendre(j,i)*dble(j-1)
        enddo
        grad_legendre(:,i) = grad_legendre(:,i) + temp2(:)
    enddo
end subroutine build_basis





subroutine build_mass()
    use m_input
    use m_flow

    implicit none

    integer :: i,j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (allocated(mass)) deallocate(mass)
    allocate(mass(ndof, ndof)); mass = 0.d0

    do i=1,ndof
        do j =1,ndof
            if(i.eq.j) then
                mass(j,i) = 2.d0/dble(2*i-1)
            endif
        enddo
    enddo

end subroutine build_mass





subroutine build_gauss_quadrature()
    use m_input
    use m_flow

    implicit none

    external :: gauss_legendre_rule

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nquad = ndof+1
    if (nquad.lt.1) then
        write(*,*) 'ERROR: nquad must be positive. current nquad = ', nquad
        stop
    endif

    allocate(gauss_x(nquad), gauss_w(nquad))
    call gauss_legendre_rule(nquad, gauss_x, gauss_w)

end subroutine build_gauss_quadrature




subroutine gauss_legendre_rule(n, x, w)
    implicit none

    integer, intent(in) :: n
    double precision, intent(out) :: x(n), w(n)
    integer :: i, j, m, iter
    logical :: converged
    double precision :: z, dz, p_nm1, p_n, p_np1, dp_n
    double precision, parameter :: tolerance = 1.d-14
    integer, parameter :: max_iterations = 100

    ! An n-point Gauss-Legendre rule exactly integrates polynomials through
    ! degree 2*n-1.  Find roots of P_n with Newton iteration and use symmetry.
    if (n.lt.1) then
        write(*,*) 'ERROR: Gauss-Legendre rule requires n >= 1. n = ', n
        stop
    endif

    m = (n+1)/2
    do i=1,m
        z = cos(acos(-1.d0)*(dble(i)-0.25d0)/(dble(n)+0.5d0))
        converged = .false.

        do iter=1,max_iterations
            p_nm1 = 1.d0
            p_n = z
            do j=2,n
                p_np1 = ((2.d0*dble(j)-1.d0)*z*p_n - dble(j-1)*p_nm1)/dble(j)
                p_nm1 = p_n
                p_n = p_np1
            enddo

            dp_n = dble(n)*(z*p_n-p_nm1)/(z*z-1.d0)
            dz = p_n/dp_n
            z = z - dz
            if (abs(dz).le.tolerance) then
                converged = .true.
                exit
            endif
        enddo

        if (.not.converged) then
            write(*,*) 'ERROR: Gauss-Legendre Newton iteration did not converge for root ', i
            stop
        endif

        ! Re-evaluate P_n' at the converged root before forming its weight.
        p_nm1 = 1.d0
        p_n = z
        do j=2,n
            p_np1 = ((2.d0*dble(j)-1.d0)*z*p_n - dble(j-1)*p_nm1)/dble(j)
            p_nm1 = p_n
            p_n = p_np1
        enddo
        dp_n = dble(n)*(z*p_n-p_nm1)/(z*z-1.d0)

        x(i) = -z
        x(n+1-i) = z
        w(i) = 2.d0/((1.d0-z*z)*dp_n*dp_n)
        w(n+1-i) = w(i)
    enddo

end subroutine gauss_legendre_rule
