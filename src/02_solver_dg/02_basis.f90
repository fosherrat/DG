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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! gauss-legendre quadrature can approximate to 2*nquad-1 order
    ! nquad = ndof + 1 --> approximate to 2*order+3
    nquad = ndof+1
    allocate(gauss_x(nquad), gauss_w(nquad)); gauss_x = 0.d0; gauss_w = 0.d0
    
    select case(nquad)
    case(1)
        gauss_x(1) = 0.d0; gauss_w(1) = 2.d0
    case(2)
        gauss_x(1) =  1.d0/sqrt(3.d0); gauss_w(1) = 1.d0
        gauss_x(2) = -1.d0/sqrt(3.d0); gauss_w(2) = 1.d0
    case(3)
        gauss_x(1) = 0.d0; gauss_w(1) = 8.d0/9.d0
        gauss_x(2) =  sqrt(3.d0/5.d0); gauss_w(2) = 5.d0/9.d0
        gauss_x(3) = -sqrt(3.d0/5.d0); gauss_w(3) = 5.d0/9.d0
    case(4)
        gauss_x(1) =  sqrt(3.d0/7.d0 - 2.d0/7.d0*sqrt(6.d0/5.d0)); gauss_w(1) = (18.d0 + sqrt(30.d0))/36.d0
        gauss_x(2) = -sqrt(3.d0/7.d0 - 2.d0/7.d0*sqrt(6.d0/5.d0)); gauss_w(2) = (18.d0 + sqrt(30.d0))/36.d0
        gauss_x(3) =  sqrt(3.d0/7.d0 + 2.d0/7.d0*sqrt(6.d0/5.d0)); gauss_w(3) = (18.d0 - sqrt(30.d0))/36.d0
        gauss_x(4) = -sqrt(3.d0/7.d0 + 2.d0/7.d0*sqrt(6.d0/5.d0)); gauss_w(4) = (18.d0 - sqrt(30.d0))/36.d0
    case(5)
        gauss_x(1) = 0.d0; gauss_w(1) = 128.d0/225.d0
        gauss_x(2) =  sqrt(5.d0 - 2.d0*sqrt(10.d0/7.d0))/3.d0; gauss_w(2) = (322.d0 + 13.d0*sqrt(70.d0))/900.d0
        gauss_x(3) = -sqrt(5.d0 - 2.d0*sqrt(10.d0/7.d0))/3.d0; gauss_w(3) = (322.d0 + 13.d0*sqrt(70.d0))/900.d0
        gauss_x(4) =  sqrt(5.d0 + 2.d0*sqrt(10.d0/7.d0))/3.d0; gauss_w(4) = (322.d0 - 13.d0*sqrt(70.d0))/900.d0
        gauss_x(5) = -sqrt(5.d0 + 2.d0*sqrt(10.d0/7.d0))/3.d0; gauss_w(5) = (322.d0 - 13.d0*sqrt(70.d0))/900.d0
    case default
        write(*,*) 'ERROR: gauss quadrature point should be in range of 1-5. current nquad = ', nquad; stop
    end select
    
end subroutine build_gauss_quadrature
