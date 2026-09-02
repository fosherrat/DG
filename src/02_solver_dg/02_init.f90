subroutine init()
    use m_input

    implicit none

    select case(initial_condition)
    case(0)
        call init_sin()
    case(1)
        call init_const()
    case(2)
        call init_sod()
    case default
        write(*,*) 'ERROR: unsupported initial condition = ', initial_condition
        stop
    end select
end subroutine init




subroutine init_sin()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: i, j, q, nquad_init
    double precision :: xi, x, phi
    double precision :: prim_q(nvar), cons_q(nvar)
    double precision, allocatable :: gauss_x_init(:), gauss_w_init(:)
    double precision, external :: poly_value
    external :: gauss_legendre_rule

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(prim(nvar,ndof,nvol)); prim = 0.d0
    allocate(cons(nvar,ndof,nvol)); cons = 0.d0

    nquad_init = max(ndof+3, 4)
    allocate(gauss_x_init(nquad_init), gauss_w_init(nquad_init))
    call gauss_legendre_rule(nquad_init, gauss_x_init, gauss_w_init)

    ! rho = 1 + 0.2*sin(8pi*x), u = 1, p = 1
    do i = 1, nvol
        do q = 1, nquad_init

            xi = gauss_x_init(q)
            x = vol_cen(1,i) + jcb(1,1,i)*xi

            prim_q(1) = 1.d0 + 0.2d0*sin(8.d0*pi*x)
            cons_q(1) = prim_q(1)
            if (nvar.eq.3) then
                prim_q(2) = 1.d0
                prim_q(3) = 1.d0
                cons_q(2) = prim_q(1)*prim_q(2)
                cons_q(3) = prim_q(3)/(gamma-1.d0) + 0.5d0*prim_q(1)*prim_q(2)**2
            endif

            do j = 1, ndof
                phi = poly_value(legendre(:,j), ndof, xi)
                prim(:,j,i) = prim(:,j,i) + gauss_w_init(q)*prim_q(:)*phi/mass(j,j)
                cons(:,j,i) = cons(:,j,i) + gauss_w_init(q)*cons_q(:)*phi/mass(j,j)
            enddo

        enddo
    enddo

    deallocate(gauss_x_init, gauss_w_init)

end subroutine init_sin




subroutine init_const()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: i
    double precision :: prim_const(nvar), cons_const(nvar)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(prim(nvar,ndof,nvol)); prim = 0.d0
    allocate(cons(nvar,ndof,nvol)); cons = 0.d0

    prim_const(1) = 1.d0
    cons_const(1) = prim_const(1)
    if (nvar.eq.3) then
        prim_const(2) = 1.d0
        prim_const(3) = 1.d0
        cons_const(2) = prim_const(1)*prim_const(2)
        cons_const(3) = prim_const(3)/(gamma-1.d0) &
            + 0.5d0*prim_const(1)*prim_const(2)**2
    endif

    do i=1,nvol
        prim(:,1,i) = prim_const(:)
        cons(:,1,i) = cons_const(:)
    enddo

end subroutine init_const




subroutine init_sod()
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: i
    double precision :: sod_x
    double precision :: prim_state(nvar), cons_state(nvar)

    ! Standard Sod shock tube: left (rho,u,p) = (1,0,1),
    ! right (rho,u,p) = (0.125,0,0.1).  Align the x = 0.5 diaphragm to
    allocate(prim(nvar,ndof,nvol)); prim = 0.d0
    allocate(cons(nvar,ndof,nvol)); cons = 0.d0

    sod_x = sur_cen(1,1)
    do i = 2, nsur_tot
        if (abs(sur_cen(1,i)-0.5d0).lt.abs(sod_x-0.5d0)) sod_x = sur_cen(1,i)
    enddo

    do i = 1, nvol
        if (vol_cen(1,i).lt.sod_x) then
            prim_state = (/1.d0, 0.d0, 1.d0/)
        else
            prim_state = (/0.125d0, 0.d0, 0.1d0/)
        endif

        cons_state(1) = prim_state(1)
        cons_state(2) = prim_state(1)*prim_state(2)
        cons_state(3) = prim_state(3)/(gamma-1.d0) + 0.5d0*prim_state(1)*prim_state(2)**2

        prim(:,1,i) = prim_state(:)
        cons(:,1,i) = cons_state(:)
    enddo

end subroutine init_sod
