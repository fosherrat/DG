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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(legendre)) call build_basis()
    if (.not.allocated(mass)) call build_mass()

    if (allocated(prim)) deallocate(prim)
    if (allocated(cons)) deallocate(cons)
    allocate(prim(nvar,ndof,nvol)); prim = 0.d0
    allocate(cons(nvar,ndof,nvol)); cons = 0.d0

    nquad_init = max(ndof+3, 4)
    allocate(gauss_x_init(nquad_init), gauss_w_init(nquad_init))
    call gauss_legendre(nquad_init, gauss_x_init, gauss_w_init)

    ! rho = 1 + 0.2*sin(8pi*x), u = 1, p = 1
    do i = 1, nvol
        do q = 1, nquad_init

            xi = gauss_x_init(q)
            x = vol_cen(1,i) + jcb(1,1,i)*xi

            prim_q(1) = 1.d0 + 0.2d0*sin(2.d0*pi*x)
            prim_q(2) = 1.d0
            prim_q(3) = 1.d0

            cons_q(1) = prim_q(1)
            cons_q(2) = prim_q(1)*prim_q(2)
            cons_q(3) = prim_q(3)/(gamma-1.d0) + 0.5d0*prim_q(1)*prim_q(2)**2

            do j = 1, ndof
                phi = poly_value(legendre(:,j), ndof, xi)
                prim(:,j,i) = prim(:,j,i) + gauss_w_init(q)*prim_q(:)*phi/mass(j,j)
                cons(:,j,i) = cons(:,j,i) + gauss_w_init(q)*cons_q(:)*phi/mass(j,j)
            enddo

        enddo
    enddo

    deallocate(gauss_x_init, gauss_w_init)

contains

    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        double precision, intent(out) :: x(n), w(n)

        integer :: ii, jj, m
        double precision :: z, z_old, p0, p1, p2, dp
        double precision, parameter :: eps = 1.d-14

        m = (n+1)/2
        do ii=1,m
            z = cos(pi*(dble(ii)-0.25d0)/(dble(n)+0.5d0))

            do
                p0 = 1.d0
                p1 = z

                if (n.eq.1) then
                    p2 = p1
                else
                    do jj=2,n
                        p2 = ((2.d0*dble(jj)-1.d0)*z*p1 - dble(jj-1)*p0)/dble(jj)
                        p0 = p1
                        p1 = p2
                    enddo
                endif

                dp = dble(n)*(z*p1-p0)/(z*z-1.d0)
                z_old = z
                z = z_old - p1/dp
                if (abs(z-z_old).le.eps) exit
            enddo

            x(ii) = -z
            x(n+1-ii) = z
            w(ii) = 2.d0/((1.d0-z*z)*dp*dp)
            w(n+1-ii) = w(ii)
        enddo
    end subroutine gauss_legendre

    double precision function poly_value(coef, ncoef, x)
        integer, intent(in) :: ncoef
        double precision, intent(in) :: coef(ncoef), x

        integer :: ii

        poly_value = 0.d0
        do ii=ncoef,1,-1
            poly_value = poly_value*x + coef(ii)
        enddo
    end function poly_value

end subroutine init_sin
