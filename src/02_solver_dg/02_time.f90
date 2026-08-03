subroutine euler_explicit(dt, dt_limit)
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    double precision, intent(out) :: dt
    double precision, intent(in) :: dt_limit

    integer :: i, j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case(dim)
    case(1)
        call compute_dt(dt)
        dt = min(dt, dt_limit)

        call rhs_conv()
        if (eq.eq.21) call rhs_diff()

        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = cons(:,j,i) + dt*conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (eq.eq.21) then
                    cons(:,j,i) = cons(:,j,i) + dt*diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo
    case default
        write(*,*) 'ERROR: euler_explicit only supports dim = 1'
        stop
    end select

contains

    subroutine compute_dt(dt)
        double precision, intent(out) :: dt

        integer :: elem, q, npoint
        double precision :: xi, speed, max_speed, rho_min, alpha_diff, dt_diff
        double precision :: u(nvar)

        if (cfl.le.0.d0) then
            write(*,*) 'ERROR: CFL must be positive'
            stop
        endif

        dt = huge(1.d0)
        npoint = max(ndof+1,3)

        do elem=1,nvol
            max_speed = 0.d0
            rho_min = huge(1.d0)

            do q=1,npoint
                xi = -1.d0 + 2.d0*dble(q-1)/dble(npoint-1)
                call solution_at(elem,xi,u)
                speed = wave_speed_1d(u)
                max_speed = max(max_speed,speed)
                rho_min = min(rho_min,u(1))
            enddo

            if (max_speed.le.0.d0) then
                write(*,*) 'ERROR: non-positive wave speed in compute_dt'
                stop
            endif

            dt = min(dt,cfl*vol(elem)/(dble(2*ndof-1)*max_speed))

            if (eq.eq.21) then
                if (rho_min.le.0.d0) then
                    write(*,*) 'ERROR: non-positive density in compute_dt'
                    stop
                endif

                alpha_diff = max(viscosity/rho_min, &
                    gamma*viscosity/(prandtl*rho_min))
                dt_diff = cfl*vol(elem)**2 &
                    /(dble(2*ndof-1)**2*max(alpha_diff,tiny(1.d0)))
                dt = min(dt,dt_diff)
            endif
        enddo
    end subroutine compute_dt

    subroutine solution_at(elem,xi,u)
        integer, intent(in) :: elem
        double precision, intent(in) :: xi
        double precision, intent(out) :: u(nvar)

        integer :: k
        double precision :: phi

        u = 0.d0
        do k=1,ndof
            phi = poly_value(legendre(:,k),ndof,xi)
            u(:) = u(:) + cons(:,k,elem)*phi
        enddo
    end subroutine solution_at

    subroutine primitive_from_cons(u,v)
        double precision, intent(in) :: u(nvar)
        double precision, intent(out) :: v(nvar)

        double precision :: rho, vel, pressure

        rho = u(1)
        if (rho.le.0.d0) then
            write(*,*) 'ERROR: non-positive density in primitive_from_cons'
            stop
        endif

        vel = u(2)/rho
        pressure = (gamma-1.d0)*(u(3)-0.5d0*rho*vel*vel)
        if (pressure.le.0.d0) then
            write(*,*) 'ERROR: non-positive pressure in primitive_from_cons'
            stop
        endif

        v(1) = rho
        v(2) = vel
        v(3) = pressure
    end subroutine primitive_from_cons

    double precision function wave_speed_1d(u)
        double precision, intent(in) :: u(nvar)

        double precision :: v(nvar)

        call primitive_from_cons(u,v)
        wave_speed_1d = abs(v(2)) + sqrt(gamma*v(3)/v(1))
    end function wave_speed_1d

    double precision function poly_value(coef,ncoef,x)
        integer, intent(in) :: ncoef
        double precision, intent(in) :: coef(ncoef), x

        integer :: ii

        poly_value = 0.d0
        do ii=ncoef,1,-1
            poly_value = poly_value*x + coef(ii)
        enddo
    end function poly_value

end subroutine euler_explicit
