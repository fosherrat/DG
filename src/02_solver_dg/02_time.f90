subroutine time_step(dt)
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    double precision, intent(in) :: dt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case(time)
    case(0)
        call time_explicit(dt)
    case(1)
        call time_implicit(dt)
    case default
        write(*,*) 'ERROR(subroutine time): Unsupproted time type t_type = ', time
    end select
end subroutine time_step





subroutine time_explicit(dt)
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    double precision, intent(in) :: dt

    integer :: i,j
    logical :: has_diff
    double precision, allocatable :: temp1(:,:,:), rk1(:,:,:), rk2(:,:,:), rk3(:,:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! t_order determines the RK order(later expanded to Adams-Bashforth method)

    has_diff = (eq.eq.11 .or. eq.eq.21 .or. shock_capture.eq.1)

    select case(t_order)
    case(1)
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = cons(:,j,i) + dt*conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    cons(:,j,i) = cons(:,j,i) + dt*diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo
    case(2)
        if (has_diff) call rhs_diff()
        call rhs_conv()
        allocate(temp1(nvar,ndof,nvol)); temp1 = cons
        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = temp1(:,j,i) + 0.5d0*dt*conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    cons(:,j,i) = cons(:,j,i) + 0.5d0*dt*diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = temp1(:,j,i) + dt*conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    cons(:,j,i) = cons(:,j,i) + dt*diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo
        deallocate(temp1)
    case(4)
        allocate(temp1(nvar,ndof,nvol), rk1(nvar,ndof,nvol), rk2(nvar,ndof,nvol), rk3(nvar,ndof,nvol))
        temp1 = cons

        ! k1 = L(U^n)
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                rk1(:,j,i) = conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    rk1(:,j,i) = rk1(:,j,i) + diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
                cons(:,j,i) = temp1(:,j,i) + 0.5d0*dt*rk1(:,j,i)
            enddo
        enddo

        ! k2 = L(U^n + dt*k1/2)
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                rk2(:,j,i) = conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    rk2(:,j,i) = rk2(:,j,i) + diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
                cons(:,j,i) = temp1(:,j,i) + 0.5d0*dt*rk2(:,j,i)
            enddo
        enddo

        ! k3 = L(U^n + dt*k2/2)
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                rk3(:,j,i) = conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (has_diff) then
                    rk3(:,j,i) = rk3(:,j,i) + diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
                cons(:,j,i) = temp1(:,j,i) + dt*rk3(:,j,i)
            enddo
        enddo

        ! k4 = L(U^n + dt*k3), retained in conv/diff
        if (has_diff) call rhs_diff()
        call rhs_conv()
        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = temp1(:,j,i) + dt*(rk1(:,j,i) &
                    + 2.d0*rk2(:,j,i) + 2.d0*rk3(:,j,i) &
                    + conv(:,j,i)/(det_jcb(i)*mass(j,j)))/6.d0
                if (has_diff) then
                    cons(:,j,i) = cons(:,j,i) + dt*diff(:,j,i)/(6.d0*det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo

        deallocate(temp1, rk1, rk2, rk3)
    case default
        write(*,*) 'ERROR(sub time_explicit): Unsupproted RK order = ', t_order
    end select
end subroutine time_explicit





subroutine time_implicit(dt)
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none
    
    double precision, intent(in) :: dt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine time_implicit





subroutine compute_time_step(dt, dt_limit)
    use m_parameter
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    double precision, intent(out) :: dt
    double precision, intent(in) :: dt_limit

    integer :: elem, q, npoint
    double precision :: xi, speed, max_speed, rho_min, alpha_diff, dt_diff
    double precision :: u(nvar)
    double precision, external :: poly_value

    if (cfl.le.0.d0) then
        write(*,*) 'ERROR: CFL must be positive'
        stop
    endif

    ! The sensor is updated once per timestep. Residual/Jacobian evaluations then
    ! use this frozen viscosity field, as in the weakly coupled Persson approach.
    if (shock_capture.eq.1) call update_artificial_viscosity()

    dt = huge(1.d0)
    npoint = max(ndof+1,3)

    do elem=1,nvol
        max_speed = 0.d0
        rho_min = huge(1.d0)

        do q=1,npoint
            xi = -1.d0 + 2.d0*dble(q-1)/dble(npoint-1)
            call solution_at(elem,xi,u)
            if (eq.eq.10 .or. eq.eq.11) then
                speed = abs(linear_advection_speed)
            else
                speed = wave_speed_1d(u)
            endif
            max_speed = max(max_speed,speed)
            rho_min = min(rho_min,u(1))
        enddo

        if (max_speed.le.0.d0) then
            write(*,*) 'ERROR: non-positive wave speed in compute_time_step'
            stop
        endif

        dt = min(dt,cfl*vol(elem)/(dble(2*ndof-1)*max_speed))

        alpha_diff = 0.d0
        if (eq.eq.11) then
            alpha_diff = viscosity
        elseif (eq.eq.21) then
            if (rho_min.le.0.d0) then
                write(*,*) 'ERROR: non-positive density in compute_time_step'
                stop
            endif

            alpha_diff = max(viscosity/rho_min, &
                gamma*viscosity/(prandtl*rho_min))
        endif
        if (shock_capture.eq.1) then
            alpha_diff = max(alpha_diff,maxval(av_end(:,elem)))
        endif

        if (alpha_diff.gt.tiny(1.d0)) then
            dt_diff = cfl*vol(elem)**2 &
                /(dble(2*ndof-1)**2*alpha_diff)
            dt = min(dt,dt_diff)
        endif
    enddo

    dt = min(dt, dt_limit)

contains

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

end subroutine compute_time_step
