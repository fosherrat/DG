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
        call compute_time_step(dt, dt_limit)

        call rhs_conv()
        if (eq.eq.21 .or. shock_capture.eq.1) call rhs_diff()

        do i=1,nvol
            do j=1,ndof
                cons(:,j,i) = cons(:,j,i) + dt*conv(:,j,i)/(det_jcb(i)*mass(j,j))
                if (eq.eq.21 .or. shock_capture.eq.1) then
                    cons(:,j,i) = cons(:,j,i) + dt*diff(:,j,i)/(det_jcb(i)*mass(j,j))
                endif
            enddo
        enddo
    case default
        write(*,*) 'ERROR: euler_explicit only supports dim = 1'
        stop
    end select

end subroutine euler_explicit





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
            speed = wave_speed_1d(u)
            max_speed = max(max_speed,speed)
            rho_min = min(rho_min,u(1))
        enddo

        if (max_speed.le.0.d0) then
            write(*,*) 'ERROR: non-positive wave speed in compute_time_step'
            stop
        endif

        dt = min(dt,cfl*vol(elem)/(dble(2*ndof-1)*max_speed))

        alpha_diff = 0.d0
        if (eq.eq.21) then
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

    double precision function poly_value(coef,ncoef,x)
        integer, intent(in) :: ncoef
        double precision, intent(in) :: coef(ncoef), x

        integer :: ii

        poly_value = 0.d0
        do ii=ncoef,1,-1
            poly_value = poly_value*x + coef(ii)
        enddo
    end function poly_value

end subroutine compute_time_step





subroutine implicit_condition_number(dt, cond_number, matrix_size, info)
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    double precision, intent(out) :: dt, cond_number
    integer, intent(out) :: matrix_size, info

    integer :: col, n
    double precision :: eps_fd, norm_a, norm_ainv
    double precision, allocatable :: u_base(:), u_pert(:)
    double precision, allocatable :: rhs_base(:), rhs_pert(:)
    double precision, allocatable :: amat(:,:), ainv(:,:)

    if (dim.ne.1) then
        write(*,*) 'ERROR: implicit_condition_number only supports dim = 1'
        stop
    endif

    call compute_time_step(dt, t_final)
    if (dt.le.0.d0) then
        info = -1
        cond_number = huge(1.d0)
        matrix_size = 0
        return
    endif

    n = nvar*ndof*nvol
    matrix_size = n
    allocate(u_base(n), u_pert(n), rhs_base(n), rhs_pert(n))
    allocate(amat(n,n), ainv(n,n))

    call pack_solution(u_base)
    call residual_from_vector(u_base, rhs_base)

    amat = 0.d0
    do col=1,n
        u_pert(:) = u_base(:)
        eps_fd = sqrt(epsilon(1.d0))*max(1.d0,abs(u_base(col)))
        u_pert(col) = u_pert(col) + eps_fd

        call residual_from_vector(u_pert, rhs_pert)
        amat(:,col) = -dt*(rhs_pert(:)-rhs_base(:))/eps_fd
        amat(col,col) = amat(col,col) + 1.d0
    enddo

    call unpack_solution(u_base)

    call invert_matrix(amat, ainv, info)
    if (info.eq.0) then
        norm_a = matrix_norm1(amat)
        norm_ainv = matrix_norm1(ainv)
        cond_number = norm_a*norm_ainv
    else
        cond_number = huge(1.d0)
    endif

    deallocate(u_base, u_pert, rhs_base, rhs_pert, amat, ainv)

contains

    integer function vector_index(var, dof, elem)
        integer, intent(in) :: var, dof, elem

        vector_index = ((elem-1)*ndof + (dof-1))*nvar + var
    end function vector_index

    subroutine pack_solution(vec)
        double precision, intent(out) :: vec(:)

        integer :: elem, dof, var, idx

        do elem=1,nvol
            do dof=1,ndof
                do var=1,nvar
                    idx = vector_index(var,dof,elem)
                    vec(idx) = cons(var,dof,elem)
                enddo
            enddo
        enddo
    end subroutine pack_solution

    subroutine unpack_solution(vec)
        double precision, intent(in) :: vec(:)

        integer :: elem, dof, var, idx

        do elem=1,nvol
            do dof=1,ndof
                do var=1,nvar
                    idx = vector_index(var,dof,elem)
                    cons(var,dof,elem) = vec(idx)
                enddo
            enddo
        enddo
    end subroutine unpack_solution

    subroutine residual_from_vector(vec, rhs)
        double precision, intent(in) :: vec(:)
        double precision, intent(out) :: rhs(:)

        integer :: elem, dof, var, idx

        call unpack_solution(vec)
        call rhs_conv()
        if (eq.eq.21 .or. shock_capture.eq.1) call rhs_diff()

        do elem=1,nvol
            do dof=1,ndof
                do var=1,nvar
                    idx = vector_index(var,dof,elem)
                    rhs(idx) = conv(var,dof,elem)/(det_jcb(elem)*mass(dof,dof))
                    if (eq.eq.21 .or. shock_capture.eq.1) then
                        rhs(idx) = rhs(idx) &
                            + diff(var,dof,elem)/(det_jcb(elem)*mass(dof,dof))
                    endif
                enddo
            enddo
        enddo
    end subroutine residual_from_vector

    double precision function matrix_norm1(a)
        double precision, intent(in) :: a(:,:)

        integer :: j
        double precision :: col_sum

        matrix_norm1 = 0.d0
        do j=1,size(a,2)
            col_sum = sum(abs(a(:,j)))
            matrix_norm1 = max(matrix_norm1,col_sum)
        enddo
    end function matrix_norm1

    subroutine invert_matrix(a, ainv, ierr)
        double precision, intent(in) :: a(:,:)
        double precision, intent(out) :: ainv(:,:)
        integer, intent(out) :: ierr

        integer :: i, j, nloc, pivot
        double precision :: factor, pivot_abs, tmp
        double precision, allocatable :: aug(:,:)

        nloc = size(a,1)
        allocate(aug(nloc,2*nloc))

        aug(:,1:nloc) = a(:,:)
        aug(:,nloc+1:2*nloc) = 0.d0
        do i=1,nloc
            aug(i,nloc+i) = 1.d0
        enddo

        ierr = 0
        do i=1,nloc
            pivot = i
            pivot_abs = abs(aug(i,i))
            do j=i+1,nloc
                if (abs(aug(j,i)).gt.pivot_abs) then
                    pivot = j
                    pivot_abs = abs(aug(j,i))
                endif
            enddo

            if (pivot_abs.le.tiny(1.d0)) then
                ierr = i
                deallocate(aug)
                return
            endif

            if (pivot.ne.i) then
                do j=1,2*nloc
                    tmp = aug(i,j)
                    aug(i,j) = aug(pivot,j)
                    aug(pivot,j) = tmp
                enddo
            endif

            factor = aug(i,i)
            aug(i,:) = aug(i,:)/factor

            do j=1,nloc
                if (j.ne.i) then
                    factor = aug(j,i)
                    aug(j,:) = aug(j,:) - factor*aug(i,:)
                endif
            enddo
        enddo

        ainv(:,:) = aug(:,nloc+1:2*nloc)
        deallocate(aug)
    end subroutine invert_matrix

end subroutine implicit_condition_number
