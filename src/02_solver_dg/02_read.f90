subroutine read_input()
    use m_folder
    use m_input
    use m_point_dg

    implicit none

    integer :: ios
    character(len=512) :: line

    open(1,file=''//folder_input//'input02.dat',status='old')
    read(1,'(A)') line
    read(1,'(A)') line
    read(1,'(A)') line

    ! Equation and initial condition
    read(1,'(A)') line
    read(1,*,iostat=ios) eq, initial_condition
    if (ios.ne.0) goto 100

    ! Physical parameters
    read(1,'(A)') line
    read(1,'(A)') line
    read(1,*,iostat=ios) viscosity, prandtl
    if (ios.ne.0) goto 100

    ! Spatial discretization
    read(1,'(A)') line
    read(1,'(A)') line
    read(1,*,iostat=ios) order, flux, diffusion, shock_capture, av_c, av_kappa
    if (ios.ne.0) goto 100

    ! Time integration
    read(1,'(A)') line
    read(1,'(A)') line
    read(1,*,iostat=ios) time, dt_type, t_order, t_fixed, cfl, t_final
    if (ios.ne.0) goto 100

    ! Linear solver
    read(1,'(A)') line
    read(1,'(A)') line
    read(1,*,iostat=ios) matrix, pre
    if (ios.ne.0) goto 100
    close(1)
    goto 101

100 continue
    close(1)
    write(*,*) 'ERROR: invalid input02.dat format'
    stop

101 continue

    select case(dim)
    case(1)
        ndof = order + 1
        if (eq.eq.10 .or. eq.eq.11) then
            nvar = 1
        else
            nvar = 3
        endif
    case default
        write(*,*) 'ERROR: Unsupported dimension = ', dim
        stop
    end select

end subroutine read_input

subroutine read_point()
    use m_folder
    use m_point_dg

    implicit none

    integer :: i, id

    open(2,file=''//folder_point//trim(file_point),status='old')

    read(2,*) dim
    read(2,*) nver
    read(2,*) nsur_tot, nsur_int, nsur_b
    read(2,*) nvol, ntri, nqua, ntet, npri, nhex

    allocate(x_ver(dim,nver))
    allocate(vol_id(nvol), vol_nv(nvol), vol(nvol), vol_cen(dim,nvol), vol_con(8,nvol))
    allocate(sur_id(nsur_tot), sur_nv(nsur_tot), sur(nsur_tot))
    allocate(sur_cen(dim,nsur_tot), sur_vec(dim,nsur_tot), con_sur_vol(2,nsur_tot))

    do i=1,nver
        read(2,*) id, x_ver(:,i)
    enddo

    do i=1,nvol
        read(2,*) vol_id(i), vol_nv(i), vol(i), vol_cen(:,i), vol_con(:,i)
    enddo

    do i=1,nsur_tot
        read(2,*) sur_id(i), sur_nv(i), sur(i), sur_cen(:,i), sur_vec(:,i), con_sur_vol(:,i)
    enddo

    close(2)
end subroutine read_point
