subroutine read_input()
    use m_folder
    use m_input
    use m_point_dg

    implicit none

    integer :: ios
    character(len=512) :: line

    open(1,file=''//folder_input//'input02.dat',status='old')
    read(1,*); read(1,*); read(1,*); read(1,*)
    read(1,'(A)') line
    close(1)

    viscosity = 0.d0
    prandtl = 0.72d0
    read(line,*,iostat=ios) order, flux, time, cfl, t_final, eq, diffusion, viscosity, prandtl

    if (ios.ne.0) then
        read(line,*,iostat=ios) order, flux, time, cfl, t_final, eq, diffusion
        if (ios.ne.0) then
            write(*,*) 'ERROR: invalid data line in input02.dat'
            stop
        endif
        if (eq.eq.21) then
            write(*,*) 'ERROR: eq = 21 requires viscosity and Prandtl number'
            stop
        endif
    endif

    if (order.lt.0) then
        write(*,*) 'ERROR: order must be non-negative. order = ', order
        stop
    endif
    if (t_final.lt.0.d0) then
        write(*,*) 'ERROR: t_final must be non-negative. t_final = ', t_final
        stop
    endif
    if (prandtl.le.0.d0) then
        write(*,*) 'ERROR: Prandtl number must be positive. Pr = ', prandtl
        stop
    endif
    if (diffusion.ne.0 .and. diffusion.ne.1) then
        write(*,*) 'ERROR: diffusion must be 0 (BR1) or 1 (BR2). diffusion = ', diffusion
        stop
    endif

    select case(eq)
    case(20)
        viscosity = 0.d0
    case(21)
        if (viscosity.le.0.d0) then
            write(*,*) 'ERROR: viscosity must be positive for eq = 21'
            stop
        endif
    case default
        write(*,*) 'ERROR: current solver supports only eq = 20 or 21. eq = ', eq
        stop
    end select

    select case(dim)
    case(1)
        ndof = order + 1
        nvar = 3
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
    allocate(jcb(dim,dim,nvol), det_jcb(nvol), inv_jcb(dim,dim,nvol))
    allocate(sur_id(nsur_tot), sur_nv(nsur_tot), sur(nsur_tot))
    allocate(sur_cen(dim,nsur_tot), sur_vec(dim,nsur_tot), con_sur_vol(2,nsur_tot))

    do i=1,nver
        read(2,*) id, x_ver(:,i)
    enddo

    do i=1,nvol
        read(2,*) vol_id(i), vol_nv(i), vol(i), vol_cen(:,i), vol_con(:,i), &
            det_jcb(i), jcb(:,:,i), inv_jcb(:,:,i)
    enddo

    do i=1,nsur_tot
        read(2,*) sur_id(i), sur_nv(i), sur(i), sur_cen(:,i), sur_vec(:,i), &
            con_sur_vol(:,i)
    enddo

    close(2)
end subroutine read_point
