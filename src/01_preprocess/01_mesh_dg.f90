subroutine mesh_setting_dg
    use m_parameter
    use m_input
    use m_point_dg

    implicit none

    integer :: i, id, s
    integer :: temp1(4)
    !==========================================================================================
    
    !calculate volume/volume center
    allocate(vol_cen(dim,nvol),vol(nvol))
    allocate(jcb(dim,dim,nvol), det_jcb(nvol), inv_jcb(dim,dim,nvol))
    vol_cen = 0.d0; vol = 0.d0
    jcb = 0.d0; det_jcb = 0.d0; inv_jcb = 0.d0

    do i=1,nvol
        call geo_volume(i,vol_nv(i))
    enddo

    !calcualte surface info
    allocate(con_sur_vol(2,nsur_tot)); con_sur_vol = -1
    allocate(sur_cen(dim,nsur_tot)); sur_cen = 0.d0
    allocate(sur_vec(dim,nsur_tot)); sur_vec = 0.d0
    allocate(sur(nsur_tot)); sur = 0.d0
    s = 1

    ! sort boundary surface vertex id
    do i=1,nsur_b
        id = i + nsur_int
        temp1(:) = sur_con(:,id); call sort(temp1, sur_nv(id)); sur_con(:,(id)) = temp1(:)
    enddo

    ! interior surface --> sur_con, sur_nv/ total surface --> left,right volume mapping
    select case(dim)
    case(1)
        do i=1,nvol ! line 1/2
            temp1(:) = (/vol_con(1,i),0,0,0/); call match(i,temp1,1,s)
            temp1(:) = (/vol_con(2,i),0,0,0/); call match(i,temp1,1,s)
        enddo
    case(2)
        do i=1,nvol
            select case(vol_nv(i))
            case(3) ! triangle 12/23/31
                temp1(:) = (/vol_con(1,i),vol_con(2,i),0,0/); call match(i,temp1,2,s)
                temp1(:) = (/vol_con(2,i),vol_con(3,i),0,0/); call match(i,temp1,2,s)
                temp1(:) = (/vol_con(3,i),vol_con(1,i),0,0/); call match(i,temp1,2,s)
            case(4) ! quad 12/23/34/41
                temp1(:) = (/vol_con(1,i),vol_con(2,i),0,0/); call match(i,temp1,2,s)
                temp1(:) = (/vol_con(2,i),vol_con(3,i),0,0/); call match(i,temp1,2,s)
                temp1(:) = (/vol_con(3,i),vol_con(4,i),0,0/); call match(i,temp1,2,s)
                temp1(:) = (/vol_con(4,i),vol_con(1,i),0,0/); call match(i,temp1,2,s)
            end select
        enddo
    case(3)
        do i=1,nvol
            select case(vol_nv(i))
            case(4) !tetrahedron 123/124/234/314
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(3,i),0/); call match(i,temp1,3,s)
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(4,i),0/); call match(i,temp1,3,s)
                temp1(:) = (/vol_con(2,i),vol_con(3,i),vol_con(4,i),0/); call match(i,temp1,3,s)
                temp1(:) = (/vol_con(3,i),vol_con(1,i),vol_con(4,i),0/); call match(i,temp1,3,s)
            case(8) !hexahedron 1234/1265/2376/3487/4158/5678
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(3,i),vol_con(4,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(6,i),vol_con(5,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(2,i),vol_con(3,i),vol_con(7,i),vol_con(6,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(3,i),vol_con(4,i),vol_con(8,i),vol_con(7,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(4,i),vol_con(1,i),vol_con(5,i),vol_con(8,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(5,i),vol_con(6,i),vol_con(7,i),vol_con(8,i)/); call match(i,temp1,4,s)
            case(6) !prism 123/1254/2365/3146/456
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(3,i),0/); call match(i,temp1,3,s)
                temp1(:) = (/vol_con(1,i),vol_con(2,i),vol_con(5,i),vol_con(4,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(2,i),vol_con(3,i),vol_con(6,i),vol_con(5,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(3,i),vol_con(1,i),vol_con(4,i),vol_con(6,i)/); call match(i,temp1,4,s)
                temp1(:) = (/vol_con(4,i),vol_con(5,i),vol_con(6,i),0/); call match(i,temp1,3,s)
            end select
        enddo
    end select

    ! total surface --> sur_cen, sur_vec, sur
    do i=1,nsur_tot
        call geo_surface(i,sur_nv(i))
    enddo

end subroutine mesh_setting_dg





subroutine match(i,temp1,size,s)
    use m_point_dg
    implicit none
    
    integer, intent(inout) :: s, temp1(4)
    integer, intent(in) :: i, size
    integer :: j
    logical :: matched
    !==========================================================================================
    call sort(temp1,size)
    
    matched = .false.
    do j=1,nsur_tot
        if(all(sur_con(:,j).eq.temp1)) then
            matched = .true.
            if(j.gt.nsur_int) then 
                con_sur_vol(1,j) = i
                exit
            else
                con_sur_vol(2,j) = i
                exit
            endif
        end if
    enddo

    if(matched.eqv..false.) then
        sur_con(:,s) = temp1
        con_sur_vol(1,s) = i
        sur_nv(s) = size
        s = s+1
    endif
end subroutine match





subroutine sort(temp1,size)
    implicit none

    integer, intent(inout) :: temp1(4)
    integer, intent(in) :: size
    integer :: i, j, value
    !==========================================================================================
    
     do i =2,size
        value = temp1(i)
        j = i - 1
        do while (j.ge.1)
            if (temp1(j) <= value) exit
            temp1(j + 1) = temp1(j); j = j - 1
        end do
        temp1(j + 1) = value
    end do
end subroutine sort





subroutine geo_surface(i,size1)
    use m_parameter
    use m_input
    use m_point_dg

    implicit none

    integer, intent(in) :: i, size1

    integer :: j, k
    double precision :: x1(dim), x2(dim), dir
    !==========================================================================================
    !calculate surface center
    do j=1,size1
        sur_cen(:,i) = sur_cen(:,i)+x_ver(:,sur_con(j,i))
    enddo
    sur_cen(:,i) = sur_cen(:,i)/dble(size1)

    ! calculate surface normal vector(sur_vec) and surface area(sur)
    select case(dim)
    case(1)
        dir = 1.d0
        sur_vec(:,i) = (/1.d0/)
        k = con_sur_vol(1,i)
        if(dot_product(sur_vec(:,i),(sur_cen(:,i)-vol_cen(:,k))).lt.0) dir= -1.d0
        sur_vec(:,i) = dir*sur_vec(:,i)
        sur(i) = 1.d0
    case default
        write(*,*) 'ERROR: Unsupported dimension: ', dim; stop
    end select
end subroutine geo_surface





subroutine geo_volume(i,size1)
    use m_parameter
    use m_input
    use m_point_dg

    implicit none

    integer, intent(in) :: i, size1
    integer :: j
    double precision :: x1, x2
    !==========================================================================================
    !calculate element center
    do j=1,size1
        vol_cen(:,i) = vol_cen(:,i)+x_ver(:,vol_con(j,i))
    enddo
    vol_cen(:,i) = vol_cen(:,i)/dble(size1)

    !calculate volume
    select case(dim)
    case(1)
        vol(i) = abs(x_ver(1,vol_con(1,i))-x_ver(1,vol_con(2,i)))
    case default
        write(*,*) 'ERROR: Unsupported dimension: ', dim; stop
    end select

    !calculate jacobian, inverse jacobian
    select case(dim)
    case(1)
        x1 = x_ver(1,vol_con(1,i))
        x2 = x_ver(1,vol_con(2,i))

        jcb(1,1,i) = 0.5d0 * (x2 - x1)
        det_jcb(i) = abs(jcb(1,1,i))
        inv_jcb(1,1,i) = 1.d0 / jcb(1,1,i)
    case default
    end select

end subroutine geo_volume
