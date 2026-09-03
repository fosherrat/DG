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






subroutine geo_surface(i,size1)
    use m_parameter
    use m_input
    use m_point_dg

    implicit none

    integer, intent(in) :: i, size1

    integer :: j, k
    double precision :: x(3,4), x1(3), x2(3), normal(3), normal2(3)
    double precision :: norm_normal, dir
    !==========================================================================================
    !calculate surface center
    do j=1,size1
        sur_cen(:,i) = sur_cen(:,i)+x_ver(:,sur_con(j,i))
    enddo
    sur_cen(:,i) = sur_cen(:,i)/dble(size1)

    x = 0.d0
    do j=1,size1
        x(1:dim,j) = x_ver(:,sur_con(j,i))
    enddo

    ! calculate surface normal vector(sur_vec) and surface area(sur)
    select case(dim)
    case(1)
        dir = 1.d0
        sur_vec(:,i) = (/1.d0/)
        k = con_sur_vol(1,i)
        if(dot_product(sur_vec(:,i),(sur_cen(:,i)-vol_cen(:,k))).lt.0) dir= -1.d0
        sur_vec(:,i) = dir*sur_vec(:,i)
        sur(i) = 1.d0
    case(2)
        if(size1.ne.2) then
            write(*,*) 'ERROR: Unsupported 2D surface vertex count: ', size1; stop
        endif

        x1 = x(:,2)-x(:,1)
        sur(i) = sqrt(dot_product(x1,x1))
        if(sur(i).le.0.d0) then
            write(*,*) 'ERROR: Degenerate 2D surface: ', i; stop
        endif

        ! Either perpendicular direction is valid initially; orient it from left to right.
        sur_vec(:,i) = (/ -x1(2), x1(1) /)/sur(i)
        k = con_sur_vol(1,i)
        if(dot_product(sur_vec(:,i),sur_cen(:,i)-vol_cen(:,k)).lt.0.d0) then
            sur_vec(:,i) = -sur_vec(:,i)
        endif
    case(3)
        select case(size1)
        case(3) ! triangular surface
            x1 = x(:,2)-x(:,1)
            x2 = x(:,3)-x(:,1)
            call cross(x1,x2,normal)
            norm_normal = sqrt(dot_product(normal,normal))
            sur(i) = 0.5d0*norm_normal
        case(4) ! quadrilateral: vertex IDs are sorted, so use all 4 triangles
            x1 = x(:,2)-x(:,1)
            x2 = x(:,3)-x(:,1)
            call cross(x1,x2,normal)
            norm_normal = sqrt(dot_product(normal,normal))
            sur(i) = norm_normal

            x1 = x(:,2)-x(:,1)
            x2 = x(:,4)-x(:,1)
            call cross(x1,x2,normal2)
            sur(i) = sur(i) + sqrt(dot_product(normal2,normal2))

            x1 = x(:,3)-x(:,1)
            x2 = x(:,4)-x(:,1)
            call cross(x1,x2,normal2)
            sur(i) = sur(i) + sqrt(dot_product(normal2,normal2))

            x1 = x(:,3)-x(:,2)
            x2 = x(:,4)-x(:,2)
            call cross(x1,x2,normal2)
            sur(i) = 0.25d0*(sur(i) + sqrt(dot_product(normal2,normal2)))
        case default
            write(*,*) 'ERROR: Unsupported 3D surface vertex count: ', size1; stop
        end select

        if(norm_normal.le.0.d0) then
            write(*,*) 'ERROR: Degenerate 3D surface: ', i; stop
        endif
        sur_vec(:,i) = normal/norm_normal
        k = con_sur_vol(1,i)
        if(dot_product(sur_vec(:,i),sur_cen(:,i)-vol_cen(:,k)).lt.0.d0) then
            sur_vec(:,i) = -sur_vec(:,i)
        endif

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
    double precision :: x(3,8), x1(3), x2(3), x3(3), x_cross(3)
    !==========================================================================================
    !save vertex information
    x = 0.d0
    do j=1,size1
        x(1:dim,j) = x_ver(:,vol_con(j,i))
    enddo

    !calculate element center
    do j=1,size1
        vol_cen(:,i) = vol_cen(:,i) + x(:,j)
    enddo
    vol_cen(:,i) = vol_cen(:,i)/dble(size1)

    !calculate volume
    select case(dim)
    case(1)
        vol(i) = abs(x_ver(1,vol_con(1,i))-x_ver(1,vol_con(2,i)))
    case(2)
        select case(size1)
        case(3) ! triangle
            x1 = x(:,2)-x(:,1)
            x2 = x(:,3)-x(:,1)
            call cross(x1,x2,x_cross)
            vol(i) = 0.5d0*sqrt(dot_product(x_cross,x_cross))
        case(4) ! quadrilateral: triangles 123 and 134
            x1 = x(:,2)-x(:,1)
            x2 = x(:,3)-x(:,1)
            call cross(x1,x2,x_cross)
            vol(i) = 0.5d0*sqrt(dot_product(x_cross,x_cross))

            x1 = x(:,3)-x(:,1)
            x2 = x(:,4)-x(:,1)
            call cross(x1,x2,x_cross)
            vol(i) = vol(i) + 0.5d0*sqrt(dot_product(x_cross,x_cross))
        case default
            write(*,*) 'ERROR: Unsupported 2D element vertex count: ', size1; stop
        end select
    case(3)
        select case(size1)
        case(4) ! tetrahedron
            x1 = x(:,2)-x(:,1); x2 = x(:,3)-x(:,1); x3 = x(:,4)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = abs(dot_product(x1,x_cross))/6.d0
        case(6) ! prism: tetrahedra 1234, 2345, and 3456
            x1 = x(:,2)-x(:,1); x2 = x(:,3)-x(:,1); x3 = x(:,4)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,3)-x(:,2); x2 = x(:,4)-x(:,2); x3 = x(:,5)-x(:,2)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,4)-x(:,3); x2 = x(:,5)-x(:,3); x3 = x(:,6)-x(:,3)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0
        case(8) ! hexahedron: six tetrahedra sharing diagonal 1--7
            x1 = x(:,2)-x(:,1); x2 = x(:,3)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,3)-x(:,1); x2 = x(:,4)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,4)-x(:,1); x2 = x(:,8)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,8)-x(:,1); x2 = x(:,5)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,5)-x(:,1); x2 = x(:,6)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0

            x1 = x(:,6)-x(:,1); x2 = x(:,2)-x(:,1); x3 = x(:,7)-x(:,1)
            call cross(x2,x3,x_cross)
            vol(i) = vol(i) + abs(dot_product(x1,x_cross))/6.d0
        case default
            write(*,*) 'ERROR: Unsupported 3D element vertex count: ', size1; stop
        end select

    case default
        write(*,*) 'ERROR: Unsupported dimension: ', dim; stop
    end select

end subroutine geo_volume
