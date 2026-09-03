subroutine read_input()
    use m_folder
    use m_input

    implicit none

    open(1,file=''//folder_input//'input01.dat',status='old')
    read(1,*); read(1,*); read(1,*); read(1,*)
    read(1,*) dim, para, scale, mesh_type, discretization, file_mesh
    close(1)

end subroutine read_input





subroutine read_msh()
    use m_parameter
    use m_folder
    use m_input
    use m_point_dg

    implicit none

    character(len=256) :: ctemp1
    integer :: i, j
    integer :: nid
    double precision, allocatable :: temp1(:)
    integer :: nele, eid, etype, ntags, nnode
    integer, allocatable :: ele_type(:), ele_phys(:), ele_nnode(:), ele_conn(:,:)
    integer, allocatable :: tags(:), node_list(:)
    !====================================================================================================
    open(2,file=''//folder_mesh//trim(file_mesh),status='old',action='read')

    do while(.true.)
        read(2,'(a)') ctemp1; if(trim(ctemp1).eq.'$Nodes') exit
    enddo

    read(2,*) nver

    allocate(x_ver(dim,nver))

    allocate(temp1(dim))
    do i=1,nver
        read(2,*) nid,temp1(:)
        x_ver(:,nid) = temp1(:)
    enddo
    deallocate(temp1)

    do while(.true.)
        read(2,'(a)') ctemp1; if(trim(ctemp1).eq.'$Elements') exit
    enddo

    read(2,*) nele

    nsur = 0; nvol = 0
    ntri = 0; nqua = 0; ntet = 0; npri = 0; nhex = 0;
    allocate(ele_type(nele), ele_phys(nele), ele_nnode(nele))
    allocate(ele_conn(8,nele))

    select case(dim)
    case(1)
        do i=1,nele
            read(2,'(a)') ctemp1; read(ctemp1,*) eid, etype, ntags

            select case(etype)
            case(1)
                nnode = 2; nvol = nvol+1;
            case(15)
                nnode = 1; nsur = nsur+1;
            case default
                write(*,*) 'ERROR: Unsupported Gmsh element type = ', etype; stop
            end select
            allocate(tags(ntags)); allocate(node_list(nnode))
            read(ctemp1,*) eid, etype, ntags, (tags(j),j=1,ntags), (node_list(j),j=1,nnode)
            ele_type(i) = etype; ele_phys(i) = tags(1); ele_nnode(i) = nnode
            ele_conn(1:nnode,i) = node_list
            deallocate(tags, node_list)
        enddo
        !calculate number of surface
        nsur_b   = nsur
        nsur_tot = nvol + 1
        nsur_int = nsur_tot - nsur_b
    case(2)
        do i=1,nele
            read(2,'(a)') ctemp1; read(ctemp1,*) eid, etype, ntags

            select case(etype)
            case(1)
                nnode = 2; nsur = nsur+1;
            case(2)
                nnode = 3; nvol = nvol+1; ntri = ntri + 1
            case(3)
                nnode = 4; nvol = nvol+1; nqua = nqua + 1
            case default
                write(*,*) 'ERROR: Unsupported Gmsh element type = ', etype; stop
            end select
            allocate(tags(ntags)); allocate(node_list(nnode))
            read(ctemp1,*) eid, etype, ntags, (tags(j),j=1,ntags), (node_list(j),j=1,nnode)
            ele_type(i) = etype; ele_phys(i) = tags(1); ele_nnode(i) = nnode
            ele_conn(1:nnode,i) = node_list
            deallocate(tags, node_list)
        enddo
        
        !calculate number of surface
        nsur_b   = nsur
        nsur_tot = (nsur_b+3*ntri+4*nqua)/2
        nsur_int = nsur_tot - nsur_b
    case(3)
        do i=1,nele
            read(2,'(a)') ctemp1; read(ctemp1,*) eid, etype, ntags
            select case(etype)
            case(2)
                nnode = 3;  nsur = nsur+1; ntri = ntri+1
            case(3)
                nnode = 4;  nsur = nsur+1; nqua = nqua+1
            case(4)
                nnode = 4;  nvol = nvol+1; ntet = ntet+1
            case(5)
                nnode = 8;  nvol = nvol+1; nhex = nhex+1
            case(6)
                nnode = 6;  nvol = nvol+1; npri = npri+1
            case default
                write(*,*) 'ERROR: Unsupported Gmsh element type = ', etype; stop
            end select

            allocate(tags(ntags)); allocate(node_list(nnode))
            read(ctemp1,*) eid, etype, ntags, (tags(j),j=1,ntags), (node_list(j),j=1,nnode)
            ele_type(i) = etype; ele_phys(i) = tags(1); ele_nnode(i) = nnode
            ele_conn(1:nnode,i) = node_list
            deallocate(tags, node_list)
        enddo
        !calculate number of surface
        nsur_b   = nsur
        nsur_tot = (nsur_b+4*ntet+6*nhex+5*npri)/2
        nsur_int = nsur_tot - nsur_b
    case default
        write(*,*) 'ERROR: Unsupported dimension = ', dim; stop
    end select



    allocate(sur_con(4,nsur_tot),sur_nv(nsur_tot),sur_id(nsur_tot))
    allocate(vol_con(8,nvol),vol_nv(nvol),vol_id(nvol))
    sur_con = 0; sur_nv = 0; sur_id = 0;
    vol_con = 0; vol_nv = 0; vol_id = 0;

    nsur = nsur_int
    nvol = 0

    select case(dim)
    case(1)
        do i=1,nele
            select case(ele_type(i))
            case(15)
                nsur = nsur+1
                sur_con(1:ele_nnode(i),nsur) = ele_conn(1:ele_nnode(i),i)
                sur_nv(nsur) = ele_nnode(i)
                sur_id(nsur) = ele_phys(i)
            case(1)
                nvol = nvol+1
                vol_con(1:ele_nnode(i),nvol) = ele_conn(1:ele_nnode(i),i)
                vol_nv(nvol) = ele_nnode(i)
                vol_id(nvol) = ele_phys(i)
            end select
        enddo
    case(2)
        do i=1,nele
            select case(ele_type(i))
            case(1)
                nsur = nsur+1
                sur_con(1:ele_nnode(i),nsur) = ele_conn(1:ele_nnode(i),i)
                sur_nv(nsur) = ele_nnode(i)
                sur_id(nsur) = ele_phys(i)
            case(2,3)
                nvol = nvol+1
                vol_con(1:ele_nnode(i),nvol) = ele_conn(1:ele_nnode(i),i)
                vol_nv(nvol) = ele_nnode(i)
                vol_id(nvol) = ele_phys(i)
            end select
        enddo
    case(3)
        do i=1,nele
            select case(ele_type(i))
            case(2,3)
                nsur = nsur+1
                sur_con(1:ele_nnode(i),nsur) = ele_conn(1:ele_nnode(i),i)
                sur_nv(nsur) = ele_nnode(i)
                sur_id(nsur) = ele_phys(i)
            case(4,5,6)
                nvol = nvol+1
                vol_con(1:ele_nnode(i),nvol) = ele_conn(1:ele_nnode(i),i)
                vol_nv(nvol) = ele_nnode(i)
                vol_id(nvol) = ele_phys(i)
            end select
        enddo
    case default
        write(*,*) 'ERROR: Unsupported dimension = ', dim; stop
    end select

    x_ver = x_ver*scale
end subroutine read_msh
