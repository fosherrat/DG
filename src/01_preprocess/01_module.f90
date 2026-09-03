    module m_parameter
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    double precision, parameter :: pi = 3.1415926535897932384626d0

    end module m_parameter
    
    
    
    
    
    module m_folder
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    character(len=6) :: folder_input = 'input/'
    character(len=5) :: folder_mesh  = 'mesh/'
    character(len=6) :: folder_point = 'point/'
    character(len=32) :: file_mesh
        
    end module m_folder
    
    
    
    
    
    module m_input
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    integer :: dim, para, mesh_type, discretization
    double precision :: scale
    
    end module m_input
    
    
    
    
    
    module m_point_dg
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! con: connected vertex id, nv: number of vertex composing an element, id: PID
    integer :: nver, nsur, nvol
    integer :: ntri, nqua, ntet, npri, nhex
    ! sur_con stores the cyclic ordering of the left volume face.
    ! sur_key stores sorted node IDs and is used only for face matching.
    integer, allocatable :: sur_con(:,:), sur_key(:,:), sur_nv(:), sur_id(:)
    integer, allocatable :: vol_con(:,:), vol_nv(:), vol_id(:)
    double precision, allocatable :: x_ver(:,:)

    ! con_sur_vol: connected left/right volume of a surface, cen: center point
    ! '': area or volume, vec: normal vector(left --> right)
    integer :: nsur_tot, nsur_int, nsur_b
    integer, allocatable :: con_sur_vol(:,:), con_vol_sur(:,:)
    double precision, allocatable :: vol(:), vol_cen(:,:)
    double precision, allocatable :: sur(:), sur_cen(:,:), sur_vec(:,:)

    ! dg specialized variables
    ! con_sur_local: local surface number of volume
    integer, allocatable :: con_sur_local(:,:), sur_perm(:,:)

    end module m_point_dg

    
    
    
    
    module m_point_fvm
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! con: connected vertex id, nv: number of vertex composing an element, id: PID
    integer :: nver, nsur, nvol
    integer :: ntri, nqua, ntet, npri, nhex
    integer, allocatable :: sur_con(:,:), sur_nv(:), sur_id(:)
    integer, allocatable :: vol_con(:,:), vol_nv(:), vol_id(:)
    double precision, allocatable :: x_ver(:,:)

    ! con_sur_vol: connected left/right volume of a surface, cen: center point
    ! '': area or volume, vec: normal vector(left --> right)
    integer :: nsur_tot, nsur_int, nsur_b
    integer, allocatable :: con_sur_vol(:,:), con_vol_sur(:,:)
    double precision, allocatable :: vol(:), vol_cen(:,:)
    double precision, allocatable :: sur(:), sur_cen(:,:), sur_vec(:,:)

    end module m_point_fvm