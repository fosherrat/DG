    module m_parameter
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    double precision, parameter :: pi = 3.1415926535897932384626d0
    double precision, parameter :: gamma = 1.4d0

    end module m_parameter





    module m_folder
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    character(len=6) :: folder_input = 'input/'
    character(len=7) :: folder_output = 'output/'
    character(len=6) :: folder_point = 'point/'
    character(len=32) :: file_point = 'point_flow.dat'
        
    end module m_folder
    
    
    
    
    
    module m_input
    
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! input file variables
    integer :: order, flux, time, diffusion
    double precision :: cfl, t_final
    
    ! calculated variables
    integer :: ndof, nvar

    end module m_input





    module m_point_dg
    
    implicit none

    integer :: dim
    integer :: nver, nsur_tot, nsur_int, nsur_b, nvol
    integer :: ntri, nqua, ntet, npri, nhex
    integer, allocatable :: sur_con(:,:), sur_nv(:), sur_id(:)
    integer, allocatable :: vol_con(:,:), vol_nv(:), vol_id(:)
    integer, allocatable :: con_sur_vol(:,:)
    double precision, allocatable :: x_ver(:,:)
    double precision, allocatable :: vol(:), vol_cen(:,:)
    double precision, allocatable :: sur(:), sur_cen(:,:), sur_vec(:,:)
    double precision, allocatable :: jcb(:,:,:), det_jcb(:), inv_jcb(:,:,:)

    end module m_point_dg





    module m_flow

    implicit none

    double precision, allocatable :: conv(:,:,:), diff(:,:,:)
    double precision, allocatable :: cons(:,:,:), prim(:,:,:), grad_cons(:,:,:)
    
    ! under is only for 1D. need to be generalized
    integer :: nquad
    double precision, allocatable :: legendre(:,:), grad_legendre(:,:)
    double precision, allocatable :: mass(:,:)
    double precision, allocatable :: gauss_x(:), gauss_w(:)

    end module m_flow
