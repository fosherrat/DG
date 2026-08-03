    subroutine deallocate_m

    use m_folder
    use m_input
    use m_point_dg

    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(allocated(sur_con)) deallocate(sur_con)
    if(allocated(sur_id)) deallocate(sur_id)
    if(allocated(vol_con)) deallocate(vol_con)
    if(allocated(vol_id)) deallocate(vol_id)
    if(allocated(x_ver)) deallocate(x_ver)

    if(allocated(con_sur_vol)) deallocate(con_sur_vol)
    if(allocated(vol_nv)) deallocate(vol_nv)
    if(allocated(sur_nv)) deallocate(sur_nv)
    if(allocated(vol)) deallocate(vol)
    if(allocated(vol_cen)) deallocate(vol_cen)
    if(allocated(sur)) deallocate(sur)
    if(allocated(sur_cen)) deallocate(sur_cen)
    if(allocated(sur_vec)) deallocate(sur_vec)

    end subroutine deallocate_m
