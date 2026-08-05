subroutine deallocate_m()
    use m_point_dg
    use m_flow

    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(allocated(sur_con)) deallocate(sur_con)
    if(allocated(sur_nv)) deallocate(sur_nv)
    if(allocated(sur_id)) deallocate(sur_id)
    if(allocated(vol_con)) deallocate(vol_con)
    if(allocated(vol_nv)) deallocate(vol_nv)
    if(allocated(vol_id)) deallocate(vol_id)
    if(allocated(con_sur_vol)) deallocate(con_sur_vol)

    if(allocated(x_ver)) deallocate(x_ver)
    if(allocated(vol)) deallocate(vol)
    if(allocated(vol_cen)) deallocate(vol_cen)
    if(allocated(sur)) deallocate(sur)
    if(allocated(sur_cen)) deallocate(sur_cen)
    if(allocated(sur_vec)) deallocate(sur_vec)
    if(allocated(jcb)) deallocate(jcb)
    if(allocated(det_jcb)) deallocate(det_jcb)
    if(allocated(inv_jcb)) deallocate(inv_jcb)

    if(allocated(conv)) deallocate(conv)
    if(allocated(diff)) deallocate(diff)
    if(allocated(cons)) deallocate(cons)
    if(allocated(prim)) deallocate(prim)
    if(allocated(grad_cons)) deallocate(grad_cons)
    if(allocated(shock_sensor)) deallocate(shock_sensor)
    if(allocated(av_elem)) deallocate(av_elem)
    if(allocated(av_face)) deallocate(av_face)
    if(allocated(av_end)) deallocate(av_end)
    if(allocated(legendre)) deallocate(legendre)
    if(allocated(grad_legendre)) deallocate(grad_legendre)
    if(allocated(mass)) deallocate(mass)
    if(allocated(gauss_x)) deallocate(gauss_x)
    if(allocated(gauss_w)) deallocate(gauss_w)

end subroutine deallocate_m
