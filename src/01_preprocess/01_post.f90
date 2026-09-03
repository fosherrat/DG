subroutine post_dg()
    use m_parameter
    use m_folder
    use m_input
    use m_point_dg

    implicit none

    integer :: i
    !==========================================================================================

    open(4,file=''//folder_point//'point_flow.dat',status='replace')
    write(4,*) dim
    write(4,*) nver
    write(4,*) nsur_tot, nsur_int, nsur_b
    write(4,*) nvol, ntri, nqua, ntet, npri, nhex

    do i=1,nver
        write(4,*) i, x_ver(:,i)
    enddo

    do i=1,nvol
        write(4,*) vol_id(i), vol_nv(i), vol(i), vol_cen(:,i), vol_con(:,i)
    enddo

    do i=1,nsur_tot
        write(4,*) sur_id(i), sur_nv(i), sur(i), sur_cen(:,i), sur_vec(:,i), con_sur_vol(:,i)
    enddo
    close(4)
end subroutine post_dg