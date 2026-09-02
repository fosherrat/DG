subroutine post_dg()
    use m_folder
    use m_input
    use m_point_dg
    use m_flow

    implicit none

    integer :: elem, dof

    ! One row stores one conservative DG coefficient.  The physical coordinate
    ! is x = x_center + jacobian*xi and basis_id refers to the Legendre basis.
    open(4,file=trim(folder_output)//'solution_coefficients.dat',status='replace', &
        action='write')

    write(4,'(A)') '# DG conservative solution coefficients'
    write(4,'(A,I0)') '# nvar = ', nvar
    write(4,'(A,I0)') '# ndof = ', ndof
    write(4,'(A,I0)') '# nvol = ', nvol
    write(4,'(A)') '# x = x_center + jacobian*xi; basis_id uses Legendre P_(basis_id-1)'
    write(4,'(A)') '# columns: element_id basis_id x_center jacobian cons_1 ... cons_nvar'

    do elem=1,nvol
        do dof=1,ndof
            write(4,*) elem, dof, vol_cen(1,elem), jcb(1,1,elem), cons(:,dof,elem)
        enddo
    enddo

    close(4)
end subroutine post_dg
