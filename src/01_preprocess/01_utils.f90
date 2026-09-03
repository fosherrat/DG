subroutine match(ivol, ilocal, face_con, size, s)
    use m_point_dg

    implicit none

    integer, intent(in) :: ivol, ilocal, size
    integer, intent(in) :: face_con(4)
    integer, intent(inout) :: s

    integer :: j
    integer :: key(4)
    logical :: found
    !==========================================================================================

    key = face_con
    call sort(key,size)

    found = .false.

    do j=1,nsur_tot
        if(all(sur_key(:,j).eq.key)) then
            found = .true.

            if(sur_nv(j).ne.size) then
                write(*,*) 'ERROR: inconsistent face vertex count'
                stop
            endif

            if(j.gt.nsur_int) then
                ! Boundary face: its only incident element is the left element.
                if(con_sur_vol(1,j).ne.-1) then
                    write(*,*) 'ERROR: duplicated boundary face = ', j
                    stop
                endif

                con_sur_vol(1,j) = ivol
                con_sur_local(1,j) = ilocal
                con_vol_sur(ilocal,ivol) = j

                ! Use the left element's cyclic ordering for every global face.
                sur_con(:,j) = face_con

            else
                ! Existing interior face: this element is the right element.
                if(con_sur_vol(2,j).ne.-1) then
                    write(*,*) 'ERROR: non-manifold face = ', j
                    stop
                endif

                con_sur_vol(2,j) = ivol
                con_sur_local(2,j) = ilocal
                con_vol_sur(ilocal,ivol) = j
            endif

            exit
        endif
    enddo

    if(.not.found) then
        if(s.gt.nsur_int) then
            write(*,*) 'ERROR: number of interior faces exceeds estimate'
            stop
        endif

        ! New interior face: first encountered element is the left element.
        sur_con(:,s) = face_con
        sur_key(:,s) = key
        sur_nv(s) = size

        con_sur_vol(1,s) = ivol
        con_sur_local(1,s) = ilocal
        con_vol_sur(ilocal,ivol) = s

        s = s + 1
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
        j = i-1
        do while (j.ge.1)
            if(temp1(j) <= value) exit
            temp1(j+1) = temp1(j); j = j-1
        end do
        temp1(j+1) = value
    end do
end subroutine sort


subroutine cross(a,b,c)
    implicit none

    double precision, intent(in) :: a(3), b(3)
    double precision, intent(out) :: c(3)
    !==========================================================================================
    ! c = a x b
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end subroutine cross
