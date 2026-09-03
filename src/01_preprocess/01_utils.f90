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
