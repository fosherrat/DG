double precision function poly_value(coef, ncoef, x)
    implicit none

    integer, intent(in) :: ncoef
    double precision, intent(in) :: coef(ncoef), x
    integer :: i

    ! Horner evaluation of coef(1) + coef(2)*x + ... .
    poly_value = 0.d0
    do i=ncoef,1,-1
        poly_value = poly_value*x + coef(i)
    enddo

end function poly_value