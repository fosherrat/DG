subroutine poly_1d(coeff, x, val)
    use m_parameter
    use m_input

    implicit none

    double precision, intent(in) :: coeff(:), x
    double precision, intent(out) :: val
    integer :: i, n
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    val = 0.d0; n = size(coeff)
    do i=1,n
        val = val + coeff(i)*x**(i-1)
    enddo

end subroutine poly_1d