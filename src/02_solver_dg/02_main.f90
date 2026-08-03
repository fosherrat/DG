program IFS_02_solver_dg

use m_input

implicit none

integer :: iter
double precision :: tic, toc, dt, t_current
character(len=8) :: timenow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,'(4/)')
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Start 02_solver_dg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

call cpu_time(tic);         call date_and_time(timenow); write(*,*) 'Program start:                 ', timenow; write(*,*)
call read_point;            call date_and_time(timenow); write(*,*) '[Setting] read_point complete: ', timenow
call read_input;            call date_and_time(timenow); write(*,*) '[Setting] read_input complete: ', timenow
call build_basis;           call date_and_time(timenow); write(*,*) '[Flow] build_basis complete:   ', timenow
call build_mass;            call date_and_time(timenow); write(*,*) '[Flow] build_mass complete:    ', timenow
call build_gauss_quadrature call date_and_time(timenow); write(*,*) '[Flow] build_gauss complete:   ', timenow
call init_sin;              call date_and_time(timenow); write(*,*) '[Flow] init_sin complete:      ', timenow

select case(time)
case(0)
    t_current = 0.d0
    dt = 0.d0
    iter = 0

    do while (t_current.lt.t_final)
        call euler_explicit(dt, t_final-t_current)
        if (dt.le.0.d0) then
            write(*,*) 'ERROR: non-positive dt in time loop. dt = ', dt
            stop
        endif
        t_current = t_current + dt
        iter = iter + 1
    enddo

    call date_and_time(timenow)
    write(*,*) '[Time] euler_explicit complete:', timenow, ' t = ', t_current, ' dt = ', dt, ' iter = ', iter
case default
    write(*,*) 'ERROR: Unsupported time = ', time
    stop
end select

call deallocate_m;  call date_and_time(timenow); write(*,*) '[Flow] deallocate complete:    ', timenow; write(*,*) 
call cpu_time(toc); call date_and_time(timenow); write(*,*) 'Program end:                   ', timenow, toc-tic

write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End 02_solver_dg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

end program IFS_02_solver_dg
