program IFS_01_preprocessing

use m_folder
use m_input

implicit none

double precision :: tic, toc
character*8 :: timenow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,'(4/)')
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Start 01_preprocessing!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                                
call cpu_time(tic);    call date_and_time(timenow); write(*,*) 'Program start:                 ', timenow; write(*,*)
call read_input;       call date_and_time(timenow); write(*,*) '[Setting] Read_input complete: ', timenow
call read_msh;         call date_and_time(timenow); write(*,*) '[Flow] read_mesh complete:     ', timenow
call mesh_setting_dg;  call date_and_time(timenow); write(*,*) '[Flow] mesh_setting complete:  ', timenow
call post_dg;          call date_and_time(timenow); write(*,*) '[Flow] post complete           ', timenow
call deallocate_m;     call date_and_time(timenow); write(*,*) '[Flow] deallocate complete:    ', timenow
write(*,*); call cpu_time(toc); call date_and_time(timenow); write(*,*) 'Program end:                   ', timenow, toc-tic

write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End 01_preprocessing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

end program IFS_01_preprocessing
