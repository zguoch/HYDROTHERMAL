SUBROUTINE closef(flag)
  ! ... Closes and deletes files
  USE f_units
  USE bc, ONLY: unconfined
  USE control
  USE i_c
  USE units
  USE variables
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: flag
  !
  CHARACTER (LEN=6) :: st(60)
  INTEGER :: ifu
  INTEGER :: fu(20)=(/30,10,11,12,13,14,15,16,17,18,19,33,20,21,22,23,  &
       24,26,32,34/)
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.25 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  !...
  IF(flag == 1) THEN
     WRITE(fustdout,2002) '*** Job Aborted Due To Input Errors ***',  &
          'Please examine Calc_log file and Out_Probdefine file'
     WRITE(fuclog,2002) '*** Job Aborted Due to Input Errors ***',  &
          'Please examine Out_Probdefine file'
     CALL errprint(1,200)
  ELSE IF(flag == 2) THEN
     WRITE(fustdout,2002) '*** Job Aborted Due To Execution Errors ***',  &
          'Please examine Calc_log file and Out_Probdefine file'
     WRITE(fuclog,2002) '*** Job Aborted Due To Execution Errors ***',  &
          'Please examine Out_Probdefine file'
     IF (iyr == 2) THEN
        WRITE (fuclog,9075) 'Step Number',ltime,'   Simulation time = ',  &
             year(time),' (yr)'
9075    FORMAT (tr2,a,i5,a,1pg12.6,a)
        WRITE (fustdout,9075) 'Step Number',ltime,'   Simulation time = ',  &
             year(time),' (yr)'
     ELSE
        WRITE (fuclog,9075) 'Step Number',ltime,'   Simulation time = ',  &
             time,' (s)'
        WRITE (fustdout,9075) 'Step Number',ltime,'   Simulation time = ',  &
             time,' (s)'
     END IF
2002 FORMAT(//tr10,a/tr15,a/)
     CALL errprint(1,200)
  ELSE IF(flag == 99) THEN
     DO  ifu=1,20
        WRITE(fu(ifu),3003) '***** Simulation Completed ***** '
3003    FORMAT(//tr40,a)
        IF (iyr == 2) THEN
           WRITE(fu(ifu),3004) 'Last time value calculated '//dots,  &
                year(time),'(yr)', 'Last time step index '//dots,ltime
3004       FORMAT(tr15,a70,1PG11.4,tr2,a/tr15,a70,i8)
        ELSE
           WRITE(fu(ifu),3004) 'Last time value calculated '//dots,  &
                time,'(s)', 'Last time step index '//dots,ltime
        END IF
     END DO
  END IF
  ! ... Close the input file
  CLOSE(fuinc,STATUS='KEEP')
  ! ... Close and delete the stripped input file
  CLOSE(fuins,STATUS='DELETE')
  ! ... If under GUI, close the file containing stdout messages
  IF(gui) CLOSE(fustdout,STATUS='KEEP')
  CLOSE(fuclog,STATUS='KEEP')
  CLOSE(fupdef,STATUS='KEEP')
  ! ... Set disposition status for optional output files
  st = 'DELETE'
  IF(ABS(print_press) > 0.) st(fup) = 'KEEP'
  IF(unconfined) st(fuwtelev) = 'KEEP'
  IF(ABS(print_enth) > 0.) st(fuen) = 'KEEP'
  IF(ABS(print_temp) > 0.) st(fut) = 'KEEP'
  IF(ABS(print_satn) > 0.) st(fusat) = 'KEEP'
  IF(ABS(print_resid) > 0.) st(fures) = 'KEEP'
  IF(ABS(print_dimno) > 0.) st(fudimno) = 'KEEP'
  IF(ABS(print_dens) > 0.) st(fuden) = 'KEEP'
  IF(ABS(print_vis) > 0.) st(fuvis) = 'KEEP'
  IF(ABS(print_vel) > 0.) st(fuvel) = 'KEEP'
  IF(ABS(print_pot) > 0.) st(fupot) = 'KEEP'
  IF(ABS(print_poros) > 0.) st(fupor) = 'KEEP'
  IF(ABS(print_perm) > 0.) st(fuperm) = 'KEEP'
  IF(ABS(print_pmprop) > 0.) st(futhp) = 'KEEP'
!!$  IF(ABS(print_source) > 0. .AND. nsrce > 0) THEN
  IF(ABS(print_source) > 0.) st(fusrc) = 'KEEP'
  IF(ABS(print_bcflow) > 0.) THEN
     st(fubcf) = 'KEEP'
     st(fubcf2) = 'KEEP'
  END IF
  IF(ABS(print_balance) > 0.) st(fubal) = 'KEEP'
  IF(ABS(print_plotscalar) > 0.) THEN
     st(fupltsca) = 'KEEP'
     st(fupltsca2) = 'KEEP'
  END IF
  IF(ABS(print_plotvector) > 0.) st(fupltvec) = 'KEEP'
  IF(ABS(print_temporal) > 0.) st(fuplttimser) = 'KEEP'
  IF(ABS(print_dump) > 0.) st(fuorst) = 'KEEP'
  ! ... close and dispose
  CLOSE(fup,STATUS=st(fup))
  CLOSE(fuen,STATUS=st(fuen))
  CLOSE(fut,STATUS=st(fut))
  CLOSE(fusat,STATUS=st(fusat))
  CLOSE(fuwtelev,STATUS=st(fuwtelev))
  CLOSE(fuden,STATUS=st(fuden))
  CLOSE(fuvis,STATUS=st(fuvis))
  CLOSE(fupot,STATUS=st(fupot))
  CLOSE(fuvel,STATUS=st(fuvel))
  CLOSE(fubcf,STATUS=st(fubcf))
  CLOSE(fubcf2,STATUS=st(fubcf2))
  CLOSE(fusrc,STATUS=st(fusrc))
  CLOSE(futhp,STATUS=st(futhp))
  CLOSE(fupor,STATUS=st(fupor))
  CLOSE(fuperm,STATUS=st(fuperm))
  CLOSE(fubal,STATUS=st(fubal))
  CLOSE(fuorst,STATUS=st(fuorst))
  CLOSE(fures,STATUS=st(fures))
  CLOSE(fudimno,STATUS=st(fudimno))
  CLOSE(fupltsca,STATUS=st(fupltsca))
  CLOSE(fupltsca2,STATUS=st(fupltsca2))
  CLOSE(fupltvec,STATUS=st(fupltvec))
  CLOSE(fuplttimser,STATUS=st(fuplttimser))
  IF(initphi) CLOSE(fuicpp,STATUS='KEEP')
END SUBROUTINE closef
