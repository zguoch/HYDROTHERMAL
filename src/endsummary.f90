SUBROUTINE endsummary
  ! ... Purpose:  To print summary data for simulation
  USE f_units
  USE bc
  USE control
  USE fdeq
  USE mesh
  USE units
  USE variables
  IMPLICIT NONE
  INTEGER :: i2phase, ic, iawphase, icphase, isphase, iwphase, ll, lll
  !..      REAL elpsecs
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.8 $//$Date: 2007/05/08 23:45:00 $'
  !     ------------------------------------------------------------------
  WRITE(fuclog, '(/tr10,A)') '--- Summary data for completed simulation ---'
  WRITE(fuclog, '(/tr10,A,I7)')  &
       'Total number of time steps:                         ', ltime
  WRITE(fuclog, '(10X,A,I7)')  &
       'Total number of Newton-Raphson iterations:          ', lnritotl
  IF (iyr == 2) THEN
     WRITE(fuclog,9005) 'Initial Time:        ',year(timestr),' (yr)'
     WRITE(fuclog,9005) 'Ending Time:         ',year(time),' (yr)'
     IF (ltime > 0) WRITE(fuclog,9005) 'Avg Time Step Length:',  &
          year((time-timestr))/REAL(ltime),' (yr)'
     WRITE(fuclog,9005) 'Maximum Time Step:   ',year(deltmax),' (yr)'
     WRITE(fuclog,9005) 'Minimum Time Step:   ',year(deltmin),' (yr)'
9005 FORMAT(tr10,a,tr25,1pg14.7,a)
  ELSE
     WRITE(fuclog,9005) 'Initial Time:        ',timestr,' (s)'
     WRITE(fuclog,9005) 'Ending Time:         ',time,' (s)'
     IF (ltime > 0) WRITE(fuclog,9005) 'Avg Time Step Length:',  &
          (time-timestr)/REAL(ltime),' (s)'
     WRITE(fuclog,9005) 'Maximum Time Step:   ',deltmax,' (s)'
     WRITE(fuclog,9005) 'Minimum Time Step:   ',deltmin,' (s)'
  END IF
  WRITE(fuclog,'(tr10,A,I3)')  &
       'Maximum number of Newton-Raphson iterations for a time step: ', ilnrimax
  IF(ltime > 0) WRITE(fuclog,'(tr10,A,f5.1)')  &
       'Average number of Newton-Raphson iterations per time step: ',  &
       REAL(lnritotl)/REAL(ltime)
  IF(slmeth == 2) THEN
     WRITE(fuclog,'(tr10,A,I3)')  &
          'Maximum number of GMRES iterations for a solve:              ', iter_gm_max
     IF(ltime > 0) WRITE(fuclog,'(tr10,A,f5.1)')  &
          'Average number of GMRES iterations per solve:              ',  &
          REAL(tot_iter_gm)/REAL(lnritotl)
  END IF
  WRITE(fuclog,9010) itmcntl(2),itmcntl(3),itmcntl(1),itmcntl(4),itmcntl(7),  &
       itmcntl(8),itmcntl(9),itmcntl(5),itmcntl(6),itmcntl(10)
9010 FORMAT (/tr10,'Factors Controlling Time Step Length:'/ i17,  &
       ' - (P) Maximum pressure change'/ i17,  &
       ' - (H) Maximum enthalpy change'/ i17,  &
       ' - (W) Maximum water saturation change'/ i17,  &
       ' - (M) Maximum time step increase'/ i17,  &
       ' - (D) Time to read new input data'/ i17,  &
       ' - (O) Time to write output data'/ i17,  &
       ' - (R) Reset to value before reading new data'/i17,  &
       ' - (R) Reset to value before writing output data'/ i17,  &
       ' - (I) Value for initial time step'/ i17,  &
       ' - (M) Maximum time step allowed')
  WRITE(fuclog,9015) itmcntl(12),itmcntl(13),itmcntl(11),itmcntl(14),itmcntl(15),  &
       itmcntl(16),itmcntl(17),itmcntl(18),itmcntl(20),itmcntl(19),itmcntl(21)
9015 FORMAT (/tr10,'Reasons for Cutting the Time Step Length:'/ i17,  &
       ' - (P) Pressure change exceeded maximum allowed'/ i17,  &
       ' - (H) Enthalpy change exceeded maximum allowed'/ i17,  &
       ' - (W) Saturation change exceeded maximum allowed'/ i17,  &
       ' - (C) Excessive multiple nodes changed phase in a N-R iteration'/i17,  &
       ' - (C) Excessive multiple nodes changed phase in a time step' / i17,  &
       ' - (I) Maximum N-R iterations reached without convergence',  &
       /, i17, ' - (N) Newton-Raphson iterations stalled out',  &
       /, i17, ' - (T) Pressure or enthalpy exceeded limits of enthalpy tables',  &
       /, i17, ' - (T) Pressure below limits of enthalpy tables',  &
       /, i17, ' - (L) Pressure exceeded lithostatic limit',  &
       /, i17, ' - (G) GMRES solver failed to converge')
  WRITE(fuclog, 9020)
9020 FORMAT(/tr10,'Convergence Criteria Satisfied:')
  WRITE(fuclog,9022) (iconvrg(ll),ll=1,3)
9022 FORMAT(i17,' - (R) Mass and energy residuals (input)'/  &
       i17,' - (C) Absolute changes in pressure and enthalpy (hard coded)'/  &
       i17,' - (%) Percent changes in pressure and enthalpy from previous step (input)')
  WRITE(fuclog,9021) 'Mass:  ', resmaxmall, ltimresmaxm
  WRITE(fuclog,9021) 'Energy:', resmaxeall, ltimresmaxe
9021 FORMAT(tr18,'Maximum residual for ', a,tr1,1PE9.2, '; at Time Step', i5)
  ! ... Determine the number of nodes in region of each phase type
  iwphase = 0
  i2phase = 0
  isphase = 0
  icphase = 0
  iawphase = 0
  DO  lll = 1, npiccons
     ic = npic(lll)
     IF (ind(ic) == 1) iwphase = iwphase + 1
     IF (ind(ic) == 2) i2phase = i2phase + 1
     IF (ind(ic) == 3) isphase = isphase + 1
     IF (ind(ic) == 4) icphase = icphase + 1
     IF (ind(ic) == 5) iawphase = iawphase + 1
  END DO
  WRITE(fuclog, '(/tr10,A)') 'Distribution of Phases at End of Simulation: '
  WRITE(fuclog, 9030) 'Number of Compressed Water Nodes ....', iwphase
  WRITE(fuclog, 9030) 'Number of 2-phase Water/Steam Nodes .', i2phase
  WRITE(fuclog, 9030) 'Number of Superheated Steam Nodes ...', isphase
  WRITE(fuclog, 9030) 'Number of Supercritical Fluid Nodes .', icphase
  WRITE(fuclog, 9030) 'Number of Air-Water Nodes ...........', iawphase
9030 FORMAT (tr12,a,i5)
  WRITE(fuclog,'(/tr10,A,I8/tr15,a)')  &
       'Number of nodes changing phase during simulation:', iphschgtot,  &
       '(to/from 2-phase or air-water phases only)' 
  IF(unconfined) WRITE(fuclog,'(tr10,A,I8)')  &
       'Number of values of pressure less than 0.5 atm:  ', itmcntl(20)
!****This is done better with a Fortran 90 cpu timing routine
!!$  ! ...  CALL TIMER(elpsecs)
!!$  ! ...  WRITE(fuclog, 9025) elpsecs, elpsecs/60.0
!!$  ! ...  WRITE(*, 9025) elpsecs, elpsecs/60.0
!!$  !.. 9025 FORMAT (/, 10X, 'Total system time used =', 1P, G12.5, ' (s);',
!!$  !..     &        G14.4, ' (min)')
END SUBROUTINE endsummary
