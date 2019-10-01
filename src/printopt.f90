SUBROUTINE printopt(ichgpr)
  ! ... Purpose:  To read and write input data selecting print/plot options
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE mesh
  USE variables, ONLY: time
  USE units
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ichgpr       ! ... number of printer options to change
                                      ! ...    99 - read in complete set of printer options
                                      ! ...    -n - [1-7] change the first n print
                                      ! ...    parameters. For format versions 1 and 2 only.
  ! ...   Ioptprta() - print option for arrays =0 - skip; =1 - print
  ! ...       (1)-pressure, (2)-enthalpy, (3)-temperature,
  ! ...            (4)-water saturation, (5)-residual mass, (6)-residual engy
  ! ...       (7)-water density,  (8)-water viscosity,  (9)-water X velocity,
  ! ...           (10)-water Y velocity, (11)-water Z velocity,
  ! ...           (12)-water hydraulic potential
  ! ...       (13)-steam density, (14)-steam viscosity, (15)-steam X velocity,
  ! ...           (16)-steam Y velocity, (17)-steam Z velocity,
  ! ...           (18)-steam hydraulic potential
  ! ...       (19)-porosity, (20)-X-permeability, (21)-X-permeability,
  ! ...            (22)-X-permeability, (23)-thermal conductivity,
  ! ...            (24)-specific heat, (25)-rock density (26)-compressibility
  ! ...       (27)-Peclet number, (28)-Nusselt number
  !
  CHARACTER (LEN=80) :: ufile
  INTEGER :: idum, iprtasum, k1, k2, kdum, ichopt, plotfile_type
  LOGICAL :: lprwater, lprsteam
  ! ...              ----- print/ploting controls -----
  ! ...    Prfreq - time interval between printing output
  ! ...    Ioptpr(1) 0 - no plotfile, 1 - Explorer plotfile,
  ! ...              4 - GNU plotfile, 5 - IDL plotfile,
  ! ...              8 - ASCII B2 plotfile; 9 - binary B2 plotfile
  ! ...    Ioptpr(2) 0 - 80 column printout; 1 - 132 column printout; 2 - 159 col
  ! ...    Ioptpr(3) 0 - short; 1 - med; 2 - long output to calc_log and screen
  ! ...    Ioptpr(4) 0 - average values not printed; 1 - average values printed
  ! ...    Ioptpr(5) 0 - don't print source/sink and specified flux; 1 - do print
  ! ...    Ioptpr(6)  - number of nodes for which data are printed to Plot_timeseries
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.7 $//$Date: 2007/12/20 20:32:47 $'
  !     ------------------------------------------------------------------
  !...
  IF(informt <= 2) THEN
    CALL read_print_control_2
  ELSE IF(informt == 3) THEN
    CALL read_print_control_3
  END IF
  ! write the print controls selected by user
     IF (print_press < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Pressure:','every ',  &
             INT(ABS(print_press)),'th time step'
9105    FORMAT(tr10,a,t53,a,i5,a)
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Pressure: ',year(print_press),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Pressure: ',print_press,' (s)'
        END IF
9205    FORMAT(tr10,a,t53,1PG14.5,a)
     END IF
     IF (print_enth < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Enthalpy:','every ',  &
             INT(ABS(print_enth)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Enthalpy: ',year(print_enth),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Enthalpy: ',print_enth,' (s)'
        END IF
     END IF
     IF (print_temp < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Temperature:','every ',  &
             INT(ABS(print_temp)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Temperature: ',year(print_temp),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Temperature: ',print_temp,' (s)'
        END IF
     END IF
     IF (print_satn < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Saturation:','every ',  &
             INT(ABS(print_satn)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Saturation: ',year(print_satn),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Saturation: ',print_satn,' (s)'
        END IF
     END IF
     IF (print_dens < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Density:','every ',  &
             INT(ABS(print_dens)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Density: ',year(print_dens),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Density: ',print_dens,' (s)'
        END IF
     END IF
     IF (print_vis < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Viscosity:','every ',  &
             INT(ABS(print_vis)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Viscosity: ',year(print_vis),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Viscosity: ',print_vis,' (s)'
        END IF
     END IF
     IF (print_pot < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Potential:','every ',  &
             INT(ABS(print_pot)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Potential: ',year(print_pot),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Potential: ',print_pot,' (s)'
        END IF
     END IF
     IF (print_vel < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Velocity:','every ',  &
             INT(ABS(print_vel)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Velocity: ',year(print_vel),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Velocity: ',print_vel,' (s)'
        END IF
     END IF
     IF (print_bcflow < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, B.C. flow rates:','every ',  &
             INT(ABS(print_bcflow)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, B.C. flow rates: ',year(print_bcflow),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, B.C. flow rates: ',print_bcflow,' (s)'
        END IF
     END IF
     IF (print_source < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Source flow rates:','every ',  &
             INT(ABS(print_source)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Source flow rates: ',year(print_source),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Source flow rates: ',print_source,' (s)'
        END IF
     END IF
     IF (print_pmprop < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Porous Media Properties:','every ',  &
             INT(ABS(print_pmprop)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Porous media properties: ',year(print_pmprop),  &
                ' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Porous media properties: ',print_pmprop,' (s)'
        END IF
     END IF
     IF (print_poros < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Porosity:','every ',  &
             INT(ABS(print_poros)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Porosity: ',year(print_poros),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Porosity: ',print_poros,' (s)'
        END IF
     END IF
     IF (print_perm < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Permeability:','every ',  &
             INT(ABS(print_perm)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Permeability: ',year(print_perm),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Permeability: ',print_perm,' (s)'
        END IF
     END IF
     IF (print_balance < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Balance Tables:','every ',  &
             INT(ABS(print_balance)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Balance tables: ',year(print_balance),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Balance tables: ',print_balance,' (s)'
        END IF
     END IF
     IF (print_dump < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Restart data:','every ',  &
             INT(ABS(print_dump)),'th time step'
     ELSE
        IF (informt < 3) THEN
           WRITE(fupdef,'(tr10,a,t50,a)') 'Output Interval, Restart data: ','at end of simulation period'
        ELSE IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Restart data: ',year(print_dump),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Restart data: ',print_dump,' (s)'
        END IF
     END IF
     IF (print_resid < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Residuals:','every ',  &
             INT(ABS(print_resid)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Residuals: ',year(print_resid),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Residuals: ',print_resid,' (s)'
        END IF
     END IF
     IF (print_dimno < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Dimensionless Numbers:','every ',  &
             INT(ABS(print_dimno)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Dimensionless numbers: ',year(print_dimno),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Dimensionless numbers: ',print_dimno,' (s)'
        END IF
     END IF
     IF (print_plotscalar < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Scalar Plot Data:','every ',  &
             INT(ABS(print_plotscalar)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Scalar plot data: ',year(print_plotscalar),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Scalar plot data: ',print_plotscalar,' (s)'
        END IF
     END IF
     IF (print_plotvector < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Vector Plot Data:','every ',  &
             INT(ABS(print_plotvector)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Vector plot data: ',year(print_plotvector),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Vector plot data: ',print_plotvector,' (s)'
        END IF
     END IF
     IF (print_temporal < 0._kdp) THEN
        WRITE (fupdef,9105)  &
             'Output Interval, Time Series Plot Data:','every ',  &
             INT(ABS(print_temporal)),'th time step'
     ELSE
        IF (iyr == 2) THEN
           WRITE(fupdef,9205) 'Output Interval, Time series plot data: ',  &
                year(print_temporal),' (yr)'
        ELSE
           WRITE(fupdef,9205) 'Output Interval, Time series plot data: ',print_temporal,' (s)'
        END IF
     END IF
     IF (ABS(print_source) > 0._kdp) WRITE (fuclog, 9015)  &
          'Fluxes for source/sink nodes printed to Out_source'
     IF (ABS(print_bcflow) > 0._kdp) WRITE (fuclog, 9115)  &
          'Fluxes for specified value nodes, specified flux nodes, ',  &
          'and seepage nodes printed to Out_bcflow'
     ! ... Set print flag stop times to next time after current time
!!$     tchg = 1.e99_kdp     ! ... Set to large number since it has not yet been read
     ! ...   already set in init_ht
     IF(print_press > 0._kdp) THEN
        timprp=(1._kdp+INT(time/print_press))*print_press
     ELSE
        timprp=tchg
     ENDIF
     IF(print_enth > 0._kdp) THEN
        timprh=(1._kdp+INT(time/print_enth))*print_enth
     ELSE
        timprh=tchg
     ENDIF
     IF(print_temp > 0._kdp) THEN
        timprt=(1._kdp+INT(time/print_temp))*print_temp
     ELSE
        timprt=tchg
     ENDIF
     IF(print_satn > 0._kdp) THEN
        timprsat=(1._kdp+INT(time/print_satn))*print_satn
     ELSE
        timprsat=tchg
     ENDIF
     IF(print_dens > 0._kdp) THEN
        timprd=(1._kdp+INT(time/print_dens))*print_dens
     ELSE
        timprd=tchg
     ENDIF
     IF(print_vis > 0._kdp) THEN
        timprv=(1._kdp+INT(time/print_vis))*print_vis
     ELSE
        timprv=tchg
     ENDIF
     IF(print_pot > 0._kdp) THEN
        timprpot=(1._kdp+INT(time/print_pot))*print_pot
     ELSE
        timprpot=tchg
     ENDIF
     IF(print_vel > 0._kdp) THEN
        timprvel=(1._kdp+INT(time/print_vel))*print_vel
     ELSE
        timprvel=tchg
     ENDIF
     IF(print_bcflow > 0._kdp) THEN
        timprbcf=(1._kdp+INT(time/print_bcflow))*print_bcflow
     ELSE
        timprbcf=tchg
     ENDIF
     IF(print_source > 0._kdp) THEN
        timprsrc=(1._kdp+INT(time/print_source))*print_source
     ELSE
        timprsrc=tchg
     ENDIF
     IF(print_pmprop > 0._kdp) THEN
        timprprop=(1._kdp+INT(time/print_pmprop))*print_pmprop
     ELSE
        timprprop=tchg
     ENDIF
     IF(print_poros > 0._kdp) THEN
        timprpor=(1._kdp+INT(time/print_poros))*print_poros
     ELSE
        timprpor=tchg
     ENDIF
     IF(print_perm > 0._kdp) THEN
        timprperm=(1._kdp+INT(time/print_perm))*print_perm
     ELSE
        timprperm=tchg
     ENDIF
     IF(print_balance > 0._kdp) THEN
        timprbal=(1._kdp+INT(time/print_balance))*print_balance
     ELSE
        timprbal=tchg
     ENDIF
     IF(print_dump > 0._kdp) THEN
        timprdump=(1._kdp+INT(time/print_dump))*print_dump
     ELSE
        timprdump=tchg
     ENDIF
     IF(print_resid > 0._kdp) THEN
        timprresid=(1._kdp+INT(time/print_resid))*print_resid
     ELSE
        timprresid=tchg
     ENDIF
     IF(print_dimno > 0._kdp) THEN
        timprdimno=(1._kdp+INT(time/print_dimno))*print_dimno
     ELSE
        timprdimno=tchg
     ENDIF
     IF(print_plotscalar > 0._kdp) THEN
        timprpscal=(1._kdp+INT(time/print_plotscalar))*print_plotscalar
     ELSE
        timprpscal=tchg
     ENDIF
     IF(print_plotvector > 0._kdp) THEN
        timprpvec=(1._kdp+INT(time/print_plotvector))*print_plotvector
     ELSE
        timprpvec=tchg
     ENDIF
     IF(print_temporal > 0._kdp) THEN
        timprtemp=(1._kdp+INT(time/print_temporal))*print_temporal
     ELSE
        timprtemp=tchg
     ENDIF
     ! ... Set the time for next printout
     prtimenxt = MIN(tchg, timprp, timprh, timprt, timprsat, timprd, timprv, timprpot,  &
          timprvel, timprbcf, timprsrc, timprprop, timprpor, timprperm, timprbal, timprdump,  &
          timprresid, timprdimno, timprpscal, timprpvec, timprtemp)
  WRITE (fuclog, 9010) 'Time series data written to '//  &
       'Plot_timeseries file for', n_locations,' nodes'
9010 FORMAT (tr10,a,i3,a)
  ! ... *****the following will go away someday****
  iprtasum = 0
  DO  kdum = 1, 26
     iprtasum = iprtasum + ioptprta(kdum)
  END DO
  WRITE (fuclog,'(10X,A,I2)')  &
       'Number of Variable Arrays Selected to be Printed: ',iprtasum
9015 FORMAT(tr10,a)
9115 FORMAT(tr10,a/tr15,a)
  RETURN
CONTAINS
  SUBROUTINE read_print_control_2
     IMPLICIT NONE
     IF (informt == 1) READ(fuins,'(1X)')
     ichopt = ABS(ichgpr)
     IF (ichgpr == 99) THEN
        READ(fuins,*) prfreq, ioptpr(1), ioptpr(2), ioptpr(3), ioptpr(4), ioptpr(5), ioptpr(6)
     ELSE IF (ichopt == 1) THEN
        READ(fuins,*) prfreq
     ELSE
        READ(fuins,*) prfreq, (ioptpr(idum), idum=1, ichopt-1)
     END IF
     IF (prfreq == 0.0_kdp) prfreq = tchg       ! ... init at 1.0e99_kdp; supress printing
     ! ...       if PRFREQ is a neg number, print results every timestep
     ! ...         beginning with timestep number = -PRFREQ
     IF (prfreq < 0._kdp) THEN
        prtimenxt = tchg
     ELSE
        IF (iyr == 2) prfreq = sec(prfreq)
     END IF
     WRITE(fupdef,'(/A)') '--- Printing Controls ---'
     WRITE(fuclog,'(/A)') '--- Printing Controls ---'
     IF (ichopt >= 2) THEN
        IF (ioptpr(1) == 1) THEN
           WRITE (fuclog, 9015) 'ASCII Explorer plotfile not presently available'
           9015 FORMAT(tr10,a)
        ELSE IF (ioptpr(1) == 2) THEN
           WRITE (fuclog, 9015) 'Binary gnu plotfile written'
        ELSE IF (ioptpr(1) == 3) THEN
           WRITE (fuclog, 9015) 'ASCII gnu plotfile written'
        ELSE IF (ioptpr(1) == 8) THEN
           WRITE (fuclog, 9015) 'ASCII B2 plotfile written'
        ELSE IF (ioptpr(1) == 9) THEN
           WRITE (fuclog, 9015) 'Binary B2 plotfile written'
        END IF
     END IF
     IF (ichopt >= 3) THEN
        IF(ioptpr(2) == 0) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 80 columns'
        ELSE IF(ioptpr(2) == 1) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 132 columns'
        ELSE IF(ioptpr(2) == 2) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 159 columns'
        END IF
     END IF
     IF (ichopt >= 4) THEN
        IF (ioptpr(3) == 0) THEN
           WRITE(fuclog, 9015) 'Minimal output written to Calc_log file and terminal screen'
        ELSE IF (ioptpr(3) == 1) THEN
           WRITE(fuclog, 9015) 'Partial output written to Calc_log file and terminal screen'
        ELSE IF (ioptpr(3) == 2) THEN
           WRITE(fuclog, 9015) 'Complete output written to Calc_log file and terminal screen'
        END IF
     END IF
     IF (ichopt >= 5 .AND. ioptpr(4) == 1) WRITE (fuclog, 9015)  &
          'Vertically averaged values of P,H,T,Sw printed'
     ! ...        read I,J,K indices of nodes to print transient results
     print_temporal = 0
     IF (ichopt >= 7) THEN
        IF (ioptpr(6) > 0) THEN
           n_locations = ioptpr(6)
           print_temporal = prfreq
           DO  k1=1,ioptpr(6)
              READ(fuins,*) (ijktp(k1,k2), k2=1, 3)
              ! ...               check for input errors
              ! ...           check that node values to print are within range
              IF (ijktp(k1,1) > nx .OR. ijktp(k1,2) > ny .OR.  &
                   ijktp(k1,3) > nz) THEN
                 WRITE (fustdout, 9020) (ijktp(k1,k2), k2=1,3)
                 9020             FORMAT (/  &
                  ' *** WARNING, Input error for nodes for time series ***',  &
                  /6X,'Indices of node ', 3I4, ' out of range')
                 WRITE (fuclog, 9020) (ijktp(k1,k2), k2=1,3)
                 IF (ijktp(k1,1) > nx) ijktp(k1,1) = 1
                 IF (ijktp(k1,2) > ny) ijktp(k1,2) = 1
                 IF (ijktp(k1,3) > nz) ijktp(k1,3) = 1
                 WRITE (fustdout, '(5X,A,3I4)') 'Indicies reset to',  &
                      (ijktp(k1,k2), k2=1,3)
                 WRITE (fuclog, '(5X,A,3I4,/)') 'Indicies reset to',  &
                      (ijktp(k1,k2), k2=1,3)
              END IF
           END DO
        END IF
     END IF
     IF (ichgpr < 0) THEN
        !Solely changing print frequency
         IF (print_plotscalar /= 0) print_plotscalar = prfreq
         IF (print_source /= 0) print_source = prfreq
         IF (print_bcflow /= 0) print_bcflow = prfreq
         IF (print_press /= 0) print_press = prfreq
         IF (print_enth /= 0) print_enth = prfreq
         IF (print_temp /= 0) print_temp = prfreq
         IF (print_satn /= 0) print_satn = prfreq
         IF (print_resid /= 0) print_resid = prfreq
         IF (print_dens /= 0) print_dens = prfreq
         IF (print_vis /= 0) print_vis = prfreq
         IF (print_vel /= 0) print_vel = prfreq
         IF (print_pot /= 0) print_pot = prfreq
         IF (print_poros /= 0) print_poros = prfreq
         IF (print_perm /= 0) print_perm = prfreq
         IF (print_pmprop /= 0) print_pmprop = prfreq
     ELSE
        ! ...              ----- arrays to be printed -----
        DO  kdum = 1, 26
           ioptprta(kdum) = 0
        END DO
        IF (informt == 1) READ(fuins, '(1X)')
        READ(fuins,*) (ioptprta(k1), k1=1,6)
        IF (informt == 1) READ(fuins, '(1X)')
        READ(fuins,*) (ioptprta(k1), k1=7,12)
        IF (informt == 1) READ(fuins, '(1X)')
        READ(fuins,*) (ioptprta(k1),k1=13,18)
        ! ...         for input format version 2 - read in print opts for rock parms
        IF (informt == 2) READ(fuins,*) (ioptprta(k1),k1=19,26)
        ! convert to print control parameters used by version 3
        ! plots
        IF (ioptpr(1) == 1 .OR. ioptpr(1) == 2 .OR. ioptpr(1) == 3 .OR. ioptpr(1) == 8 .OR. ioptpr(1) == 9) then
           print_plotscalar = prfreq
        ELSE
           print_plotscalar = 0
        END IF
        ! source/sink and constant-block flux data (printed to HTout5 by version 2.2)
        IF (ioptpr(5) == 1) THEN
           print_source = prfreq
           print_bcflow = prfreq
        ELSE
           print_source = 0
           print_bcflow = 0
        END IF
        iprsource = 0
        iprbcflow = 0           ! not used
        ! pressure
        IF (ioptprta(1) == 1) THEN
           print_press = prfreq
        ELSE
           print_press = 0
        END IF
        ! enthalpy
        IF (ioptprta(2) == 1) THEN
           print_enth = prfreq
        ELSE
           print_enth = 0
        END IF
        ! temperature
        IF (ioptprta(3) == 1) THEN
           print_temp = prfreq
        ELSE
           print_temp = 0
        END IF
        ! saturation
        IF (ioptprta(4) >= 1 .AND. ioptprta(4) <= 3) THEN
           print_satn = prfreq
           iprsat = ioptprta(4)
        ELSE
           print_satn = 0
           iprsat = 0
        END IF
         ! residual mass and energy flux
        IF (ioptprta(5) == 1 .OR. ioptprta(6) == 1) THEN
           print_resid = prfreq
        ELSE
           print_resid = 0
        END IF
        ! water and steam density
        IF (ioptprta(7) == 1 .OR. ioptprta(13) == 1) THEN
           print_dens = prfreq
           IF (ioptprta(13) /= 1) THEN
              iprden = 10                 ! print water density only
           ELSE IF (ioptprta(7) /= 1) THEN
              iprden =  1                 ! print steam density only
           ELSE
              iprden = 11                 ! print both water and steam density
           END IF
        ELSE
           print_dens = 0
           iprden = 0
        END IF
        ! water and steam viscosity
        IF (ioptprta(8) == 1 .OR. ioptprta(14) == 1) THEN
           print_vis = prfreq
           IF (ioptprta(14) /= 1) THEN
              iprvis = 10                 ! print water viscosity only
           ELSE IF (ioptprta(8) /= 1) THEN
              iprvis =  1                 ! print steam viscosity only
           ELSE
              iprvis = 11                 ! print both water and steam viscosity
           END IF
        ELSE
           print_vis = 0
           iprvis = 0
        END IF
        ! water and steam velocity
        lprwater = .FALSE.
        IF (( ioptprta(9) >= 1 .AND.  ioptprta(9) <= 3) .OR.   &
            (ioptprta(10) >= 1 .AND. ioptprta(10) <= 3) .OR.   &
            (ioptprta(11) >= 1 .AND. ioptprta(11) <= 3)) lprwater = .TRUE.
        lprsteam = .FALSE.
        IF ((ioptprta(15) >= 1 .AND. ioptprta(15) <= 3) .OR.   &
            (ioptprta(16) >= 1 .AND. ioptprta(16) <= 3) .OR.   &
            (ioptprta(17) >= 1 .AND. ioptprta(17) <= 3)) lprsteam = .TRUE.
        IF (lprwater .OR. lprsteam) THEN
           print_vel = prfreq
           IF (.NOT. lprsteam) THEN
              iprvel = 30
           ELSE IF (.NOT. lprwater) THEN
              iprvel = 3
           ELSE
              iprvel = 33
           END IF
        ELSE
           print_vel = 0
           iprvel = 0
        END IF
        ! water and steam potential
        IF (ioptprta(12) == 1 .OR. ioptprta(18) == 1) THEN
           print_pot = prfreq
           IF (ioptprta(18) /= 1) THEN
              iprpot = 10                 ! print water hydraulic pressure only
           ELSE IF (ioptprta(12) /= 1) THEN
              iprpot = 1                  ! print steam hydraulic pressure only
           ELSE
              iprpot = 11                 ! print both water and steam hydraulic pressure
           END IF
        ELSE
           print_pot = 0
           iprpot = 0
        END IF
        ! porosity
        IF (ioptprta(19) == 1) THEN
           print_poros = prfreq
        ELSE
           print_poros = 0
        END IF
        ! permeability
        IF (ioptprta(20) == 1 .OR. ioptprta(21) == 1 .OR. ioptprta(22) == 1) THEN
           print_perm = prfreq
        ELSE
           print_perm = 0
        END IF
        ! rock properties: thermal conductivity, specific head, density, compressibility
        IF (ioptprta(23) == 1 .OR. ioptprta(24) == 1 .OR. ioptprta(25) == 1 .OR. ioptprta(26) == 1) THEN
           print_pmprop = prfreq
           IF (ioptprta(25) /= 1 .AND. ioptprta(26) /= 1) THEN
              iprpmprop = 1
           ELSE IF (ioptprta(23) /= 1 .AND. ioptprta(24) /= 1) THEN
              iprpmprop = 10
           ELSE
              iprpmprop = 11
           END IF
        ELSE
           print_pmprop = 0
           iprpmprop = 0
        END IF
        print_balance = 0
        print_plotvector = 0
        print_dump = 1.0e40
        print_dimno = 0
        iprdimno = 0
     END IF
  END SUBROUTINE read_print_control_2
  !
  SUBROUTINE read_print_control_3
     IMPLICIT NONE
     ichopt = ABS(ichgpr)
     ! Note that for data input format VER3 or higher, ichgpr is either 99 or 1.
     ! ichgrp is 99 when printopt is called for the first time
     ! ichgrp is 1 when printopt is call at the start of a new simulation period
     prfreq = 1.e99_kdp 
     IF (ichgpr == 99) THEN
        READ(fuins,*) ioptpr(2), ioptpr(3), ioptpr(4)
        WRITE(fupdef,'(/A)') '--- Printing Controls ---'
        WRITE(fuclog,'(/A)') '--- Printing Controls ---'
        IF (ioptpr(2) == 0) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 80 columns'
           9015 FORMAT(tr10,a)
        ELSE IF(ioptpr(2) == 1) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 132 columns'
        ELSE IF(ioptpr(2) == 2) THEN
           WRITE(fuclog, 9015) 'Printout formatted for 159 columns'
        END IF
        IF (ioptpr(3) == 0) THEN
           WRITE (fuclog, 9015) 'Minimal output written to Calc_log file and terminal screen'
        ELSE IF (ioptpr(3) == 1) THEN
           WRITE (fuclog, 9015) 'Partial output written to Calc_log file and terminal screen'
        ELSE IF (ioptpr(3) == 2) THEN
           WRITE (fuclog, 9015) 'Complete output written to Calc_log file and terminal screen'
        END IF
        IF(ioptpr(4) == 1) WRITE(fuclog,9015)  &
             'Vertically averaged values of P,H,T,Sw printed'
     END IF
     ! ... Read print intervals and flags for arrays to be printed
     ioptprta = 0     ! ... Initialize print option flags
     READ(fuins,*) print_press, print_enth, print_temp, print_satn, iprsat
     READ(fuins,*) print_dens, iprden, print_vis, iprvis, print_pot, iprpot
     READ(fuins,*) print_vel, iprvel, print_bcflow, iprbcflow,  &
          print_source, iprsource
     READ(fuins,*) print_pmprop, iprpmprop, print_poros, print_perm
     READ(fuins,*) print_balance, print_dimno, iprdimno, print_resid, print_dump
     READ(fuins,*) print_plotscalar, print_plotvector, plotfile_type,  &
          print_temporal
     n_locations = 0
     ioptpr(6) = 0
     ioptpr(1) = plotfile_type
     IF (ABS(print_temporal) > 0._kdp) THEN
        READ(fuins,*) n_locations, ((ijktp(k1,k2),k2=1,3),k1=1,n_locations)
        ioptpr(6) = n_locations
        DO  k1 = 1, n_locations
          ! ... Check that node values are within mesh
          IF (ijktp(k1,1) > nx .OR. ijktp(k1,2) > ny .OR.  &
                ijktp(k1,3) > nz) THEN
              WRITE (fustdout, 9020) (ijktp(k1,k2), k2=1,3)
              WRITE (fuclog, 9020) (ijktp(k1,k2), k2=1,3)
9020             FORMAT (/  &
                  ' *** WARNING, Input error for nodes for time series ***',  &
                  /6X,'Indices of node ', 3I4, ' out of range')
              ierr(127) = .TRUE.
          END IF
        END DO
     END IF
     ! ... Set the print control flags
     IF(ABS(print_press) > 0._kdp) ioptprta(1) = 1
     IF(ABS(print_enth) > 0._kdp) ioptprta(2) = 1
     IF(ABS(print_temp) > 0._kdp) ioptprta(3) = 1
     ioptprta(4) = iprsat
     IF(ABS(print_resid) > 0._kdp) THEN
        ioptprta(5) = 1
        ioptprta(6) = 1
     END IF
     IF(nx > 1) ioptprta(9) = iprvel/10
     IF(ny > 1) ioptprta(10) = iprvel/10
     IF(nz > 1) ioptprta(11) = iprvel/10
     IF(nx > 1) ioptprta(15) = MOD(iprvel,10)
     IF(ny > 1) ioptprta(16) = MOD(iprvel,10)
     IF(nz > 1) ioptprta(17) = MOD(iprvel,10)
     ioptprta(7) = iprden/10
     ioptprta(13) = MOD(iprden,10)
     ioptprta(8) = iprvis/10
     ioptprta(14) = MOD(iprvis,10)
     ioptprta(12) = iprpot/10
     ioptprta(18) = MOD(iprpot,10)
     ioptprta(23) = MOD(iprpmprop,10)
     ioptprta(24) = MOD(iprpmprop,10)
     ioptprta(25) = iprpmprop/10
     ioptprta(26) = iprpmprop/10
     ioptprta(27) = iprdimno/10
     ioptprta(28) = MOD(iprdimno,10)
     IF(ABS(print_poros) > 0._kdp) ioptprta(19) = 1
     IF(ABS(print_perm) > 0._kdp) THEN
        IF(nx > 1) ioptprta(20) = 1
        IF(ny > 1) ioptprta(21) = 1
        IF(nz > 1) ioptprta(22) = 1
     END IF
     IF(ioptpr(1) == 0) ioptpr(1) = plotfile_type     ! ... only set once
     IF (ioptpr(1) == 1) THEN
        WRITE (fuclog, 9015) 'ASCII Explorer plotfile written'
     ELSE IF (ioptpr(1) == 2) THEN
        WRITE (fuclog, 9015) 'Binary gnu plotfile written'
     ELSE IF (ioptpr(1) == 3) THEN
        WRITE (fuclog, 9015) 'ASCII gnu plotfile written'
     ELSE IF (ioptpr(1) == 4) THEN
        WRITE (fuclog, 9015) 'Binary IDL plotfile written'
     ELSE IF (ioptpr(1) == 5) THEN
        WRITE (fuclog, 9015) 'ASCII IDL plotfile written'
     ELSE IF (ioptpr(1) == 6) THEN
        WRITE (fuclog, 9015) 'ASCII tab-separated plotfile written'
     ELSE IF (ioptpr(1) == 8) THEN
        WRITE (fuclog, 9015) 'ASCII B2 plotfile written'
     ELSE IF (ioptpr(1) == 9) THEN
        WRITE (fuclog, 9015) 'Binary B2 plotfile written'
     ELSE
        ioptpr(1) = 0
     END IF
     ioptpr(5) = 0
     IF(ABS(print_source) > 0._kdp) ioptpr(5) = 1
     IF(ABS(print_bcflow) > 0._kdp) ioptpr(5) = 1
     IF(iyr == 2) THEN     ! ... convert to internal time units (s)
        IF (print_press > 0._kdp) print_press = sec(print_press)
        IF (print_enth > 0._kdp) print_enth = sec(print_enth)
        IF (print_temp > 0._kdp) print_temp = sec(print_temp)
        IF (print_satn > 0._kdp) print_satn = sec(print_satn)
        IF (print_dens > 0._kdp) print_dens = sec(print_dens)
        IF (print_vis > 0._kdp) print_vis = sec(print_vis)
        IF (print_pot > 0._kdp) print_pot = sec(print_pot)
        IF (print_vel > 0._kdp) print_vel = sec(print_vel)
        IF (print_bcflow > 0._kdp) print_bcflow = sec(print_bcflow)
        IF (print_source > 0._kdp) print_source = sec(print_source)
        IF (print_pmprop > 0._kdp) print_pmprop = sec(print_pmprop)
        IF (print_poros > 0._kdp) print_poros = sec(print_poros)
        IF (print_perm > 0._kdp) print_perm = sec(print_perm)
        IF (print_balance > 0._kdp) print_balance = sec(print_balance)
        IF (print_dump > 0._kdp) print_dump = sec(print_dump)
        IF (print_resid > 0._kdp) print_resid = sec(print_resid)
        IF (print_dimno > 0._kdp) print_dimno = sec(print_dimno)
        IF (print_plotscalar > 0._kdp) print_plotscalar = sec(print_plotscalar)
        IF (print_plotvector > 0._kdp) print_plotvector = sec(print_plotvector)
        IF (print_temporal > 0._kdp) print_temporal = sec(print_temporal)
     END IF
  END SUBROUTINE read_print_control_3
END SUBROUTINE printopt
