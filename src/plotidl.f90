SUBROUTINE plotidl
  ! ... Purpose:  To write results for the current time-level on disk for
  ! ...       later plotting using the plotting program IDL with the
  ! ...       HTView IDL widget viewer interface.
  ! ...       File Plot_scalar contains pressure, enthalpy, temperature, 
  ! ...            saturation, water mass fraction data.
  ! ...       File Plot_vector contains X, Y, and Z fluid velocity vectors
  ! ...            and mass flux vectors.
  ! ... note: IOPTPR(1) =4 for binary plotfile;  =5 for ASCII text plotfile;
  USE machine_constants, ONLY: kdp
  USE f_units, ONLY: fupltsca, fupltsca2, fupltvec
  USE bc, ONLY: ibc, unconfined
  USE control
  USE mesh
  USE variables
  USE units
  IMPLICIT NONE
  INTEGER :: i, ic, ik, icdum, ipltkodd, j, k, ms
  !  INTEGER :: lchar1, lchar2, lchar3
  LOGICAL :: lprint
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: satnw2, dumarry
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: aprint
  INTEGER :: alloc_stat
  INTERFACE
     SUBROUTINE wr(funit, dumarray, ikodd)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: ikodd
     END SUBROUTINE wr
     SUBROUTINE wrnew(funit,dumarray,kodd,convfac)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: kodd
       REAL(kind=kdp), INTENT(IN OUT) :: convfac
     END SUBROUTINE wrnew
  END INTERFACE
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.5.1.9 $//$Date: 2007/04/10 23:14:34 $'
  ALLOCATE(satnw2(nxyz), dumarry(nxyz), aprint(nxy), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: plotidl'
    ierr(184) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  ipltkodd = 50 + ioptpr(1)  ! ... 54 for unformated write, 55 for formatted write
  ! ...     for first plotfile open unit and write header
  lpltexpl = lpltexpl + 1
  IF (lpltexpl == 1) THEN
     IF (ioptpr(1) == 4) THEN
        ! ...      header information for plot files
        ! ... unformatted
        WRITE(fupltsca) title1
        WRITE(fupltsca) title2
        WRITE(fupltsca) usuff
        WRITE(fupltsca) nx, nz, ny, nxy, nxz, nxyz
        WRITE(fupltsca) (xx(i)/unitfac(12), i=1,nx)
        WRITE(fupltsca) (yy(j)/unitfac(13), j=1,ny)
        WRITE(fupltsca) (zz(k)/unitfac(14), k=1,nz)
        WRITE(fupltsca) ((ibc(ik,j),ik=1,nxz),j=1,ny)
        IF(unconfined) THEN
           WRITE(fupltsca2) title1
           WRITE(fupltsca2) title2
           WRITE(fupltsca2) usuff
           WRITE(fupltsca2) nx, nz, ny, nxy, nxz, nxyz
           WRITE(fupltsca2) (xx(i)/unitfac(12), i=1,nx)
           WRITE(fupltsca2) (yy(j)/unitfac(13), j=1,ny)
           WRITE(fupltsca2) (zz(k)/unitfac(14), k=1,nz)
           WRITE(fupltsca2) ((ibc(ik,j),ik=1,nxz),j=1,ny)
        ENDIF
        WRITE(fupltvec) title1
        WRITE(fupltvec) title2
        WRITE(fupltvec) usuff
        WRITE(fupltvec) nx, nz, ny, nxy, nxz, nxyz
        WRITE(fupltvec) (xx(i)/unitfac(12), i=1,nx)
        WRITE(fupltvec) (yy(j)/unitfac(13), j=1,ny)
        WRITE(fupltvec) (zz(k)/unitfac(14), k=1,nz)
        WRITE(fupltvec) ((ibc(ik,j),ik=1,nxz),j=1,ny)
     ELSE
        ! ... formatted
        WRITE(fupltsca,'(A)') title1
        WRITE(fupltsca,'(A)') title2
        WRITE(fupltsca,'(A)') usuff
        WRITE(fupltsca,'(6I6)') nx, nz, ny, nxy, nxz, nxyz
        WRITE(fupltsca,9005) (xx(i)/unitfac(12), i=1,nx)
        WRITE(fupltsca,9005) (yy(j)/unitfac(13), j=1,ny)
        WRITE(fupltsca,9005) (zz(k)/unitfac(14), k=1,nz)
        WRITE(fupltsca,9006) ((ibc(ik,j),ik=1,nxz),j=1,ny)
        IF(unconfined) THEN
           WRITE(fupltsca2,'(A)') title1
           WRITE(fupltsca2,'(A)') title2
           WRITE(fupltsca2,'(A)') usuff
           WRITE(fupltsca2,'(6I6)') nx, nz, ny, nxy, nxz, nxyz
           WRITE(fupltsca2,9005) (xx(i)/unitfac(12), i=1,nx)
           WRITE(fupltsca2,9005) (yy(j)/unitfac(13), j=1,ny)
           WRITE(fupltsca2,9005) (zz(k)/unitfac(14), k=1,nz)
           WRITE(fupltsca2,9006) ((ibc(ik,j),ik=1,nxz),j=1,ny)
        ENDIF
        WRITE(fupltvec,'(A)') title1
        WRITE(fupltvec,'(A)') title2
        WRITE(fupltvec,'(A)') usuff
        WRITE(fupltvec,'(6I6)') nx, nz, ny, nxy, nxz, nxyz
        WRITE(fupltvec,9005) (xx(i)/unitfac(12), i=1,nx)
        WRITE(fupltvec,9005) (yy(j)/unitfac(13), j=1,ny)
        WRITE(fupltvec,9005) (zz(k)/unitfac(14), k=1,nz) 
        WRITE(fupltvec,9006) ((ibc(ik,j),ik=1,nxz),j=1,ny)
9005    FORMAT (14(1pe11.3))
9006    FORMAT (14I6)
     END IF
  END IF
  IF(ABS(print_plotscalar) > 0._kdp) THEN
     CALL print_control(print_plotscalar,time,ltime,tchg,timprpscal,lprint)
     IF(lprint) THEN
        ! ... Write the header line for each time step for the plot output file
        IF (ioptpr(1) == 4) THEN
           IF (iyr == 2) THEN
              WRITE(fupltsca) 'Time',year(time),' (yr)'
           ELSE
              WRITE(fupltsca) 'Time',time,' (s)'
           END IF
        ELSE
           IF (iyr == 2) THEN
              WRITE(fupltsca,9015) 'Time',year(time),' (yr)'
9015          FORMAT (tr20,a,tr2,1PG14.7,a)
           ELSE
              WRITE(fupltsca,9015) 'Time',time,' (s)'
           END IF
        END IF
        ! ...       ----- pressure -----
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltsca) 'Pressure', unitlabel(11)
        ELSE
           WRITE(fupltsca,'(TR20,2A)') 'Pressure', unitlabel(11)
        END IF
        CALL wrnew(fupltsca, p, ipltkodd, unitfac(11))
        ! ...       ----- enthalpy -----
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltsca) 'Enthalpy', unitlabel(10)
        ELSE
           WRITE(fupltsca,'(TR20,2A)') 'Enthalpy', unitlabel(10)
        END IF
        CALL wrnew(fupltsca, h, ipltkodd, unitfac(10))
        ! ...       ----- temperature -----
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltsca) 'Temperature', unitlabel(15)
        ELSE
           WRITE(fupltsca,'(TR20,2A)') 'Temperature', unitlabel(15)
        END IF
        CALL wrnew(fupltsca, tc, ipltkodd, unitfac(15))
        ! ...       ----- water saturation -----
        ! ... create a modified saturation array SATNW2 so that all supercritical
        ! ...     nodes (P>Pc) have a saturation of -1 instead of 0.5
        DO  ic = 1, nxyz
           IF (p(ic) > 220.55e6_kdp) THEN
              ! ...          use the following if supercritical includes T>Tc
!!$              ! ...   IF ( (P(IC).GT.220.55D6.AND.SATNW(IC).EQ.0.5D0).OR.
!!$              !    &       (TC(IC).GE.374.0D0.AND.SATNW(IC).EQ.0.0D0) ) THEN
              satnw2(ic) = -1.0_kdp
           ELSE
              satnw2(ic) = satnw(ic)
           END IF
        END DO
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltsca) 'Water Saturation', '(volumetric fraction)'
        ELSE
           WRITE(fupltsca,'(TR20,2A)') 'Water Saturation', '(volumetric fraction)'
        END IF
        CALL wr(fupltsca, satnw2, ipltkodd)
        ! ...       ----- water mass fraction -----
        DO  icdum = 1, npiccons
           ic = npic(icdum)
           dumarry(ic) = satnw(ic)*denw(ic)/  &
                (satnw(ic)*denw(ic)+(1._kdp-satnw(ic))*dens(ic))
        END DO
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltsca) 'Mass Fraction Water', '(-)'
        ELSE
           WRITE(fupltsca,'(TR20,2A)') 'Mass Fraction Water', '(-)'
        END IF
        CALL wr(fupltsca, dumarry, ipltkodd)
        IF(unconfined) THEN
           ! ... Load the print array with water-table elevations
           aprint = 0._kdp
           DO  j=1,ny
              DO  i=1,nx
                 DO  k=1,nz
                    ic = (j-1)*nxz + (k-1)*nx + i
                    ms =  (j-1)*nx + i
                    IF(fs_elev(ic) > 0._kdp) THEN
                       aprint(ms) = fs_elev(ic)/unitfac(14)
                    END IF
                 END DO
              END DO
           END DO
           ! ... Write the header line for each time step for the second plot output file
           IF (ioptpr(1) == 4) THEN
              IF (iyr == 2) THEN
                 WRITE(fupltsca2) 'Time',year(time),' (yr)'
              ELSE
                 WRITE(fupltsca2) 'Time',time,' (s)'
              END IF
           ELSE
              IF (iyr == 2) THEN
                 WRITE(fupltsca2,9015) 'Time',year(time),' (yr)'
              ELSE
                 WRITE(fupltsca2,9015) 'Time',time,' (s)'
              END IF
           END IF
           ! ... water-table elevation
           IF (ioptpr(1) == 4) THEN
              WRITE(fupltsca2) 'Water-Table Elevation', unitlabel(14)
              WRITE(fupltsca2) aprint
           ELSE
              WRITE(fupltsca2,'(TR20,2A)') 'Water-Table Elevation', unitlabel(14)
              WRITE(fupltsca2,8001) aprint
8001          FORMAT(9(1pe14.6))
           END IF
        END IF
     END IF
  END IF
  IF(ABS(print_plotvector) > 0._kdp) THEN
     CALL print_control(print_plotvector,time,ltime,tchg,timprpvec,lprint)
     IF(lprint) THEN
        !.....Write the header line for each time step for the plot output file
        IF (ioptpr(1) == 4) THEN
           IF (iyr == 2) THEN
              WRITE(fupltvec) 'Time',year(time),' (yr)'
           ELSE
              WRITE(fupltvec) 'Time',time,' (s)'
           END IF
        ELSE
           IF (iyr == 2) THEN
              WRITE(fupltvec,9015) 'Time',year(time),' (yr)'
           ELSE
              WRITE(fupltvec,9015) 'Time',time,' (s)'
           END IF
        END IF
        ! ...       X velocity of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Velocity X-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Velocity X-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, xwvel, ipltkodd)
        ! ...       Y velocity of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Velocity Y-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Velocity Y-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, ywvel, ipltkodd)
        ! ...       Z velocity of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Velocity Z-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Velocity Z-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, zwvel, ipltkodd)
        ! ...       X mass flux of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Mass Flux X-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Mass Flux X-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, xwmflx, ipltkodd)
        ! ...       Y mass flux of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Mass Flux Y-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Mass Flux Y-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, ywmflx, ipltkodd)
        ! ...       Z mass flux of water
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Water Mass Flux Z-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Water Mass Flux Z-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, zwmflx, ipltkodd)
        ! ...       X velocity of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Velocity X-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Velocity X-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, xsvel, ipltkodd)
        ! ...       Y velocity of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Velocity Y-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Velocity Y-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, ysvel, ipltkodd)
        ! ...       Z velocity of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Velocity Z-direction','(cm/s)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Velocity Z-direction','(cm/s)'
        END IF
        CALL wr(fupltvec, zsvel, ipltkodd)
        ! ...       X mass flux of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Mass Flux X-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Mass Flux X-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, xsmflx, ipltkodd)
        ! ...       Y mass flux of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Mass Flux Y-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Mass Flux Y-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, ysmflx, ipltkodd)
        ! ...       Z mass flux of steam
        IF (ioptpr(1) == 4) THEN
           WRITE(fupltvec) 'Steam Mass Flux Z-direction', '(g/s-cm^2)'
        ELSE
           WRITE(fupltvec,'(TR20,2A)') 'Steam Mass Flux Z-direction', '(g/s-cm^2)'
        END IF
        CALL wr(fupltvec, zsmflx, ipltkodd)
     END IF
  END IF
  DEALLOCATE(satnw2, dumarry, aprint)
END SUBROUTINE plotidl
