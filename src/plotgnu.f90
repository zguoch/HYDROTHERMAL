SUBROUTINE plotgnu
  ! ... Purpose:  To write results for the current time-level on disk for
  ! ...       later plotting using the plotting program gnuplot with the
  ! ...       HTgnu interface.
  ! ...       File Plot_scalar.xxx contains
  ! ...       File Plot_vector.xxx contains X and Z fluid velocity vectors.
  ! ...       pressure, enthalpy, temperature, etc. data for each node.
  ! ... note: IOPTPR(1) =2 for binary plotfile;  =3 for ASCII text plotfile
  USE machine_constants, ONLY: kdp
  USE f_units, ONLY: fupltsca, fupltvec, fustdout, fuclog
  USE control
  USE mesh
  USE variables
  USE units
  IMPLICIT NONE
  CHARACTER(LEN=5) :: unittime
  INTEGER :: i, ic, icdum, ipltkodd, j, k
  ! ...   ipltkodd = 52 for binary plotfile  = 53 for ASCII plotfile;
  INTEGER :: lchar1, lchar2, lchar3
  LOGICAL :: lprint
  REAL(KIND=kdp) :: timedisply
  REAL(KIND=kdp), DIMENSION(:), ALLOCATABLE :: satnw2, dumarry
  INTEGER :: alloc_stat
  INTERFACE
     SUBROUTINE wr(funit, dumarray, ikodd)
       USE machine_constants, ONLY: kdp
       INTEGER, INTENT(IN) :: funit
       REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: dumarray
       INTEGER, INTENT(IN) :: ikodd
     END SUBROUTINE wr
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4.1.4 $//$Date: 2004/08/05 18:34:42 $'
  ALLOCATE(satnw2(nxyz), dumarry(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: plotgnu'
    ierr(185) = .TRUE.
    RETURN
  END IF
  !     ------------------------------------------------------------------
  !...
!!$  !---set counter for "relative time step" for filename suffix
!!$  ! ... ***suspended***
!!$  lpltexpl = lpltexpl + 1
!!$  IF (lpltexpl >= 100) THEN
!!$     IF (lpltexpl <= 110) WRITE(fustdout, '(5X,A)')  &
!!$          '**** Limit of 100 plotfiles exceeded ****'
!!$     WRITE(fuclog, '(10X,A)')  &
!!$          'Limit of 100 plotfiles exceeded, increase limit in plotgnu.f'
!!$     RETURN
!!$  END IF
!!$  !---convert value of LPLTEXPL to a character
!!$  lchar1 = IDINT(DBLE(lpltexpl)/100.0_kdp)
!!$  lchar2 = IDINT(DBLE(lpltexpl-lchar1*100)/10._kdp)
!!$  lchar3 = lpltexpl - (lchar1*100+lchar2*10)
!!$  ufile = CHAR(lchar1+48)//CHAR(lchar2+48)//CHAR(lchar3+48)
  ! ...     set IPLTKODD = 52 for unformated WR, 53 for formatted WR
  ipltkodd = 50 + ioptpr(1)
  ! ...     for first plotfile write header
  IF (lpltexpl == 1) THEN
     IF (ioptpr(1) == 2) THEN
        ! ...      header information for plot files
        WRITE (fupltsca) 'HYDROTHERM  Version:'//version_name
        WRITE (fupltsca) usuff
        WRITE (fupltsca) nx, nz, ny
        WRITE (fupltsca) (xx(i), i=1, nx)
        WRITE (fupltsca) (zz(k), k=1, nz)
        WRITE (fupltsca) (yy(j), j=1, ny)
        WRITE (fupltvec) 'HYDROTHERM  Version:'//version_name
        WRITE (fupltvec) usuff
        WRITE (fupltvec) nx, nz, ny
        WRITE (fupltvec) (xx(i), i=1, nx)
        WRITE (fupltvec) (zz(k), k=1, nz)
        WRITE (fupltvec) (yy(j), j=1, ny)
     ELSE
        WRITE (fupltsca,'(A)') 'HYDROTHERM  Version:'//version_name
        WRITE (fupltsca,'(A)') usuff
        WRITE (fupltsca,'(3I5)') nx, nz, ny
        WRITE (fupltsca,9005) (xx(i), i=1, nx)
        WRITE (fupltsca,9005) (zz(k), k=1, nz)
        WRITE (fupltsca,9005) (yy(j), j=1, ny)
        WRITE (fupltvec,'(A)') 'HYDROTHERM  Version:'//version_name
        WRITE (fupltvec,'(A)') usuff
        WRITE (fupltvec,'(3I5)') nx, nz, ny
        WRITE (fupltvec,9005) (xx(i), i=1, nx)
        WRITE (fupltvec,9005) (zz(k), k=1, nz)
        WRITE (fupltvec,9005) (yy(j), j=1, ny)
9005    FORMAT (14(1PE11.3))
     END IF
  END IF
  ! ... Set value for display of time in seconds or years
  IF (iyr == 2) THEN
     timedisply = year(time)
     unittime = ' (yr)'
  ELSE
     timedisply = time
     unittime = ' (s)'
  END IF
  IF(ABS(print_plotscalar) > 0.) THEN
     CALL print_control(print_plotscalar,time,ltime,tchg,timprpscal,lprint)
     IF(lprint) THEN
        ! ... write the time stamp line 
        IF (ioptpr(1) == 2) THEN
           WRITE (fupltsca) 'Time: ', timedisply, unittime
        ELSE
           WRITE(fupltsca,'(A6,1PG14.7)') 'Time: ', timedisply, unittime
        END IF
        ! ...       ----- pressure -----
        CALL wr(fupltsca, p, ipltkodd)
        ! ...       ----- enthalpy -----
        CALL wr(fupltsca, h, ipltkodd)
        ! ...       ----- temperature -----
        CALL wr(fupltsca, tc, ipltkodd)
        ! ...       ----- water saturation -----
        ! ...     create a modified sat'n array SATNW2 so that all supercrital
        ! ...       nodes (P>Pc) have a saturation of -1 instead of 0.5
        DO  ic = 1,nxyz
           IF (p(ic) > 220.55e6_kdp) THEN
              ! ...          use the following if supercritical includes T>Tc
!!$                IF ( (P(IC).GT.220.55e6_kdp.AND.SATNW(IC).EQ.0.5_kdp).OR.
!!$               &       (TC(IC).GE.374.0D0.AND.SATNW(IC).EQ.0.0D0) ) THEN
              satnw2(ic) = -1.0_kdp
           ELSE
              satnw2(ic) = satnw(ic)
           END IF
        END DO
        CALL wr(fupltsca, satnw2, ipltkodd)
        ! ...       ----- water mass fraction -----
        DO  icdum = 1,npiccons
           ic = npic(icdum)
           dumarry(ic) = satnw(ic)*denw(ic)/  &
                (satnw(ic)*denw(ic)+(1.0_kdp-satnw(ic))*dens(ic))
        END DO
        CALL wr(fupltsca, dumarry, ipltkodd)
     END IF
  END IF
  IF(ABS(print_plotvector) > 0.) THEN
     CALL print_control(print_plotvector,time,ltime,tchg,timprpvec,lprint)
     IF(lprint) THEN
     ! ...       write the time stamp line 
        IF (ioptpr(1) == 2) THEN
           WRITE (fupltsca) 'Time: ', timedisply, unittime
        ELSE
           WRITE(fupltsca,'(A6,1PG14.7)') 'Time: ', timedisply, unittime
        END IF
        ! ...       X velocity of water
        CALL wr(fupltvec, xwvel, ipltkodd)
        ! ...       Z velocity of water
        CALL wr(fupltvec, zwvel, ipltkodd)
        ! ...       X mass flux of water
        CALL wr(fupltvec, xwmflx, ipltkodd)
        ! ...       Z mass flux of water
        CALL wr(fupltvec, zwmflx, ipltkodd)
        ! ...       X velocity of steam
        CALL wr(fupltvec, xsvel, ipltkodd)
        ! ...       Z velocity of steam
        CALL wr(fupltvec, zsvel, ipltkodd)
        ! ...       X mass flux of steam
        CALL wr(fupltvec, xsmflx, ipltkodd)
        ! ...       Z mass flux of steam
        CALL wr(fupltvec, zsmflx, ipltkodd)
     END IF
  END IF
  DEALLOCATE(satnw2, dumarry)
END SUBROUTINE plotgnu
