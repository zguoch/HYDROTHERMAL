SUBROUTINE plotxyz
  ! ... Purpose:  To write results for the current time-level on disk for
  ! ...       later plotting. Format is tab separated columns of x,y,z,t and scalar
  ! ...       variables. Independent variables in user units.
  ! ...       File Plot_scalar contains pressure, temperature, 
  ! ...            saturation, phase index, Nuselt number.
  ! ...       File Plot_scalar2 contains water-table elevation
  ! ...       File Plot_vector contains vector components of water and steam fluxes at nodes
  ! ...            (cgs units) 
  ! ..
  ! ... note: IOPTPR(1) = 6 for this columnar ASCII text plotfile
  USE machine_constants, ONLY: kdp
  USE f_units, ONLY: fupltsca, fupltsca2, fupltvec
  USE bc, ONLY: ibc, unconfined
  USE control
  USE mesh
  USE variables
  USE units
  IMPLICIT NONE
  CHARACTER(LEN=80) :: ufile
  CHARACTER(LEN=4) :: unitlabel_time
  INTEGER :: i, ic, ik, j, k
  LOGICAL :: lprint
  INTEGER, DIMENSION(:), ALLOCATABLE :: indprt
  REAL(kind=kdp), DIMENSION(:), ALLOCATABLE :: satnw2
  INTEGER :: alloc_stat
!!$  REAL(kind=kdp), DIMENSION(nxy) :: aprint
  REAL(kind=kdp) :: cnvtmi
  !
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.6 $//$Date: 2006/02/09 21:40:31 $'
  ALLOCATE(indprt(nxyz), satnw2(nxyz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
    PRINT *, 'ARRAY ALLOCATION FAILED: plotxyz'
    ierr(183) = .TRUE.
    RETURN
  END IF    
  !     ------------------------------------------------------------------
  IF(iyr == 1) THEN
     cnvtmi = 1._kdp
     unitlabel_time = '(s)'
  ELSEIF(iyr == 2) THEN
     cnvtmi = 3.16880878e-8_kdp 
     unitlabel_time = '(yr)'
  ENDIF
  ! ... For first access of plotfiles, write headers
  lpltexpl = lpltexpl + 1
  IF (lpltexpl == 1) THEN
     ! ... formatted output only
     WRITE(fupltsca,'(A)') title1
     WRITE(fupltsca,'(A)') title2
     WRITE(fupltsca,'(A)') usuff
     ! ... Write header to file 'fupltsca' for xyz scalar field data
     WRITE(fupltsca,5011) 'x'//achar(9)//'y'//achar(9)//'z'//achar(9)//'time'//  &
          achar(9)//'temperature'//achar(9)//'pressure'//achar(9)//'saturation'//  &
          achar(9)//'phase index'//achar(9)//'Cell Nusselt No.'//achar(9) 
     WRITE(fupltsca,5011) unitlabel(12)//achar(9)//unitlabel(13)//achar(9)//  &
          unitlabel(14)//achar(9)//unitlabel_time//  &
          achar(9)//unitlabel(15)//achar(9)//unitlabel(11)//achar(9)//'(-)'//  &
          achar(9)//' '//achar(9)//'(-)'//achar(9) 
     IF(unconfined) THEN
        WRITE(fupltsca2,'(A)') title1
        WRITE(fupltsca2,'(A)') title2
        WRITE(fupltsca2,'(A)') usuff
        ! ... Write header to file 'fupltsca2' for xyz scalar field data
        WRITE(fupltsca2,5011) 'x'//achar(9)//'y'//achar(9)//'time'//achar(9)//  &
             'water-table elevation'//achar(9)
        WRITE(fupltsca2,5011) unitlabel(12)//achar(9)//unitlabel(13)//achar(9)//  &
             unitlabel_time//achar(9)//unitlabel(14)//achar(9)
     END IF
     WRITE(fupltvec,'(A)') title1
     WRITE(fupltvec,'(A)') title2
     WRITE(fupltvec,'(A)') usuff
     ! ... Write header to file 'fupltvec' for xyz vector field data
     WRITE(fupltvec,5011) 'x'//achar(9)//'y'//achar(9)//'z'//achar(9)//'time'//  &
          achar(9)//  &
          'x water mass flux'//achar(9)//'y water mass flux'//achar(9)//  &
          'z water mass flux'//achar(9)//  &
          'x steam mass flux'//achar(9)//'y steam mass flux'//achar(9)//  &
          'z steam mass flux'//achar(9)
     WRITE(fupltvec,5011) unitlabel(12)//achar(9)//unitlabel(13)//achar(9)//  &
          unitlabel(14)//achar(9)//unitlabel_time//achar(9)//  &
          '(g/s-cm^2)'//achar(9)//'(g/s-cm^2)'//achar(9)//'(g/s-cm^2)'//achar(9)//  &
          '(g/s-cm^2)'//achar(9)//'(g/s-cm^2)'//achar(9)//'(g/s-cm^2)'//achar(9)
5011 format(tr5,a)
  END IF
  IF(ABS(print_plotscalar) > 0._kdp) THEN
     CALL print_control(print_plotscalar,time,ltime,tchg,timprpscal,lprint)
     IF(lprint) THEN
        ! ... Compute and print Peclet and Nusselt numbers for all active cells, if necessary
        IF(ABS(print_dimno) <= 0._kdp) CALL dimnum     ! ... if not already computed
        DO ic=1,nxyz
           CALL ictoijk(ic,i,j,k,nx,nz)
           ik = (k-1)*nx + i
           IF(ibc(ik,j) /= -1) THEN
              ! ... create a modified saturation array SATNW2 so that all supercritical
              ! ...     nodes (P>Pc) have a saturation of -1 instead of 0.5
              IF (ind(ic) == 4) THEN
                 satnw2(ic) = -1.0_kdp
              ELSE
                 satnw2(ic) = satnw(ic)
              END IF
              ! ... change phase index of air-water from 5 to 0 for printout
              IF (ind(ic) == 5) THEN
                 indprt(ic) = 0
              ELSE
                 indprt(ic) = ind(ic)
              END IF
      WRITE(fupltsca,5003) xx(i)/unitfac(12),achar(9),yy(j)/unitfac(13),achar(9),  &
                   zz(k)/unitfac(14),achar(9),  &
                   cnvtmi*time,achar(9),  &
                   tc(ic),achar(9),p(ic)/unitfac(11),achar(9),satnw2(ic),achar(9),  &
                   indprt(ic),achar(9),nu(ic),achar(9)
5003          FORMAT(7(1pg14.6,a),i4,a,1pg14.6,a)
           ENDIF
        END DO
        IF(unconfined) THEN
           DO ic=1,nxyz
              CALL ictoijk(ic,i,j,k,nx,nz)
              IF(abs(fs_elev(ic)) > 0._kdp) THEN
                 WRITE(fupltsca2,5004) xx(i)/unitfac(12),achar(9),yy(j)/unitfac(13),achar(9),  &
                      cnvtmi*time,achar(9),  &
                      fs_elev(ic)/unitfac(14)
5004             FORMAT(4(1pg14.6,a))
              ENDIF
           END DO
        END IF
     END IF
  END IF
  IF(ABS(print_plotvector) > 0._kdp) THEN
     CALL print_control(print_plotvector,time,ltime,tchg,timprpvec,lprint)
     IF(lprint) THEN
        DO ic=1,nxyz
           CALL ictoijk(ic,i,j,k,nx,nz)
           ik = (k-1)*nx + i 
           IF(ibc(ik,j) == -1) CYCLE
           WRITE(fupltvec,5005) xx(i)/unitfac(12),achar(9),yy(j)/unitfac(13),achar(9),  &
                zz(k)/unitfac(14),achar(9),  &
                cnvtmi*time,achar(9),  &
                xwmflx(ic),achar(9),ywmflx(ic),achar(9),zwmflx(ic),achar(9),  &
                xsmflx(ic),achar(9),ysmflx(ic),achar(9),zsmflx(ic),achar(9)
5005       FORMAT(10(1pg14.6,a))
        END DO
     END IF
  END IF
  DEALLOCATE(indprt, satnw2)
END SUBROUTINE plotxyz
