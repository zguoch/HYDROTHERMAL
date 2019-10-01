SUBROUTINE sinkwell
  ! ... Purpose:  To read source/sink data for wells
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE mesh
  USE source 
  USE steam_table_funct
  USE units
  USE variables
  IMPLICIT NONE
  !
  INTEGER :: i, ic, ij, j, k, kdum, mdum, nlayers, rdstat, iwell
  INTEGER, DIMENSION(:), ALLOCATABLE :: ktmp
  INTEGER :: alloc_stat
  REAL(KIND=kdp) :: dum1, qq, qqh, qqp, qqt
  !
  ! ...         I,J - Map view location of well
  ! ...         QQ = QT(IJ) - total strength of well
  ! ...                (negative for production, positive for injection)
  ! ...         QQH = QHI(IJ) - enthalpy of injection fluid
  ! ...         NLAYERS - number of open intervals (cells) for this well
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.10 $//$Date: 2006/12/28 21:10:33 $'
  ALLOCATE(ktmp(nz), STAT=alloc_stat)
  IF (alloc_stat /= 0) THEN
     PRINT *, 'ARRAY ALLOCATION FAILED: sink'
     ierr(178) = .TRUE.
     RETURN
  END IF
  !     ------------------------------------------------------------------
  IF(informt == 1) READ(fuins,*)
  i = -1
  j = -1
  iwell = 0
  DO
     READ(fuins,*,IOSTAT=rdstat) i, j, qq, qqh, nlayers
     IF(rdstat /= 0) THEN
        WRITE (fustdout,'(a)') 'ERROR - Problem with WELL input for well source/sink '
        WRITE (fuclog,'(a)') 'ERROR - Problem with WELL input for well source/sink '
        ierr(114) = .TRUE.
     END IF
     IF(i == 0 .AND. j == 0) EXIT     ! ... Terminate well source input
     IF(irad .AND. j == 0) j = 1
     ! ... check for input errors
     IF (qq < 0._kdp .AND. ABS(qqh) > 0._kdp) THEN
!!$     WRITE(fuclog, '(a)')  &
!!$          'ERROR - Input enthalpy of fluid must be 0 for discharging wells'
!!$     ierr(107) = .TRUE.
        qqh = 0._kdp          ! ... set to zero
     END IF
     IF (i <= 0 .OR. i > nx .OR. j<= 0 .OR. j > ny .OR. nlayers <= 0 .OR. nlayers > nz) THEN
        WRITE(fustdout, 9215) 'ERROR - Invalid i,j indices or nlayers of well: ',  &
             'i, j, nlayers: ',i,j,nlayers
        WRITE(fuclog, 9215) 'ERROR - Invalid i,j indices or nlayers of well: ',  &
             'i, j, nlayers: ',i,j,nlayers
9215    FORMAT(a/tr5,a,3i5)
        ierr(108) = .TRUE.
     END IF
     WRITE(fupdef,'(/tr18,a,i4,a,i4)') 'Well location indices    i: ',i,', j: ',j
     ij = (j-1)*nx + i
     IF(unitlabel(17) == '(Deg.C)      ') THEN   ! ... temperature input, calculate enthalpy
        qqt = qqh
        qqp = 225.0e6_kdp
        ! ... first guess at pressure
        IF (qqt < 95._kdp) THEN
           qqp = 1.0e6_kdp
        ELSE IF (qqt < 170._kdp) THEN
           qqp = 10.0e6_kdp
        ELSE IF (qqt < 260._kdp) THEN
           qqp = 50.0e6_kdp
        ELSE IF (qqt < 300._kdp) THEN
           qqp = 100.e6_kdp
        ELSE
           qqp = 225.e6_kdp
        END IF
        CALL tempress(qqt, qqp, dum1, qqh)
        qhi(ij) = qqh
     ELSE               ! ... enthalpy input
        qqt = 0._kdp
        qhi(ij) = qqh*unitfac(17)    ! ... Convert enthalpy from input units to cgs units
     END IF
     qt(ij) = qq*unitfac(16)     ! ... Convert flow rate from input units to cgs units
     IF (qq < 0._kdp) THEN          ! ... production well
        WRITE(fupdef,9305) 'Production Well Flow Rate ......... ',qq,unitlabel(16)
9305    FORMAT (tr18,a,1pg12.4,tr1,a)
     ELSE                           ! ... injection well
        WRITE(fupdef,9305) 'Injection Well Flow Rate .......... ',qq,unitlabel(16)
        IF (qqt == 0._kdp) THEN     ! ... enthalpy input
           WRITE(fupdef,9005) 'Injection Fluid Enthalpy .......... ',qqh,unitlabel(17)
        ELSE                        ! ... temperature input
           WRITE(fupdef,9005) 'Injection Fluid Temperature ....... ',qqt,'(Deg.C)'
           WRITE(fupdef,9005) 'Calculated Injection Fluid Enthalpy ',  &
                qqh/unitfac(10),unitlabel(10)
9005       FORMAT (tr18,a,1pg12.4,tr1,a)
        END IF
     END IF
     ktmp = 0         ! ... re-initialize before reading in K's for open intervals
     ! ... Read indices of layers over which flow of well is applied
     READ(fuins,*) (ktmp(kdum),kdum=1,nlayers)
     WRITE(fupdef,'(tr23,a,9i4/(tr46,9i4))') 'Well open interval indices   k: ',  &
          (ktmp(kdum),kdum=1,nlayers)
     ! ...         assign IQ(IC); 1 - open layers, 0 - cased layers
     DO  mdum = 1,nlayers
        k = ktmp(mdum)
        ic = (k-1)*nx + i + (j-1)*nxz
        ! ... check for input errors
        IF(k <= 0 .OR. k > nz) THEN
           WRITE (fuclog,9215) 'ERROR - Index for well open interval out of range ',  &
                'k: ',k
           ierr(109) = .TRUE.
        ELSEIF(ABS(qnodsrc(ic)) > 0._kdp) THEN
           WRITE(fustdout,9020) 'ERROR - Well source in a cell column with a point source',  &
                'i,j,k = ',i,j,k
           WRITE(fuclog,9020) 'ERROR - Well source in a cell column with a point source',  &
                'i,j,k = ',i,j,k
9020       FORMAT (/a/tr10,a,3i4)
           ierr(113) = .TRUE.
        ELSE
           iq(ic) = 1      ! ...    communicate the source to the cell
        END IF
     END DO
     iwell = iwell + 1
     IF (informt<=2 .AND. iwell==nsrce) EXIT
  END DO
  DEALLOCATE(ktmp)
END SUBROUTINE sinkwell
