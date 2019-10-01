SUBROUTINE sinknode
  ! ... Purpose:  To read point source/sink data for nodes
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE mesh
  USE source
  USE steam_table_funct
  USE units
  IMPLICIT NONE
  !
  INTEGER :: i, ic, ij, j, k, rdstat
  REAL(KIND=kdp) :: dum1, qq, qqh, qqp, qqt
  !
  ! ...         I,J,K - Location of point source
  ! ...         QQ - total strength of source 
  ! ...                (negative for production, positive for injection)
  ! ...         QQH - enthalpy of injection fluid
  !
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.12 $//$Date: 2006/12/28 21:10:33 $'
  !     ------------------------------------------------------------------
  i = -1
  j = -1
  k = -1
  DO 
     READ(fuins,*,IOSTAT=rdstat) i, j, k, qq, qqh
     IF(rdstat /= 0) THEN
        WRITE (fustdout,'(a)') 'ERROR - Problem with NODE input for point source/sink '
        WRITE (fuclog,'(a)') 'ERROR - Problem with NODE input for point source/sink '
        ierr(106) = .TRUE.
     END IF
     IF(i == 0 .AND. j == 0 .AND. k == 0) EXIT     ! ... Terminate point source input
     IF(irad .AND. j == 0) j = 1
     ! ... check for input errors
     IF (qq < 0._kdp .AND. ABS(qqh) > 0._kdp) THEN
!!$        WRITE(fuclog,'(A)') 'Input enthalpy of fluid will be ignored for sink nodes'
!!$        ierr(104) = .TRUE.
        qqh = 0._kdp          ! ... set to zero
     END IF
     IF (i <= 0 .OR. i > nx .OR. j<= 0 .OR. j > ny .OR. k <= 0 .OR. k > nz) THEN
        WRITE(fustdout,9020) 'ERROR - Invalid i,j,k indices of source/sink node: ',  &
             'i,j,k: ',i,j,k
        WRITE(fuclog,9020) 'ERROR - Invalid i,j,k indices of source/sink node: ',  &
             'i,j,k: ',i,j,k
9020    FORMAT (a/tr5,a,3i5)
        ierr(105) = .TRUE.
     END IF
     IF(.NOT.irad) THEN
        WRITE(fupdef,'(/tr18,3(a,i4))') 'Point source location indices    i:',i,',  j:',j,',  k:',k
     ELSEIF(irad) THEN
        WRITE(fupdef,'(/tr18,2(a,i4))') 'Point source location indices    i:',i,',  k:',k
     END IF
     ij = (j-1)*nx + i
     ic = (k-1)*nx + i + (j-1)*nxz
     IF(unitlabel(17) == '(Deg.C)      ') THEN     ! ... temperature input, calculate enthalpy
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
        qhinodsrc(ic) = qqh
     ELSE                     ! ... enthalpy input
        qqt = 0._kdp
        qhinodsrc(ic) = qqh*unitfac(17)   ! ... Convert enthalpy from input units to cgs units
     END IF
     qnodsrc(ic) = qq*unitfac(16)   ! ... Convert flow rate from input units to cgs units
     IF (qq < 0._kdp) THEN       ! ... point sink node
        WRITE(fupdef,9305) 'Point Sink Node Flow Rate ......... ',qq,unitlabel(16)
9305    FORMAT (tr18,a,1pg12.4,tr1,a)
     ELSE                        ! ... point source node
        WRITE(fupdef,9305) 'Point Source Node Flow Rate ......... ',qq,unitlabel(16)
        IF (qqt == 0._kdp) THEN     ! ... enthalpy input
           WRITE(fupdef,9105) 'Point Source Fluid Enthalpy .......... ',qqh,unitlabel(17)
        ELSE                        ! ... temperature input
           WRITE(fupdef,9105) 'Point Source Fluid Temperature ....... ',qqt,'(Deg.C)'
           WRITE(fupdef,9105) 'Calculated Source Fluid Enthalpy ',  &
                qqh/unitfac(10),unitlabel(10)
9105       FORMAT (tr18,a,1pg12.4,tr1,a)
        END IF
     END IF
     IF(.NOT.irad .AND. ABS(qt(ij)) > 0._kdp) THEN
        WRITE(fustdout,9120) 'ERROR - Point source node in a cell column with a well',  &
             'i,j,k = ',i,j,k
        WRITE(fuclog,9120) 'ERROR - Point source node in a cell column with a well',  &
             'i,j,k = ',i,j,k
9120    FORMAT (/a/tr10,a,3i4)
        ierr(112) = .TRUE.
     ELSE
        iq(ic) = 1      ! ...    communicate the source to the cell
     END IF
  END DO
END SUBROUTINE sinknode
