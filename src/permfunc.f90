SUBROUTINE permfunc(dumarray,ndim,kodparm)
  !     Purpose:  To compute the permeability of selected rock types as a
  !               function of time.
  USE machine_constants, ONLY: kdp
  USE mesh
  USE parameters
  USE variables, ONLY: time
  IMPLICIT NONE
  REAL(kind=kdp), DIMENSION(ndim), INTENT(IN OUT) :: dumarray
  INTEGER, INTENT(IN) :: ndim
  INTEGER, INTENT(IN) :: kodparm
  !       DUMARRAY - array that is to be assigned values that are f(time)
  !       KODPARM  - Rock Parameter code
  !             =2 - X permeability
  !             =3 - Y permeability
  !             =4 - Z permeability
  !
  INTEGER :: ic, icdum, irx, lll
  REAL(kind=kdp) :: slope
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.1 $//$Date: 2002/03/18 23:38:32 $'
  !     ------------------------------------------------------------------
  loop30:  DO  icdum = 1, npiccons
     ic = npic(icdum)
     irx = icrxtype(ic)
     !---if all values are to be set
     IF (irxftimopt(kodparm,irx,1) == 1) THEN     ! ... option 1 - NOT a function of time
        dumarray(ic) = rxftimparm(kodparm,irx,1)
     ELSE IF (irxftimopt(kodparm,irx,1) == 2) THEN   ! ... option 2 - linear function of time
        IF (time < rxftimparm(kodparm,irx,2)) THEN    ! ...  time < first time data point
           dumarray(ic) = rxftimparm(kodparm,irx,1)
        ELSE IF (time >= rxftimparm(kodparm,irx,2*irxftimopt(kodparm,irx,2))) THEN
           ! ...  time > last time data point
           dumarray(ic) = rxftimparm(kodparm,irx,2*irxftimopt(kodparm,irx,2)-1)
        ELSE
           DO  lll=1,irxftimopt(kodparm,irx,2)-1
              IF (time >= rxftimparm(kodparm,irx,2*lll) .AND.  &
                   time < rxftimparm(kodparm,irx,2*lll+2)) THEN
                 slope = (rxftimparm(kodparm,irx,2*lll+1)  &
                      -rxftimparm(kodparm,irx,2*lll-1))/(rxftimparm(kodparm,irx,2*lll+2)  &
                      -rxftimparm(kodparm,irx,2*lll))
                 dumarray(ic) = rxftimparm(kodparm,irx,2*lll-1) + slope*(time  &
                      -rxftimparm(kodparm,irx,2*lll))
              END IF
              CYCLE loop30
           END DO
        END IF
     END IF
  END DO loop30
END SUBROUTINE permfunc
