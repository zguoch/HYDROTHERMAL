SUBROUTINE ldiaja
  ! ... Loads the IA and JA pointer arrays for the A matrix
  ! ... Used to generate the compressed row storage
  USE control, ONLY: ierr
  USE mesh
  USE linked_list_mod
  USE solver_gmres
  IMPLICIT NONE
  ! 
  INTEGER :: a_err, ic, ma
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2003/01/22 23:24:16 $'
  !     ------------------------------------------------------------------
  !...
  ALLOCATE(ia(2*mrnomax+1), ib(2*mrnomax+1),  &
       STAT=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: ldiaja"
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  CALL init_list_is
  ! ... Loop over all mesh nodes; Assemble for each active node into linked list
  DO ic=1,nxyz
     ma = mrno(ic)
     IF(ma == 0) CYCLE
     ! ...  Assembly of active cell; pressure
     !    ******* Main Diagonal Block Terms *********
     ! ... main diagonal stored first in each row
     CALL add_list_is(2*ma-1)
     ia(2*ma-1) = e_count_is
     CALL add_list_is(2*ma)
     !    ******* Off Diagonal Block Terms  ***************
     IF(ci(1,ma) > 0) THEN
        CALL add_list_is(2*ci(1,ma)-1)
        CALL add_list_is(2*ci(1,ma))
     END IF
     IF(ci(2,ma) > 0) THEN
        CALL add_list_is(2*ci(2,ma)-1)
        CALL add_list_is(2*ci(2,ma))
     END IF
     IF(ci(3,ma) > 0) THEN
        CALL add_list_is(2*ci(3,ma)-1)
        CALL add_list_is(2*ci(3,ma))
     END IF
     IF(ci(4,ma) > 0) THEN
        CALL add_list_is(2*ci(4,ma)-1)
        CALL add_list_is(2*ci(4,ma))
     END IF
     IF(ci(5,ma) > 0) THEN
        CALL add_list_is(2*ci(5,ma)-1)
        CALL add_list_is(2*ci(5,ma))
     END IF
     IF(ci(6,ma) > 0) THEN
        CALL add_list_is(2*ci(6,ma)-1)
        CALL add_list_is(2*ci(6,ma))
     END IF
     ! ...  Assembly of active cell; enthalpy
     !    ******* Main Diagonal Block Terms *********
     CALL add_list_is(2*ma)    ! ... diagonal is loaded first
     ia(2*ma) = e_count_is
     CALL add_list_is(2*ma-1)
     !    ******* Off Diagonal Block Terms  ***************
     IF(ci(1,ma) > 0) THEN
        CALL add_list_is(2*ci(1,ma)-1)
        CALL add_list_is(2*ci(1,ma))
     END IF
     IF(ci(2,ma) > 0) THEN
        CALL add_list_is(2*ci(2,ma)-1)
        CALL add_list_is(2*ci(2,ma))
     END IF
     IF(ci(3,ma) > 0) THEN
        CALL add_list_is(2*ci(3,ma)-1)
        CALL add_list_is(2*ci(3,ma))
     END IF
     IF(ci(4,ma) > 0) THEN
        CALL add_list_is(2*ci(4,ma)-1)
        CALL add_list_is(2*ci(4,ma))
     END IF
     IF(ci(5,ma) > 0) THEN
        CALL add_list_is(2*ci(5,ma)-1)
        CALL add_list_is(2*ci(5,ma))
     END IF
     IF(ci(6,ma) > 0) THEN
        CALL add_list_is(2*ci(6,ma)-1)
        CALL add_list_is(2*ci(6,ma))
     END IF
  END DO
  npaja = e_count_is
  ia(2*mrnomax+1) = e_count_is+1
  ! ... Allocate and load the ja array
  ALLOCATE(ja(npaja), jb(npaja),  &
       STAT=a_err)
  IF (a_err /= 0) THEN  
     PRINT *, "Array allocation failed: ldiaja"
     ierr(199) = .TRUE.
     RETURN
  ENDIF
  CALL array_list_is(ja)
END SUBROUTINE ldiaja
