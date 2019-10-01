SUBROUTINE solve(mtrxopp,nbby)
  !     Purpose:  SOLVE is an asymmetric band matrix solver.  It first
  !          triangularizes the band matix (A).
  !          In the second call to SOLVE, the RHS column vector (R) is
  !          modified to correspond with the earlier triangularization of
  !          the RHS (A), and then the matrix is solved by back
  !          substitution.
  !     Originally Programed by James O. Duguid
  USE machine_constants, ONLY: kdp
  USE solver
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: mtrxopp
  INTEGER, INTENT(IN) :: nbby
  !       MTRXOPP (input) - type of matrix operation
  !          =1 triangularize band matrix (ALIT)
  !          =2 solve for right side (RLIT), solution returns in (RLIT)
  !          =3 do both of the above operations
  !       NBBY (input) - number of unknown variables for the given slice J
  !
  INTEGER :: i, ibacj, iband, ihalfb, ihbp, j, jc, ji, jj, jr, k, kc,  &
       ki, kp, lim, mr, nn, nrs
  REAL(kind=kdp)  :: csolv, pivot, solvsum
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2000/12/12 21:27:01 $'
  !     ------------------------------------------------------------------
  nrs = nbby - 1
  ihalfb = (mbw-1)/2
  ihbp = ihalfb + 1
  IF (mtrxopp /= 2) THEN
     !     ******* triangularize matrix A using Doolittle method ********
     DO  j = 1, nrs
        pivot = alit(j,ihbp)
        jj = j + 1
        jc = ihbp
        DO  i = jj, MIN(nbby,jj+ihbp-2)
           jc = jc - 1
           csolv = -alit(i,jc)/pivot
           alit(i,jc) = csolv
           IF (alit(i,jc) == 0._kdp) CYCLE
           ji = jc + 1
           lim = jc + ihalfb
           DO  k = ji, lim
              kc = ihbp + k - jc
              alit(i,k) = alit(i,k) + csolv*alit(j,kc)
           END DO
        END DO
     END DO
     IF (mtrxopp == 1) RETURN
  END IF
  !     ******************** modify right hand side vector R ********************
  nn = nbby + 1
  iband = 2*ihalfb + 1
  DO  i = 2, nbby
     kc = ihbp - i + 1
     ki = 1
     IF (kc < 1) THEN
        kc = 1
        ki = i - ihbp + 1
     END IF
     solvsum = 0._kdp
     DO  k = kc, ihalfb
        solvsum = solvsum + alit(i,k)*rlit(ki)
        ki = ki + 1
     END DO
     rlit(i) = rlit(i) + solvsum
  END DO
  !     *********************** Back Solution ****************************
  rlit(nbby) = rlit(nbby)/alit(nbby,ihbp)
  DO  ibacj = 2,nbby
     i = nn - ibacj
     kp = i
     jr = ihbp + 1
     mr = MIN0(iband, ihalfb+ibacj)
     solvsum = 0._kdp
     DO  k = jr, mr
        kp = kp + 1
        solvsum = solvsum + alit(i,k)*rlit(kp)
     END DO
     rlit(i) = (rlit(i)-solvsum)/alit(i,ihbp)
  END DO
END SUBROUTINE solve
