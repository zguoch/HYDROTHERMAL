SUBROUTINE rowscale(nrow,job,norm,a,ja,ia,diag,b,jb,ib,ierr)
  ! ... Scales the rows of matrix A such that their norms are one
  ! ... Choices of norms: 1-norm, 2-norm, max-norm (infinity-norm).
  ! ... B = Diag*A
  ! Notes:
  !     algorithm could be done in place (B can take the place of A).
  !-----------------------------------------------------------------
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nrow    ! ... The row dimension of A
  INTEGER, INTENT(IN) :: job     ! ... job selector. 0: calculate B only
                                 ! ...   1: calculate B and load ib,jb
  INTEGER, INTENT(IN) :: norm     ! ... norm selector. 1: 1-norm,
                                 ! ... 2: 2-norm, 0: max-norm
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
  INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
  INTEGER, DIMENSION(:), INTENT(IN) :: ia 
  REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag     ! ... Diagonal matrix of scale factors
  REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b        ! ... Matrix B in compressed
  INTEGER, DIMENSION(:), INTENT(OUT) :: jb              ! ... sparse row format
  INTEGER, DIMENSION(:), INTENT(OUT) :: ib
  INTEGER, INTENT(OUT) :: ierr              ! ... error message. 0 : Normal return
                                            ! ...   i > 0 : Row number i is a zero row
  !
  INTEGER :: j
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  !...
  CALL rownorms(norm,a,ja,ia,diag)
  ierr = 0
  DO  j=1,nrow
     IF (diag(j) == 0.0_kdp) THEN
        ierr = j
        RETURN
     ELSE
        diag(j) = 1.0_kdp/diag(j)
     END IF
  END DO
  ! ... Scale matrix A
  CALL diamua(job,a,ja,ia,diag,b,jb,ib)

  CONTAINS

    SUBROUTINE rownorms(norm,a,ja,ia,diag)
      ! ... Calculates the norms of each row of A. (choice of three norms)
      !-----------------------------------------------------------------
!!$      USE machine_constants, ONLY: kdp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: norm     ! ... norm selector. 1:1-norm,
                                     ! ... 2:2-norm, 0: max-norm
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
      INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(IN) :: ia
      REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag     ! ... the norms
      !
      REAL(KIND=kdp) :: scal
      INTEGER :: ii, k, k1, k2
      ! ... Set string for use with RCS ident command
      !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
      !     ------------------------------------------------------------------
      !...
      ! ... nrow is known from host
      DO  ii=1,nrow
         ! ... Compute the norm of each row
         scal = 0.0_kdp
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         IF (norm == 0) THEN
            DO  k=k1, k2
               scal = MAX(scal,ABS(a(k)))
            END DO
         ELSE IF (norm == 1) THEN
            DO  k=k1, k2
               scal = scal + ABS(a(k))
            END DO
         ELSE
            DO  k=k1, k2
               scal = scal+a(k)**2
            END DO
         END IF
         IF (norm == 2) scal = SQRT(scal)
         diag(ii) = scal
      END DO
    END SUBROUTINE rownorms
    
    SUBROUTINE diamua(job,a,ja,ia,diag,b,jb,ib)
      ! ... Calculates a diagonal matrix times a matrix; B = Diag*A  (in place)
      ! Notes:
      !       algorithm in place (B can take the place of A). job=0
      !-----------------------------------------------------------------
!!$      USE machine_constants, ONLY: kdp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: job     ! ... job selector. 0: calculate B only
                                     ! ...   1: calculate B and load ib,jb
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
      INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(IN) :: ia
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: diag  ! ... diagonal matrix
      REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b    ! ... Matrix B in compressed
      INTEGER, DIMENSION(:), INTENT(OUT) :: jb          ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(OUT) :: ib
      !
      REAL(KIND=kdp) :: scal
      INTEGER :: ii, k, k1, k2
      ! ... Set string for use with RCS ident command
      !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
      !     ------------------------------------------------------------------
      !...
      ! ... nrow is known from host
      DO  ii=1,nrow
         !     normalize each row
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii)
         DO  k=k1, k2
            b(k) = a(k)*scal
         END DO
      END DO
      IF (job == 1) THEN
         DO  ii=1, nrow+1
            ib(ii) = ia(ii)
         END DO
         DO  k=ia(1),ia(nrow+1)-1
            jb(k) = ja(k)
         END DO
      END IF
    END SUBROUTINE diamua

END SUBROUTINE rowscale

SUBROUTINE colscale(nrow,job,norm,a,ja,ia,diag,b,jb,ib,ierr)
  ! ... Scales the columns of matrix A such that their norms are one
  ! ... Result matrix written on B, or overwritten on A.
  ! 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
  ! ... B = A*Diag
  !      algorithm could be done in place, B overwritten onto A
  !-----------------------------------------------------------------
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nrow    ! ... The row dimension of A
  INTEGER, INTENT(IN) :: job     ! ... job selector. 0: calculate B only
                                 ! ...   1: calculate B and load ib,jb
  INTEGER, INTENT(IN) :: norm     ! ... norm selector. 1: 1-norm,
                                 ! ... 2: 2-norm, 0: max-norm
  REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
  INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
  INTEGER, DIMENSION(:), INTENT(IN) :: ia
  REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag     ! ... Diagonal matrix of scale factors
  REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b        ! ... Matrix B in compressed
  INTEGER, DIMENSION(:), INTENT(OUT) :: jb              ! ... sparse row format
  INTEGER, DIMENSION(:), INTENT(OUT) :: ib
  INTEGER, INTENT(OUT) :: ierr              ! ... error message. 0 : Normal return
                                            ! ...   i > 0: Column number i is a zero column
  !
  INTEGER :: j
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
  !     ------------------------------------------------------------------
  !...
  CALL colnorms(norm,a,ja,ia,diag)
  ierr = 0
  DO  j=1,nrow
     IF (diag(j) == 0.0_kdp) THEN
        ierr = j
        RETURN
     ELSE
        diag(j) = 1.0_kdp/diag(j)
     END IF
  END DO
  ! ... Scale matrix A
  CALL amudia(job,a,ja,ia,diag,b,jb,ib)

  CONTAINS

    SUBROUTINE colnorms(norm,a,ja,ia,diag)
      ! ... Calculates the norms of each column of A. (choice of three norms)
      !-----------------------------------------------------------------
!!$      USE machine_constants, ONLY: kdp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: norm     ! ... norm selector. 1: 1-norm,
                                     ! ... 2: 2-norm, 0: max-norm
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
      INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(IN) :: ia
      REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: diag     ! ... the norms
      !
      INTEGER :: ii, j, k, k1, k2
      ! ... Set string for use with RCS ident command
      !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
            !     ------------------------------------------------------------------
      !...
      ! ... nrow is known from host
      diag = 0.0_kdp
      DO  ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         DO  k=k1, k2
            j = ja(k)
            ! ... update the norm of each column
            IF (norm == 0) THEN
               diag(j) = MAX(diag(j),ABS(a(k)))
            ELSE IF (norm == 1) THEN
               diag(j) = diag(j) + ABS(a(k))
            ELSE     ! ... 2-norm
               diag(j) = diag(j)+a(k)**2
            END IF
         END DO
      END DO
      IF (norm == 2) THEN
         DO  k=1,nrow
            diag(k) = SQRT(diag(k))
         END DO
      END IF
    END SUBROUTINE colnorms

    SUBROUTINE amudia(job,a,ja,ia,diag,b,jb,ib)
      ! ... Calculates a matrix times a diagonal matrix; B = A*Diag  (in place)
      ! Notes:
      !        algorithm could be done in place (B can take the place of A). job=0
      !-----------------------------------------------------------------
!!$      USE machine_constants, ONLY: kdp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: job     ! ... job selector. 0: calculate B only
                                     ! ...   1: calculate B and load ib,jb
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: a     ! ... Matrix A in compressed
      INTEGER, DIMENSION(:), INTENT(IN) :: ja           ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(IN) :: ia 
      REAL(KIND=kdp), DIMENSION(:), INTENT(IN) :: diag  ! ... diagonal matrix
      REAL(KIND=kdp), DIMENSION(:), INTENT(OUT) :: b    ! ... Matrix B in compressed
      INTEGER, DIMENSION(:), INTENT(OUT) :: jb          ! ... sparse row format
      INTEGER, DIMENSION(:), INTENT(OUT) :: ib
      !
      INTEGER :: ii, k, k1, k2
      ! ... Set string for use with RCS ident command
      !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3.1.3 $//$Date: 2007/10/16 21:17:22 $'
      !     ------------------------------------------------------------------
      !...
      ! ... nrow is known from host
      DO  ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         DO  k=k1, k2
            b(k) = a(k)*diag(ja(k))
         END DO
      END DO
      IF (job == 1) THEN
         DO  ii=1, nrow+1
            ib(ii) = ia(ii)
         END DO
         DO  k = ia(1),ia(nrow+1)-1
            jb(k) = ja(k)
         END DO
      END IF
    END SUBROUTINE amudia

END SUBROUTINE colscale
