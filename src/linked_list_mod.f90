MODULE linked_list_mod
  ! ... Subroutines for building simple linked list 
  ! ... and writing results to arrays
  ! ...
  ! ... Calling sequence:
  ! ... init_list_*
  ! ... add_list_*
  ! ... array_list_*
  ! ...     where * can be: is - integer scalar; rs - real scalar; ri - both
  ! ... Must use in calling sequence
  ! ...
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  ! ... Change ARRAY_DIM to increase or decrease data block associated with 
  ! ... each pointer in linked list.
  INTEGER, PARAMETER :: array_dim=500
  INTEGER, SAVE :: e_count_rs, e_count_is, e_count_ri  ! ... default counter names
  ! ... 
  ! ... Integer Scalar Linked List
  TYPE :: integer_scalar
     INTEGER, DIMENSION(1:array_dim) :: value
     TYPE(integer_scalar), POINTER :: next
  END TYPE integer_scalar
  ! ... optional head, tail and counter
  TYPE :: list_is
     TYPE(integer_scalar), POINTER :: head, tail
     INTEGER :: counter
     INTEGER, DIMENSION(1:array_dim) :: scratch
  END TYPE list_is
  ! ... default head and tail
  TYPE(integer_scalar), POINTER :: hd_is, tl_is
  INTEGER, DIMENSION(1:array_dim) :: scratch_is
  ! ...
  ! ... Real Scalar Linked List
  TYPE :: real_scalar
     REAL(kind=kdp), DIMENSION(1:array_dim) :: value
     TYPE(real_scalar), POINTER :: next
  END TYPE real_scalar
  ! ... optional head, tail and counter
  TYPE :: list_rs
     TYPE(real_scalar), POINTER :: head, tail
     INTEGER :: counter
     REAL(kind=kdp), DIMENSION(1:array_dim) :: scratch
  END TYPE list_rs
  ! ... default head and tail
  TYPE(real_scalar), POINTER :: hd_rs, tl_rs
  REAL(kind=kdp), DIMENSION(1:array_dim) :: scratch_rs
  ! ...
  ! ... Real/Integer Combination Linked List
  TYPE :: combo
     INTEGER :: index, i_value
     REAL(kind=kdp) :: r_value
  END TYPE combo
  ! ...
  INTEGER :: max_index, min_index
  !.....Set string for use with RCS ident command
  CHARACTER(LEN=80), PRIVATE :: ident_string='$Revision: 1.2 $//$Date: 2003/01/22 23:25:29 $'

CONTAINS

  ! ... LINKED LIST FOR INTEGER SCALARS
  ! ... 
  ! ... Optional STATUS for error detection in linked list.
  ! ... If using ARRAY_LIST_IS to create data arrays, memory for arrays
  ! ... must be allocated in calling program prior to CALL statement.
  ! ... 
  ! ... Optional NAME:
  ! ... Default (no NAME) can only be used with sequential 
  ! ... linked list creation and deallocation-write.
  ! ... For simultaneous creation of two or more linked lists, 
  ! ... second and subsequent list(s) must use NAME option.
  ! ... NAME must be declared as type(LIST_IS) in calling routine.
  ! ... NAME option consists of an alternate name for list 
  ! ... head, tail and counter.
  ! ... 
  ! ... ARRAY_DIM is size of array associated with every pointer in
  ! ... linked list
  ! ... 

  SUBROUTINE init_list_is(name, status)
    ! ... initialize head and tail for integer_scalar list and 
    ! ... general entry counter
    ! ...
    TYPE(list_is), INTENT(INOUT), OPTIONAL :: name
    INTEGER, INTENT(OUT), OPTIONAL :: status
    !     ------------------------------------------------------------------
    !...
    IF (PRESENT(status)) status=0
    IF (PRESENT(name)) THEN
       IF (ASSOCIATED(name%head)) THEN
          ! ... Associated list head indicates that list already exists.
          IF (PRESENT(status)) THEN
             status=2
             RETURN
          ELSE
             PRINT *,'Associated list head indicates that list integer_scalar &
                  &already exists.'
             STOP
          ENDIF
       ENDIF
       NULLIFY(name%head,name%tail)
       name%counter=0
       name%scratch=0
    ELSE           ! default
       IF (ASSOCIATED(hd_is)) THEN
          ! ... Associated list head indicates that list already exists.
          IF (PRESENT(status)) THEN
             status=2
             RETURN
          ELSE
             PRINT *,'Associated list head indicates that list integer_scalar &
                  &already exists.'
             STOP
          ENDIF
       ENDIF
       ! ...
       NULLIFY(hd_is, tl_is)
       e_count_is=0
       scratch_is=0
    ENDIF
  END SUBROUTINE init_list_is

  SUBROUTINE add_list_is(item, name, status)
    ! ... Construct linked list
    ! ...
    INTEGER, INTENT(in) :: item
    TYPE(list_is), INTENT(inout), OPTIONAL :: name
    INTEGER, INTENT(out), OPTIONAL :: status
    !
    INTEGER :: count
    !     ------------------------------------------------------------------
    !...
    IF (PRESENT(status)) status=0
    ! ... build scratch array
    IF (PRESENT(name)) THEN
       count=MOD(name%counter,array_dim)+1
       name%counter=name%counter+1        ! count entries
       name%scratch(count)=item
    ELSE
       count=MOD(e_count_is,array_dim)+1
       e_count_is=e_count_is+1          ! count entries
       scratch_is(count)=item
    ENDIF
    ! ... When scratch full, place in linked list
    IF (count==array_dim) CALL flush_to_list_is(name=name, status=status)
  END SUBROUTINE add_list_is

  SUBROUTINE flush_to_list_is(name, status)
    ! ... Create linked list of scratch arrays
    TYPE(list_is), INTENT(inout), OPTIONAL :: name
    INTEGER, INTENT(inout), OPTIONAL :: status
    ! ...
    INTEGER :: err
    TYPE(integer_scalar), POINTER :: new
    !     ------------------------------------------------------------------
    !...
    ALLOCATE(new, stat=err)              ! allocate new storage location
    IF (err /= 0) THEN
       ! ... error in allocation of space
       IF (PRESENT(status)) THEN
          status=1
          RETURN
       ELSE
          PRINT*,'allocation error in flush_to_list_is: stat=',err
          STOP
       ENDIF
    ENDIF
    NULLIFY(new%next)            ! no successor to new item: not associated
    IF (PRESENT(name)) THEN
       new%value=name%scratch    ! Set scratch space to value in linked list
       name%scratch=0
       ! ...
       IF (ASSOCIATED(name%head)) THEN
          ! ... primary list is not empty
          name%tail%next => new  ! point next to new storage location 
          name%tail => new       ! make tail equivalent to new storage location
       ELSE
          ! ... primary list is empty: add first item
          name%head => new
          name%tail => new
       ENDIF
    ELSE ! default
       new%value=scratch_is      ! Set scratch space to value in linked list
       scratch_is=0
       ! ...
       IF (ASSOCIATED(hd_is)) THEN
          ! ... list is not empty
          tl_is%next => new    ! point next to new storage location 
          tl_is => new         ! make tail equivalent to new storage location
       ELSE
          ! ... list is empty: add first item
          hd_is => new
          tl_is => new
       ENDIF
    ENDIF
  END SUBROUTINE flush_to_list_is

  SUBROUTINE array_list_is(array, name, status)
    ! ... Write linked list to an array
    INTEGER, INTENT(OUT), DIMENSION(:) :: array
    TYPE(list_is), INTENT(INOUT), OPTIONAL :: name
    INTEGER, INTENT(OUT), OPTIONAL :: status
    ! ...
    TYPE(integer_scalar), POINTER :: point
    INTEGER :: i, count
    !     ------------------------------------------------------------------
    !...
    IF (PRESENT(status)) status=0
    IF (PRESENT(name)) THEN 
       ! ... flush last array into linked list
       IF (MOD(name%counter,array_dim)/=0) CALL flush_to_list_is(name=name, status=status)
       IF (ASSOCIATED(name%head)) THEN
          point => name%head
          DO i=1, name%counter
             ! ...
             ! ... Take information from storage site, 
             ! ... then deallocate storage site.
             count=MOD(i-1,array_dim)+1
             array(i) = point%value(count)
             IF (count<array_dim) CYCLE
             ! ... Go to next pointer in list
             point => point%next
             CALL deallocate_is('array_list_is', name=name, status=status)
             name%head=>point
          ENDDO
          IF (ASSOCIATED(name%head)) CALL &
               deallocate_is('array_list_is', name=name, status=status)
       ELSE
          ! ... probable programming error; linked list array_list_is empty
          IF (PRESENT(status)) THEN
             status=-2
             RETURN
          ELSE
             PRINT*,'probable programming error; linked list integer_scalar &
                  &empty'
             STOP
          ENDIF
       ENDIF
    ELSE           ! default: no name
       ! ... flush last array into linked list
       IF (MOD(e_count_is,array_dim)/=0) CALL flush_to_list_is(name=name, status=status)
       IF (ASSOCIATED(hd_is)) THEN
          point => hd_is
          DO i=1, e_count_is
             ! ... Take information from storage site, 
             ! ... then deallocate storage site.
             count = MOD(i-1,array_dim)+1
             array(i) = point%value(count)
             IF (count<array_dim) CYCLE
             ! ... Go to next pointer in list
             point => point%next
             CALL deallocate_is('array_list_is', name=name, status=status)
             hd_is=>point
          ENDDO
          IF (ASSOCIATED(hd_is)) CALL &
               deallocate_is('array_list_is', name=name, status=status)
       ELSE
          ! ... probable programming error; linked list array_list_is empty
          IF (PRESENT(status)) THEN
             status=-2
             RETURN
          ELSE
             PRINT*,'probable programming error; linked list integer_scalar &
                  &empty'
             STOP
          ENDIF
       ENDIF
    ENDIF              ! end, if (present(name) . . .
  END SUBROUTINE array_list_is

  SUBROUTINE deallocate_is(sub_name,name,status)
    ! ... deallocate integer_scalar linked list head
    TYPE(list_is), INTENT(INOUT), OPTIONAL :: name
    CHARACTER(LEN=*), INTENT(IN) ::sub_name
    INTEGER, INTENT(OUT), OPTIONAL :: status
    INTEGER :: err
    !     ------------------------------------------------------------------
    !...
    IF (PRESENT(name)) THEN
       DEALLOCATE(name%head, stat=err)
    ELSE
       DEALLOCATE(hd_is, stat=err)
    ENDIF
    IF (err /= 0) THEN
       ! ... error in deallocation of space
       IF (PRESENT(status)) THEN
          status=-1
          RETURN
       ELSE
          PRINT*,'deallocation error ',TRIM(sub_name),': stat=',err
          STOP
       ENDIF
    ENDIF
  END SUBROUTINE deallocate_is

END MODULE linked_list_mod
