SUBROUTINE max_ddv_NR(pconvrgmx, pconvrgpc, hconvrgmx, hconvrgpc,  &
     icpdiffmax, ichdiffmax, icpdiffpc, ichdiffpc)
  ! ... Purpose:  To compute the maximum change in P and H for each active
  ! ...        node from the most recent Newton-Raphson iteration. L-inf norm.
  ! ... Also saves the node locations of the maximum changes.
  USE machine_constants, ONLY: kdp
  USE mesh
  USE variables
  IMPLICIT NONE
  REAL(KIND=kdp), intent(out) :: hconvrgmx, hconvrgpc, pconvrgmx, pconvrgpc
  INTEGER, INTENT(OUT) :: icpdiffmax, ichdiffmax, icpdiffpc, ichdiffpc
  ! ...   ICPDIFFMAX - IC index of node with max change in P
  ! ...   ICHDIFFMAX - IC index of node with max change in H
  ! ...   ICPDIFFPC  - IC index of node with max % change in P
  ! ...   ICHDIFFPC  - IC index of node with max % change in H
  !
  INTEGER :: ic, lll
  REAL(KIND=kdp) :: maxdelh, maxdelhpc, maxdelppc, maxdelp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.1 $//$Date: 2003/09/11 16:10:15 $'
  !     ------------------------------------------------------------------
  ! ... Initialize the IC index of the max P & H Diff to the first active node
  icpdiffmax = npic(1)
  ichdiffmax = npic(1)
  icpdiffpc = npic(1)
  ichdiffpc = npic(1)
  pconvrgmx = p(icpdiffmax) - poldnr(icpdiffmax)
  pconvrgpc = (p(icpdiffpc) - poldnr(icpdiffpc))*100._kdp/p(icpdiffpc)
  hconvrgmx = h(ichdiffmax) - holdnr(ichdiffmax)
  hconvrgpc = (h(ichdiffpc) - holdnr(ichdiffpc))*100._kdp/h(ichdiffpc)
  DO  lll = 2, npicmax  ! ... check all the active nodes
     ic = npic(lll)
     ! ... Find the maximum change in P and H and determine the IC indices
     maxdelp = p(ic) - poldnr(ic) 
     IF(ABS(maxdelp) > ABS(pconvrgmx)) THEN
        pconvrgmx = maxdelp
        icpdiffmax = ic
     END IF
     maxdelh = h(ic) - holdnr(ic)
     IF(ABS(maxdelh) > ABS(hconvrgmx)) THEN
        hconvrgmx = maxdelh
        ichdiffmax = ic
     END IF
     ! ... Find the maximum percent change in P and H and determine the IC indices
     maxdelppc = 100._kdp*maxdelp/p(ic)
     IF(ABS(maxdelppc) > ABS(pconvrgpc)) THEN
        pconvrgpc = maxdelppc
        icpdiffpc = ic
     END IF
     maxdelhpc = 100._kdp*maxdelh/h(ic)
     IF(ABS(maxdelhpc) > ABS(hconvrgpc)) THEN
        hconvrgpc = maxdelhpc
        ichdiffpc = ic
     END IF
  END DO
END SUBROUTINE max_ddv_NR
