SUBROUTINE damper
  ! ... Computes the under-relaxation factor for changes in
  ! ...      the dependent variables over a Newton iteration
  ! ... Uses the Cooley algorithm WRR,19,5,p.1274
  ! ... Only for unconfined flow
  USE machine_constants, ONLY: kdp
  USE control
  USE mesh
  USE variables
  IMPLICIT NONE
  INTEGER :: ic, lll
  REAL(kind=kdp) :: s, wstar, maxdp, maxdh, delpmxnr, delhmxnr,  &
       absmaxdp, absmaxdh, maxdpold, maxdhold
  SAVE maxdp, maxdh
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.3 $//$Date: 2004/09/20 21:51:00 $'
  !-----------------------------------------------------------------------
  ! ...  Maximum allowed changes over a NR iteration based on max over time step
  delpmxnr = pchgmxts          ! .. unconfined pressure change ts limit is absolute
  delhmxnr = (hchgmxts/100._kdp)*maxenthic   ! .. enthalpy change is relative (percent)
  !..
  maxdpold = maxdp
  maxdhold = maxdh
  maxdp = dpn(npic(1))
  maxdh = dhn(npic(1))
  absmaxdp = 0._kdp
  absmaxdh = 0._kdp
  DO  lll=1,npicmax
     ic = npic(lll)
     IF(ABS(dpn(ic)) > absmaxdp) THEN
        maxdp = dpn(ic)
        absmaxdp = ABS(maxdp)
     END IF
     IF(ABS(dhn(ic)) > absmaxdh) THEN
        maxdh = dhn(ic)
        absmaxdh = ABS(maxdh)
     END IF
  END DO
  ! ... For pressure
  ! ... Step 1
  s = 1._kdp
  IF(lnri > 1) s = maxdp/(wrelaxp*maxdpold)
  ! ... Step 2
  IF(s >= -1._kdp) THEN
     wstar = (3._kdp+s)/(3._kdp+ABS(s))
  ELSE
     wstar = 1._kdp/(2._kdp*ABS(s))
  END IF
  ! ... Step 3
  IF(wstar*ABS(maxdp) <= delpmxnr) THEN
     wrelaxp = wstar
  ELSE
     wrelaxp = delpmxnr/absmaxdp
  END IF
  ! ... For enthalpy
  ! ... Step 1
  s = 1._kdp
  IF(lnri > 1) s = maxdh/(wrelaxh*maxdhold)
  ! ... Step 2
  IF(s >= -1._kdp) THEN
     wstar = (3._kdp+s)/(3._kdp+ABS(s))
  ELSE
     wstar = 1._kdp/(2._kdp*ABS(s))
  END IF
  ! ... Step 3
  IF(wstar*ABS(maxdh) <= delhmxnr) THEN
     wrelaxh = wstar
  ELSE
     wrelaxh = delhmxnr/absmaxdh
  END IF
END SUBROUTINE damper
