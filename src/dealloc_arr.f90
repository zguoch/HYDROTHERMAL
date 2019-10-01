SUBROUTINE dealloc_arr
  ! ... Deallocates the array space used for the simulation
  USE bc
  USE parameters          !!$*** for flint
  USE control
  USE i_c
  USE fdeq
  USE mesh
!!$  USE parameters
  USE solver
  USE solver_gmres
  USE source
  USE units
  USE variables
  IMPLICIT NONE
  !
  INTEGER :: da_err
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.13 $//$Date: 2007/12/19 20:55:47 $'
  !     ------------------------------------------------------------------
  IF(ALLOCATED(ax)) DEALLOCATE(ax, ay, az, dx, dy, dz, ci, mrno,  &
       nb, npic, np, npp, lprntxif, lprntyif, lprntzif,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #1"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(ptopa)) DEALLOCATE(ptopa, cdtn, qprecip, tflux, denflux, ehflux, qhflux,  &
       land_surf_area, qseep, ibc, prflag, seeping, ehassoc, tcassoc,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #2"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(c)) DEALLOCATE(c, d, f, g,  &
       dhrsh, dhrsp, dhrwh, dhrwp, dqhh,  &
       dqhp, dth, dtp, drsh, drsp, drwh, drwp,  &
       dxvdsh, dxvdsp, dxvdwh, dxvdwp,  &
       dyvdsh, dyvdsp, dyvdwh, dyvdwp,  &
       dzvdsh, dzvdsp, dzvdwh, dzvdwp,  &
       dzvdgsh, dzvdgsp, dzvdgwh, dzvdgwp,  &
       reside, residm,  &
       hrs, hrw,  &
       tx, txk, ty, tyk, tz, tzk,  &
       usx, usy, usz, uwx, uwy, uwz,  &
       xvds, xvdw, xvs, xvw, xds, xdw,  &
       yvds, yvdw, yvs,  & 
       yvw, yds, ydw,  &
       zvdgs, zvdgw, zvds, zvdw, zvs, zvw,  &
       zds, zdw,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #3"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(pinit)) DEALLOCATE(pinit, plith, phiinit,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #4"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(rs)) DEALLOCATE(rs, rw, phi, xkc, ykc, zkc,  &
       xk, yk, zk,  &
       beta, df, phfwt, phfwtdt,  &
       icrxtype,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #5"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(q)) DEALLOCATE(q, qh, qhi, qhcondflx,  &
       qhinodsrc, qhnodsrc, qhwelsrc, qnodsrc, qwelsrc, qt, qv, iq,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #6"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(ALLOCATED(ind)) DEALLOCATE(ind, indoldnr, indoldt,  &
       dens, denw, viss, visw, satnw, wpot, spot,  &
       swoldt, hspot, hwpot, fs_elev,  &
       tc,  &
       dhn, dpn,  &
       en, enoldt, enoldt2, xm, xmoldt, xmoldt2,  &
       h, holdt, holdnr, p, pe, poldt, poldnr, nu,  &
       xsvel, xwvel, ysvel, ywvel, zsvel, zwvel,  &
       xsmflx, xwmflx, ysmflx, ywmflx, zsmflx, zwmflx,  &
       hrock,  &
       xsflux, xwflux, ysflux, ywflux, zsflux,  &
       zwflux,  &
       rsq, xx, yy, zz, zzz, zztop, zls,  &
       ktop,  &
       STAT=da_err)
  IF (da_err /= 0) THEN  
     PRINT *, "Array deallocation failed: dealloc_arr, #7"
     ierr(200) = .TRUE.
     RETURN
  ENDIF
  IF(slmeth == 1) THEN
     DEALLOCATE(alit, rlit, cscal,  &
          STAT=da_err)
     IF (da_err /= 0) THEN  
        PRINT *, "Array deallocation failed: dealloc_arr, #8.1"
        ierr(200) = .TRUE.
        RETURN
     ENDIF
  ELSEIF(slmeth == 2) THEN
     DEALLOCATE(av, ia, ja, rhs, xxs, bv, ib, jb, diagc, diagr,  &
          ailu, jailu, iilu, levs,  &
          STAT=da_err)
     IF (da_err /= 0) THEN
        PRINT *, "Array deallocation failed: dealloc_arr, #8.2"
        ierr(200) = .TRUE.
        RETURN
     ENDIF
  END IF
END SUBROUTINE dealloc_arr
