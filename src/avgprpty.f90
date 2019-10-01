SUBROUTINE avgprpty
  ! ... Purpose:  To compute the average density and viscosity of
  ! ...      steam and water at the cell boundary between two nodes.
  ! ...      Also to determine the derivatives of those properties
  ! ... Method:  Compute the average P & H for the cell boundary between two nodes.
  ! ...      Then determine the density and viscosity.  If the boundary
  ! ...      lies in the two phase region, determine the
  ! ...      density and viscosity for the coexisting steam and water for
  ! ...      that average P & H. When the two cells being averaged lie on
  ! ...      different sides of the two-phase region or in the supercritical
  ! ...      region, then special techniques are applied to determine the
  ! ...      average properties.  These techniques attempt to make smooth
  ! ...      transitions, so that there are no marked jumps in density or
  ! ...      viscosity, when a different technique is applied.
  USE machine_constants, ONLY: kdp
  USE fdeq
  USE mesh
  USE parameters
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, icxf, icxl, icxm, icyf, icyl, icym, iczf, iczl,  &
       iczm, ik, ikxm, ikzm, j, k, mxf, mxl, mxm, myf
  INTEGER :: myl, mym, mzf, mzl, mzm
  REAL(KIND=kdp) :: dum1, dum2, dum3, dum4, dum5, dum6
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.4 $//$Date: 2006/09/06 21:04:28 $'
  !     ------------------------------------------------------------------
  !...
  ! ...   Loop to determine the weighted arithmetic averages of
  ! ...        density and viscosity terms at the cell face and an arithmetic weighting
  ! ...        factor to apply to the derivatives
  ! ... Compute arithmetic mean of density/viscosity terms
  ! ...        --------  X-direction ----------
  DO  j = 1, ny
     DO  k = 1, nz
        DO  i = 2, nx
           ik = (k-1)*nx + i
           ikxm = ik - 1
           ! ...           if the node or its adjacent node (i-1) is not in domain
           ! ...              then skip to end of loop
           IF (np(ik,j) == 0 .OR. np(ikxm,j) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nx*nz
           icxm = ic - 1
           mxm = (k-1)*nxx + i + (j-1)*nxx*nz
!!$       write (fustdout,*) 'in avgprpty IC, ICXM=', ic, icxm
           CALL avgvalue(ic, icxm, xdw(mxm), xvw(mxm),  &
                xvdw(mxm), dxvdwp(mxm), dxvdwh(mxm), dum1, dum2, dum3,  &
                xds(mxm), xvs(mxm), xvds(mxm), dxvdsp(mxm), dxvdsh(mxm),  &
                dum4, dum5, dum6)
        END DO
     END DO
  END DO
  ! ...         ----------- Z-direction ---------------
  DO  j = 1, ny
     DO  i = 1, nx
        DO  k = 2, nz
           ik = (k-1)*nx + i
           ikzm = ik - nx
           ! ...         if the node or its adjacent node (Z-1) is not in domain
           ! ...              then skip to end of loop
           IF (np(ik,j) == 0 .OR. np(ikzm,j) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nx*nz
           iczm = ic - nx
           mzm = (k-1)*nx + i + (j-1)*nx*nzz
!!$       write (fustdout,*) 'in avgprpty IC, ICZM=', ic, iczm
           CALL avgvalue(ic, iczm, zdw(mzm), zvw(mzm),  &
                zvdw(mzm), dzvdwp(mzm), dzvdwh(mzm),  &
                zvdgw(mzm), dzvdgwp(mzm), dzvdgwh(mzm), zds(mzm), zvs(mzm),  &
                zvds(mzm), dzvdsp(mzm), dzvdsh(mzm),  &
                zvdgs(mzm), dzvdgsp(mzm), dzvdgsh(mzm))
        END DO
     END DO
  END DO
  ! ...        ------------ Y-direction --------------
  ! ...     if only one block deep, skip calculations
  DO  i = 1, nx
     DO  k = 1, nz
        DO  j = 2, ny
           ik = (k-1)*nx + i
           ! ...         if the node or its adjacent node (Y-1) is not in domain
           ! ...              then skip to end of loop
           IF (np(ik,j) == 0 .OR. np(ik,j-1) == 0) CYCLE
           ic = (k-1)*nx + i + (j-1)*nx*nz
           icym = ic - nx*nz
           mym = ic
           CALL avgvalue(ic, icym, ydw(mym), yvw(mym),  &
                yvdw(mym), dyvdwp(mym), dyvdwh(mym), dum1, dum2, dum3,  &
                yds(mym), yvs(mym), yvds(mym), dyvdsp(mym), dyvdsh(mym),  &
                dum4, dum5, dum6)
        END DO
     END DO
  END DO
  ! ...   set density/viscosity of boundary cells to value of nearest node
  ! ...     X boundary cells, ICXF & MXF: left,  ICXL & MXL: right
  DO  j = 1, ny
     DO  k = 1, nz
        icxf = (k-1)*nx + 1 + (j-1)*nx*nz
        icxl = k*nx + (j-1)*nx*nz
        mxf = (k-1)*nxx + 1 + (j-1)*nxx*nz
        mxl = k*nxx + (j-1)*nxx*nz
        ! ...               water
        xvw(mxf) = visw(icxf)
        xvw(mxl) = visw(icxl)
        xdw(mxf) = denw(icxf)
        xdw(mxl) = denw(icxl)
        xvdw(mxf) = denw(icxf)/visw(icxf)
        xvdw(mxl) = denw(icxl)/visw(icxl)
        ! ...                steam
        xvs(mxf) = viss(icxf)
        xvs(mxl) = viss(icxl)
        xds(mxf) = dens(icxf)
        xds(mxl) = dens(icxl)
        xvds(mxf) = dens(icxf)/viss(icxf)
        xvds(mxl) = dens(icxl)/viss(icxl)
     END DO
  END DO
  ! ...     Z boundary cells, ICZF & MZF: bottom,  ICZL & MZL: top
  DO  j = 1, ny
     DO  i = 1, nx
        iczf = i + (j-1)*nz*nx
        iczl = (nz-1)*nx + i + (j-1)*nz*nx
        mzf = i + (j-1)*nzz*nx
        mzl = nz*nx + i + (j-1)*nzz*nx
        ! ...                water
        zvw(mzf) = visw(iczf)
        zvw(mzl) = visw(iczl)
        zdw(mzf) = denw(iczf)
        zdw(mzl) = denw(iczl)
        zvdw(mzf) = denw(iczf)/visw(iczf)
        zvdw(mzl) = denw(iczl)/visw(iczl)
        zvdgw(mzf) = grav*denw(iczf)**2/visw(iczf)
        zvdgw(mzl) = grav*denw(iczl)**2/visw(iczl)
        ! ...                steam
        zvs(mzf) = viss(iczf)
        zvs(mzl) = viss(iczl)
        zds(mzf) = dens(iczf)
        zds(mzl) = dens(iczl)
        zvds(mzf) = dens(iczf)/viss(iczf)
        zvds(mzl) = dens(iczl)/viss(iczl)
        zvdgs(mzf) = grav*dens(iczf)**2/viss(iczf)
        zvdgs(mzl) = grav*dens(iczl)**2/viss(iczl)
     END DO
  END DO
  ! ...     Y boundary cells, ICYF & MYF: front,  ICYL & MYL: back
  IF (ny > 1) THEN
     DO  k = 1, nz
        DO  i = 1, nx
           icyf = (k-1)*nx + i
           icyl = (k-1)*nx + i + (ny-1)*nx*nz
           myf = (k-1)*nx + i
           myl = myf + nx*nz*ny
           ! ...               water
           yvw(myf) = visw(icyf)
           yvw(myl) = visw(icyl)
           ydw(myf) = denw(icyf)
           ydw(myl) = denw(icyl)
           yvdw(myf) = denw(icyf)/visw(icyf)
           yvdw(myl) = denw(icyl)/visw(icyl)
           ! ...               steam
           yvs(myf) = viss(icyf)
           yvs(myl) = viss(icyl)
           yds(myf) = dens(icyf)
           yds(myl) = dens(icyl)
           yvds(myf) = dens(icyf)/viss(icyf)
           yvds(myl) = dens(icyl)/viss(icyl)
        END DO
     END DO
  END IF
END SUBROUTINE avgprpty
