SUBROUTINE velocity
  !     Purpose:  To compute the X,Y, and Z velocities of water and steam
  !        for each node.  Velocities are used for plotting.
  !     Units:  Flux variables (e.g. XWFLUX) are volumetric (Darcy) fluxes in
  !        cm^3/cm^2-sec or cm/sec for the cell face.  Interstitial velocity
  !        variables (e.g. ZSVEL) are in cm/sec and are calculated at the
  !        node points. Mass fluxes are calculated at the nodes.
  USE machine_constants, ONLY: kdp
  USE fdeq
  USE mesh
  USE parameters
  USE math_constants
  USE variables
  IMPLICIT NONE
  INTEGER :: i, ic, icxp, icyp, iczp, ik, ikxp, ikzp, j, k, mxf, mxl,  &
       mxm, mxp, myf, myl, mym, myp, mzf, mzl
  INTEGER :: mzm, mzp
  REAL(KIND=kdp) :: rface
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2007/04/10 23:14:34 $'
  !     ------------------------------------------------------------------
  ! ... Calculate specific discharge (Darcy flux) at cell boundaries
  !                  ------------ X-direction ------------
  IF (nx > 1) THEN
     DO  j = 1, ny
        DO  k = 1, nz
           DO  i = 1, nx - 1
              ! ... indices in x-z slice
              ik = (k-1)*nx + i
              ikxp = ik + 1
              ! ... node indices
              ic = (k-1)*nx + i + (j-1)*nx*nz
              icxp = ic + 1
              ! ... indices of cell boundaries X direction
              mxm = (k-1)*nxx + i + (j-1)*nxx*nz
              mxp = mxm + 1
              ! ... set flux to 0 on both sides of cells not in the domain
              IF (np(ik,j) == 0) THEN
                 xwflux(mxm) = 0._kdp
                 xsflux(mxm) = 0._kdp
                 xwflux(mxp) = 0._kdp
                 xsflux(mxp) = 0._kdp
                 ! ... set flux to 0 if next cell (i+1) is not in the domain
              ELSE IF (np(ikxp,j) == 0) THEN
                 xwflux(mxp) = 0._kdp
                 xsflux(mxp) = 0._kdp
              ELSE
                 IF(irad) THEN
                    rface = sqrt(rsq(i+1))
                    ! ... water Darcy flux at cell boundary
                    xwflux(mxp) = -tx(mxp)/(dz(k)*twopi*rface*xvw(mxp))*  &
                         (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))
                    ! ... steam Darcy flux at cell boundary
                    xsflux(mxp) = -tx(mxp)/(dz(k)*twopi*rface*xvs(mxp))*  &
                         (rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*(p(icxp)-p(ic))
                 ELSE     ! ... Cartesian
                    ! ... water Darcy flux at cell boundary
                    xwflux(mxp) = -tx(mxp)/(dz(k)*dy(j)*xvw(mxp))*  &
                         (rw(icxp)*uwx(mxp)+rw(ic)*(1._kdp-uwx(mxp)))*(p(icxp)-p(ic))
                    ! ... steam Darcy flux at cell boundary
                    xsflux(mxp) = -tx(mxp)/(dz(k)*dy(j)*xvs(mxp))*  &
                         (rs(icxp)*usx(mxp)+rs(ic)*(1._kdp-usx(mxp)))*(p(icxp)-p(ic))
                 END IF
              END IF
           END DO
           ! ...        set fluxes at edges of active region to 0
           ! ...        X boundary cell boundaries, MXF - leftmost, MXL - rightmost
           mxf = (k-1)*nxx + 1 + (j-1)*nxx*nz
           mxl = k*nxx + (j-1)*nxx*nz
           xwflux(mxf) = 0._kdp
           xsflux(mxf) = 0._kdp
           xwflux(mxl) = 0._kdp
           xsflux(mxl) = 0._kdp
           ! ...         if last cell in row is not in active region, set interface
           ! ...             flux on left side to 0.
           i = nx
           ik = (k-1)*nx + i
           IF(np(ik,j) == 0) THEN
              mxm = (k-1)*nxx + i + (j-1)*nxx*nz
              xwflux(mxm) = 0._kdp
              xsflux(mxm) = 0._kdp
           END IF
        END DO
     END DO
  END IF
  ! ...              ------------ Y-direction ------------
  IF (ny > 1) THEN
     DO  k = 1, nz
        DO  i = 1, nx
           DO  j = 1, ny - 1
              ! ... indices in x-z slice
              ik = (k-1)*nx + i
              ! ... node indices
              ic = (k-1)*nx + i + (j-1)*nx*nz
              icyp = ic + nx*nz
              ! ... indices of cell boundaries Y direction
              mym = (k-1)*nx + i + (j-1)*nx*nz
              myp = mym + nx*nz
              ! ...        set flux to 0 on both sides of cells not in the active region
              IF (np(ik,j) == 0) THEN
                 ywflux(mym) = 0._kdp
                 ysflux(mym) = 0._kdp
                 ywflux(myp) = 0._kdp
                 ysflux(myp) = 0._kdp
                 ! ...         set flux to 0 if next cell (Y+1) is not in the active region
              ELSE IF (np(ik,j+1) == 0) THEN
                 ywflux(myp) = 0._kdp
                 ysflux(myp) = 0._kdp
              ELSE
                 ! ...                water flux at cell boundary
                 ywflux(myp) = -ty(myp)/(dz(k)*dx(i)*yvw(myp))*  &
                      (rw(icyp)*uwy(myp)+rw(ic)*(1._kdp-uwy(myp)))*(p(icyp)-p(ic))
                 ! ...                steam flux at cell boundary
                 ysflux(myp) = -ty(myp)/(dz(k)*dx(i)*yvs(myp))*  &
                      (rs(icyp)*usy(myp)+rs(ic)*(1._kdp-usy(myp)))*(p(icyp)-p(ic))
              END IF
           END DO
           ! ...        set fluxes at edges of active region to 0
           ! ...        Y boundary cell boundaries, MYF - front, MYL - back
           myf = (k-1)*nx + i
           myl = myf + nx*nz*ny
           ywflux(myf) = 0._kdp
           ysflux(myf) = 0._kdp
           ywflux(myl) = 0._kdp
           ysflux(myl) = 0._kdp
           ! ...         if back node in row is not in active region then set interface
           ! ...             flux on front side to 0.
           j = ny
           ik = (k-1)*nx + i
           IF (np(ik,j) == 0) THEN
              mym = (k-1)*nx + i + (j-1)*nx*nz
              ywflux(mym) = 0._kdp
              ysflux(mym) = 0._kdp
           END IF
        END DO
     END DO
  END IF
  ! ...              ------------ Z-direction ------------
  IF (nz > 1) THEN
     DO  j = 1, ny
        DO  i = 1, nx
           DO  k = 1, nz - 1
              ! ... indices in x-z slice
              ik = (k-1)*nx + i
              ikzp = ik + nx
              ! ... node indices
              ic = (k-1)*nx + i + (j-1)*nx*nz
              iczp = ic + nx
              ! ... indices of cell boundaries Z direction
              mzm = (k-1)*nx + i + (j-1)*nx*nzz
              mzp = mzm + nx
              ! ...        set flux to 0 on both sides of cells not in the active region
              IF (np(ik,j) == 0) THEN
                 zwflux(mzm) = 0._kdp
                 zsflux(mzm) = 0._kdp
                 zwflux(mzp) = 0._kdp
                 zsflux(mzp) = 0._kdp
                 ! ...         set flux to 0 if node above (z+1) is not in the active region
              ELSE IF (np(ikzp,j) == 0) THEN
                 zwflux(mzp) = 0._kdp
                 zsflux(mzp) = 0._kdp
              ELSE
                 IF(irad) THEN
                          ! ... water Darcy flux at cell boundary
                    zwflux(mzp) = -tz(mzp)/(pi*(rsq(i+1)-rsq(i))*zvw(mzp))*  &
                         (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*  &
                         ((p(iczp)-p(ic))+zdw(mzp)*grav*zzz(k))
                    ! ... steam Darcy flux at cell boundary
                    zsflux(mzp) = -tz(mzp)/(pi*(rsq(i+1)-rsq(i))*zvs(mzp))*  &
                         (rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*  &
                         ((p(iczp)-p(ic))+zds(mzp)*grav*zzz(k))
                 ELSE
                    ! ... water Darcy flux at cell boundary
                    zwflux(mzp) = -tz(mzp)/(dx(i)*dy(j)*zvw(mzp))*  &
                         (rw(iczp)*uwz(mzp)+rw(ic)*(1._kdp-uwz(mzp)))*  &
                         ((p(iczp)-p(ic))+zdw(mzp)*grav*zzz(k))
                    ! ... steam Darcy flux at cell boundary
                    zsflux(mzp) = -tz(mzp)/(dx(i)*dy(j)*zvs(mzp))*  &
                         (rs(iczp)*usz(mzp)+rs(ic)*(1._kdp-usz(mzp)))*  &
                         ((p(iczp)-p(ic))+zds(mzp)*grav*zzz(k))
                 END IF
              END IF
           END DO
           ! ...        set fluxes at edges of active region to 0
           ! ...          Z boundary cell boundaries, MZF - bottom,  MZL - top
           mzf = i + (j-1)*nzz*nx
           mzl = nz*nx + i + (j-1)*nzz*nx
           zwflux(mzf) = 0._kdp
           zsflux(mzf) = 0._kdp
           zwflux(mzl) = 0._kdp
           zsflux(mzl) = 0._kdp
           ! ...         if top cell in column is not in active region then set interface
           ! ...             flux on bottom side to 0.
           k = nz
           ik = (k-1)*nx + i
           IF (np(ik,j) == 0) THEN
              mzm = (k-1)*nx + i + (j-1)*nx*nzz
              zwflux(mzm) = 0._kdp
              zsflux(mzm) = 0._kdp
           END IF
        END DO
     END DO
  END IF
  ! ... Calculate X, Y and Z velocities and mass fluxes at nodes
  xwvelmax = 0._kdp
  ywvelmax = 0._kdp
  zwvelmax = 0._kdp
  xsvelmax = 0._kdp
  ysvelmax = 0._kdp
  zsvelmax = 0._kdp
  xwmflxmax = 0._kdp
  ywmflxmax = 0._kdp
  zwmflxmax = 0._kdp
  xsmflxmax = 0._kdp
  ysmflxmax = 0._kdp
  zsmflxmax = 0._kdp
  ! ... Zero the velocity arrays
  xwvel = 0._kdp
  ywvel = 0._kdp
  zwvel = 0._kdp
  xsvel = 0._kdp
  ysvel = 0._kdp
  zsvel = 0._kdp
  xwmflx = 0._kdp
  ywmflx = 0._kdp
  zwmflx = 0._kdp
  xsmflx = 0._kdp
  ysmflx = 0._kdp
  zsmflx = 0._kdp
  DO  j = 1, ny
     DO  k = 1, nz
        DO  i = 1, nx
           ic = (k-1)*nx + i + (j-1)*nx*nz
           ik = (k-1)*nx + i
           ! ...       cell boundaries X direction
           mxm = (k-1)*nxx + i + (j-1)*nxx*nz
           mxp = mxm + 1
           ! ...        cell boundaries Y direction
           mym = (k-1)*nx + i + (j-1)*nx*nz
           myp = mym + nx*nz
           ! ...       cell boundaries Z direction
           mzm = (k-1)*nx + i + (j-1)*nx*nzz
           mzp = mzm + nx
           ! ... For all active nodes in mesh with fluid content
           IF (np(ik,j) /= 0 .AND. phi(ic) > 0._kdp) THEN
              ! ...       Average the fluxes and divide by porosity to get velocity
              ! ...       If either of the fluxes is 0 then the velocity is just 
              ! ...       the finite flux divided by porosity. This eliminates funny
              ! ...       looking velocity fields.
              ! ...       Go to the appropriate thermo region for the cell
              IF(ind(ic) == 1 .OR. satnw(ic) >= 1._kdp) THEN    ! ... Compressed Water Region
                 ! ...                 X direction
                 IF (nx > 1) THEN
                    IF (xwflux(mxm) == 0._kdp .OR. xwflux(mxp) == 0._kdp) THEN
                       xwvel(ic) = (xwflux(mxm)+xwflux(mxp))/phi(ic)
                       xwmflx(ic) = xwflux(mxm)*xdw(mxm) + xwflux(mxp)*xdw(mxp)
                    ELSE
                       xwvel(ic) = 0.5_kdp*(xwflux(mxm)+xwflux(mxp))/phi(ic)
                       xwmflx(ic) = 0.5_kdp*(xwflux(mxm)*xdw(mxm)+xwflux(mxp)*xdw(mxp))
                    END IF
                 END IF
                 ! ...                 Y direction
                 IF (ny > 1) THEN
                    IF (ywflux(mym) == 0._kdp .OR. ywflux(myp) == 0._kdp) THEN
                       ywvel(ic) = (ywflux(mym)+ywflux(myp))/phi(ic)
                       ywmflx(ic) = ywflux(mym)*ydw(mym) + ywflux(myp)*ydw(myp)
                    ELSE
                       ywvel(ic) = 0.5_kdp*(ywflux(mym)+ywflux(myp))/phi(ic)
                       ywmflx(ic) = 0.5_kdp*(ywflux(mym)*ydw(mym)+ywflux(myp)*ydw(myp))
                    END IF
                 END IF
                 ! ...                 Z direction
                 IF (nz > 1) THEN
                    IF (zwflux(mzm) == 0._kdp .OR. zwflux(mzp) == 0._kdp) THEN
                       zwvel(ic) = (zwflux(mzm)+zwflux(mzp))/phi(ic)
                       zwmflx(ic) = zwflux(mzm)*zdw(mzm) + zwflux(mzp)*zdw(mzp)
                    ELSE
                       zwvel(ic) = 0.5_kdp*(zwflux(mzm)+zwflux(mzp))/phi(ic)
                       zwmflx(ic) = 0.5_kdp*(zwflux(mzm)*zdw(mzm)+zwflux(mzp)*zdw(mzp))
                    END IF
                 END IF
              ELSE IF (ind(ic) == 2 .OR. ind(ic) == 4) THEN
                 ! ...     ****  Two-Phase Region or Super-Critical Fluid
                 ! ...                 X direction
                 IF (nx > 1) THEN
                    ! ...                  water
                    IF (xwflux(mxm) == 0._kdp .OR. xwflux(mxp) == 0._kdp) THEN
                       xwvel(ic) = (xwflux(mxm)+xwflux(mxp))/(phi(ic)*satnw(ic))
                       xwmflx(ic) = xwflux(mxm)*xdw(mxm) + xwflux(mxp)*xdw(mxp)
                    ELSE
                       xwvel(ic) = 0.5_kdp*(xwflux(mxm)+xwflux(mxp))/(phi(ic)*satnw(ic))
                       xwmflx(ic) = 0.5_kdp*(xwflux(mxm)*xdw(mxm)+xwflux(mxp)*xdw(mxp))
                    END IF
                    ! ...                  steam
                    IF (xsflux(mxm) == 0._kdp .OR. xsflux(mxp) == 0._kdp) THEN
                       xsvel(ic) = (xsflux(mxm)+xsflux(mxp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       xsmflx(ic) = xsflux(mxm)*xds(mxm) + xsflux(mxp)*xds(mxp)
                    ELSE
                       xsvel(ic) = 0.5_kdp*(xsflux(mxm)+xsflux(mxp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       xsmflx(ic) = 0.5_kdp*(xsflux(mxm)*xds(mxm)+xsflux(mxp)*xds(mxp))
                    END IF
                 END IF
                 ! ...                 Y direction
                 IF (ny > 1) THEN
                    ! ...                  water
                    IF (ywflux(mym) == 0._kdp .OR. ywflux(myp) == 0._kdp) THEN
                       ywvel(ic) = (ywflux(mym)+ywflux(myp))/(phi(ic)*satnw(ic))
                       ywmflx(ic) = ywflux(mym)*ydw(mym) + ywflux(myp)*ydw(myp)
                    ELSE
                       ywvel(ic) = 0.5_kdp*(ywflux(mym)+ywflux(myp))/(phi(ic)*satnw(ic))
                       ywmflx(ic) = 0.5_kdp*(ywflux(mym)*ydw(mym)+ywflux(myp)*ydw(myp))
                    END IF
                    ! ...                  steam
                    IF (ysflux(mym) == 0._kdp .OR. ysflux(myp) == 0._kdp) THEN
                       ysvel(ic) = (ysflux(mym)+ysflux(myp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       ysmflx(ic) = ysflux(mym)*yds(mym) + ysflux(myp)*yds(myp)
                    ELSE
                       ysvel(ic) = 0.5_kdp*(ysflux(mym)+ysflux(myp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       ysmflx(ic) = 0.5_kdp*(ysflux(mym)*yds(mym)+ysflux(myp)*yds(myp))
                    END IF
                 END IF
                 ! ...                 Z direction
                 IF (nz > 1) THEN
                    ! ...                  water
                    IF (zwflux(mzm) == 0._kdp .OR. zwflux(mzp) == 0._kdp) THEN
                       zwvel(ic) = (zwflux(mzm)+zwflux(mzp))/(phi(ic)*satnw(ic))
                       zwmflx(ic) = zwflux(mzm)*zdw(mzm) + zwflux(mzp)*zdw(mzp)
                    ELSE
                       zwvel(ic) = 0.5_kdp*(zwflux(mzm)+zwflux(mzp))/(phi(ic)*satnw(ic))
                       zwmflx(ic) = 0.5_kdp*(zwflux(mzm)*zdw(mzm)+zwflux(mzp)*zdw(mzp))
                    END IF
                    ! ...                  steam
                    IF (zsflux(mzm) == 0._kdp .OR. zsflux(mzp) == 0._kdp) THEN
                       zsvel(ic) = (zsflux(mzm)+zsflux(mzp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       zsmflx(ic) = zsflux(mzm)*zds(mzm) + zsflux(mzp)*zds(mzp)
                    ELSE
                       zsvel(ic) = 0.5_kdp*(zsflux(mzm)+zsflux(mzp))/  &
                            (phi(ic)*(1._kdp-satnw(ic)))
                       zsmflx(ic) = 0.5_kdp*(zsflux(mzm)*zds(mzm)+zsflux(mzp)*zds(mzp))
                    END IF
                 END IF
                 IF (ind(ic) == 4) THEN
                    ! ...       For supercritical fluid combine calculated water and
                    ! ...       steam velocities into the water velocity
                    ! ...       Note: supercritical nodes behave as though they have both steam
                    ! ...          and water with identical properties, the relative
                    ! ...          permeabilities of both phases are set to 0.5
                    xwvel(ic) = xwvel(ic) + xsvel(ic)
                    xsvel(ic) = 0._kdp
                    ywvel(ic) = ywvel(ic) + ysvel(ic)
                    ysvel(ic) = 0._kdp
                    zwvel(ic) = zwvel(ic) + zsvel(ic)
                    zsvel(ic) = 0._kdp
                    xwmflx(ic) = xwmflx(ic) + xsmflx(ic)
                    xsmflx(ic) = 0._kdp
                    ywmflx(ic) = ywmflx(ic) + ysmflx(ic)
                    ysmflx(ic) = 0._kdp
                    zwmflx(ic) = zwmflx(ic) + zsmflx(ic)
                    zsmflx(ic) = 0._kdp
                 END IF
              ELSE IF (ind(ic) == 3) THEN
                 ! ...     **  Super-Heated Steam Region
                 ! ...                 X direction
                 IF (nx > 1) THEN
                    IF (xsflux(mxm) == 0._kdp .OR. xsflux(mxp) == 0._kdp) THEN
                       xsvel(ic) = (xsflux(mxm)+xsflux(mxp))/phi(ic)
                       xsmflx(ic) = xsflux(mxm)*xds(mxm) + xsflux(mxp)*xds(mxp)
                    ELSE
                       xsvel(ic) = 0.5_kdp*(xsflux(mxm)+xsflux(mxp))/phi(ic)
                       xsmflx(ic) = 0.5_kdp*(xsflux(mxm)*xds(mxm)+xsflux(mxp)*xds(mxp))
                    END IF
                 END IF
                 ! ...                 Y direction
                 IF (ny > 1) THEN
                    IF (ysflux(mym) == 0._kdp .OR. ysflux(myp) == 0._kdp) THEN
                       ysvel(ic) = (ysflux(mym)+ysflux(myp))/phi(ic)
                       ysmflx(ic) = ysflux(mym)*yds(mym) + ysflux(myp)*yds(myp)
                    ELSE
                       ysvel(ic) = 0.5_kdp*(ysflux(mym)+ysflux(myp))/phi(ic)
                       ysmflx(ic) = 0.5_kdp*(ysflux(mym)*yds(mym)+ysflux(myp)*yds(myp))
                    END IF
                 END IF
                 ! ...                 Z direction
                 IF (nz > 1) THEN
                    IF (zsflux(mzm) == 0._kdp .OR. zsflux(mzp) == 0._kdp) THEN
                       zsvel(ic) = (zsflux(mzm)+zsflux(mzp))/phi(ic)
                       zsmflx(ic) = zsflux(mzm)*zds(mzm) + zsflux(mzp)*zds(mzp)
                    ELSE
                       zsvel(ic) = 0.5_kdp*(zsflux(mzm)+zsflux(mzp))/phi(ic)
                       zsmflx(ic) = 0.5_kdp*(zsflux(mzm)*zds(mzm)+zsflux(mzp)*zds(mzp))
                    END IF
                 END IF
              ELSE IF(ind(ic) == 5) THEN
                 ! ...  Air-Water Zone
                 ! ...                 X direction
                 IF (nx > 1) THEN
                    IF (xwflux(mxm) == 0._kdp .OR. xwflux(mxp) == 0._kdp) THEN
                       xwvel(ic) = (xwflux(mxm)+xwflux(mxp))/(phi(ic)*satnw(ic))
                       xwmflx(ic) = xwflux(mxm)*xdw(mxm) + xwflux(mxp)*xdw(mxp)
                    ELSE
                       xwvel(ic) = 0.5_kdp*(xwflux(mxm)+xwflux(mxp))/(phi(ic)*satnw(ic))
                       xwmflx(ic) = 0.5_kdp*(xwflux(mxm)*xdw(mxm)+xwflux(mxp)*xdw(mxp))
                    END IF
                 END IF
                 ! ...                 Y direction
                 IF (ny > 1) THEN
                    IF (ywflux(mym) == 0._kdp .OR. ywflux(myp) == 0._kdp) THEN
                       ywvel(ic) = (ywflux(mym)+ywflux(myp))/(phi(ic)*satnw(ic))
                       ywmflx(ic) = ywflux(mym)*ydw(mym) + ywflux(myp)*ydw(myp)
                    ELSE
                       ywvel(ic) = 0.5_kdp*(ywflux(mym)+ywflux(myp))/(phi(ic)*satnw(ic))
                       ywmflx(ic) = 0.5_kdp*(ywflux(mym)*ydw(mym)+ywflux(myp)*ydw(myp))
                    END IF
                 END IF
                 ! ...                 Z direction
                 IF (nz > 1) THEN
                    IF (zwflux(mzm) == 0._kdp .OR. zwflux(mzp) == 0._kdp) THEN
                       zwvel(ic) = (zwflux(mzm)+zwflux(mzp))/(phi(ic)*satnw(ic))
                       zwmflx(ic) = zwflux(mzm)*zdw(mzm) + zwflux(mzp)*zdw(mzp)
                    ELSE
                       zwvel(ic) = 0.5_kdp*(zwflux(mzm)+zwflux(mzp))/(phi(ic)*satnw(ic))
                       zwmflx(ic) = 0.5_kdp*(zwflux(mzm)*zdw(mzm)+zwflux(mzp)*zdw(mzp))
                    END IF
                 END IF
              END IF
           END IF
           ! ...     determine the maximum magnitudes of X,Y & Z velocities for water
           xwvelmax = MAX(xwvelmax,ABS(xwvel(ic)))
           ywvelmax = MAX(ywvelmax,ABS(ywvel(ic)))
           zwvelmax = MAX(zwvelmax,ABS(zwvel(ic)))
           ! ...     determine the maximum magnitudes of X,Y & Z velocities for steam
           xsvelmax = MAX(xsvelmax,ABS(xsvel(ic)))
           ysvelmax = MAX(ysvelmax,ABS(ysvel(ic)))
           zsvelmax = MAX(zsvelmax,ABS(zsvel(ic)))
           ! ...     determine the maximum magnitudes of X,Y & Z mass fluxes for water
           xwmflxmax = MAX(xwmflxmax,ABS(xwmflx(ic)))
           ywmflxmax = MAX(ywmflxmax,ABS(ywmflx(ic)))
           zwmflxmax = MAX(zwmflxmax,ABS(zwmflx(ic)))
           ! ...     determine the maximum magnitudes of X,Y & Z mass fluxes for steam
           xsmflxmax = MAX(xsmflxmax,ABS(xsmflx(ic)))
           ysmflxmax = MAX(ysmflxmax,ABS(ysmflx(ic)))
           zsmflxmax = MAX(zsmflxmax,ABS(zsmflx(ic)))
        END DO
     END DO
  END DO
END SUBROUTINE velocity
