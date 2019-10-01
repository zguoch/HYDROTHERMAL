SUBROUTINE tcalc
  ! ... Purpose:  To compute the harmonic mean of the permeability and
  ! ...           thermal conductivity coefficients at the cell faces.
  ! ...           These coefficients contain the cross-sectional area
  ! ...           and the distance between the nodes.
  ! ... Internal Units:  cm^3 for TX,TY,TZ
  ! ...         erg/(s-deg.C) for TXK,TYK,TZK
  USE machine_constants, ONLY: kdp
  USE math_constants
  USE fdeq
  USE mesh
  USE parameters
  USE variables, ONLY: rsq, xx
  IMPLICIT NONE
  INTEGER :: i, ic, icxm, icym, iczm, ik, j, k, mxf, mxl, mxm, mxp,  &
       myf, myl, mym, myp, mzf, mzl, mzm, mzp
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.8 $//$Date: 2006/10/03 23:39:43 $'
  !     ------------------------------------------------------------------
  ! ... Isotropic conduction for now; make ykc=zkc=xkc
  ykc = xkc
  zkc = xkc
  IF (.NOT.irad) THEN
     ! ... Rectangular coordinates
     ! ...      X-direction 
     IF (nx > 1) THEN 
        DO  j = 1,ny
           DO  k = 1,nz
              DO  i = 2,nx
                 ic = (k-1)*nx + i + (j-1)*nxz
                 icxm = ic - 1
                 mxm = (k-1)*nxx + i + (j-1)*nxx*nz
                 ax(mxm) = dx(i)/(dx(i)+dx(i-1))     ! ... fraction of internode distance to face
                 ! ...        compute harmonic mean of permeability term if both Xperms not eq 0
                 IF (xk(icxm) > 0.0 .AND. xk(ic) > 0.0) THEN
                    tx(mxm) = 2._kdp*dz(k)*dy(j)*xk(icxm)*xk(ic)/  &
                         (xk(icxm)*dx(i)+xk(ic)*dx(i-1))
                 ELSE
                    tx(mxm) = 0._kdp
                 END IF
                 ! ...        compute harmonic mean of heat conduction term if both XKC's ne. 0
                 IF (xkc(icxm) > 0.0 .AND. xkc(ic) > 0.0) THEN
                    txk(mxm) = 2._kdp*dz(k)*dy(j)*xkc(icxm)*xkc(ic)/  &
                         (xkc(icxm)*dx(i)+xkc(ic)*dx(i-1))
                 ELSE
                    txk(mxm) = 0._kdp
                 END IF
              END DO
           END DO
        END DO
     END IF
     ! ...  Z-direction
     IF (nz > 1) THEN
        DO  j = 1,ny
           DO  i = 1,nx
              DO  k = 2,nz
                 ic = (k-1)*nx + i + (j-1)*nxz
                 iczm = ic - nx
                 mzm = (k-1)*nx + i + (j-1)*nx*nzz
                 az(mzm) = dz(k)/(dz(k)+dz(k-1))     ! ... fraction of internode distance to face
                 ! ...        compute harmonic mean of permeability term if both Zperms not = 0
                 IF (zk(iczm) > 0.0 .AND. zk(ic) > 0.0) THEN
                    tz(mzm) = 2._kdp*dx(i)*dy(j)*zk(iczm)*zk(ic)/  &
                         (zk(iczm)*dz(k)+zk(ic)*dz(k-1))
                 ELSE
                    tz(mzm) = 0._kdp
                 END IF
                 ! ...        compute harmonic mean of heat conduction term if both ZKC's ne. 0
                 IF (zkc(iczm) > 0.0 .AND. zkc(ic) > 0.0) THEN
                    tzk(mzm) = 2._kdp*dx(i)*dy(j)*zkc(iczm)*zkc(ic)/  &
                         (zkc(iczm)*dz(k)+zkc(ic)*dz(k-1))
                 ELSE
                    tzk(mzm) = 0._kdp
                 END IF
              END DO
           END DO
        END DO
     END IF
     ! ... Y-direction
     IF (ny > 1) THEN
        DO  i = 1,nx
           DO  k = 1,nz
              DO  j = 2,ny
                 ic = (k-1)*nx + i + (j-1)*nxz
                 icym = ic - nxz
                 mym = ic
                 ay(mym) = dy(j)/(dy(j)+dy(j-1))     ! ... fraction of internode distance to face
                 ! ...        compute harmonic mean of permeability term if both Yperms not = 0
                 IF (yk(icym) > 0.0 .AND. yk(ic) > 0.0) THEN
                    ty(mym) = 2._kdp*dx(i)*dz(k)*yk(icym)*yk(ic)/  &
                         (yk(icym)*dy(j)+yk(ic)*dy(j-1))
                 ELSE
                    ty(mym) = 0._kdp
                 END IF
                 ! ...        compute harmonic mean of heat conduction term if both YKC's ne. 0
                 IF (ykc(icym) > 0.0 .AND. ykc(ic) > 0.0) THEN
                    tyk(mym) = 2._kdp*dx(i)*dz(k)*ykc(icym)*ykc(ic)/  &
                         (ykc(icym)*dy(j)+ykc(ic)*dy(j-1))
                 ELSE
                    tyk(mym) = 0._kdp
                 END IF
              END DO
           END DO
        END DO
     END IF
  ELSE
     ! ...  Radial-direction
     j = 1
     DO  k = 1,nz
        DO  i = 2,nx
           ic = (k-1)*nx + i + (j-1)*nxz
           icxm = ic - 1
           mxm = (k-1)*nxx + i
           ax(mxm) = (xx(i)**2-rsq(i))/(xx(i)**2-xx(i-1)**2)     ! ... fraction of internode
           ! ...    distance to face
           ! ...        compute harmonic mean of permeability term if both Xperms not = 0
           IF (xk(icxm) > 0.0 .AND. xk(ic) > 0.0) THEN
              tx(mxm) = 4._kdp/(xx(i)**2-xx(i-1)**2)*pi*rsq(i)*dz(k)*  &
                   (2._kdp*xk(icxm)*xk(ic))/(xk(icxm)+xk(ic))
           ELSE
              tx(mxm) = 0._kdp
           END IF
           ! ...        compute harmonic mean of heat conduction term if both XKC's ne. 0
           IF (xkc(icxm) > 0.0 .AND. xkc(ic) > 0.0) THEN
              txk(mxm) = 4._kdp/(xx(i)**2-xx(i-1)**2)*pi*rsq(i)*dz(k)*  &
                   (2._kdp*xkc(icxm)*xkc(ic))/(xkc(icxm)+xkc(ic))
           ELSE
              txk(mxm) = 0._kdp
           END IF
        END DO
     END DO
     ! ... Z-direction
     IF (nz > 1) THEN
        DO  i = 1,nx
           DO  k = 2,nz
              ic = (k-1)*nx + i
              iczm = ic - nx
              mzm = (k-1)*nx + i
              !..                m=(k-1)*NXY + (j-1)*NX+ i
              az(mzm) = dz(k)/(dz(k)+dz(k-1))     ! ... fraction of internode distance to face
              ! ...        compute harmonic mean of permeability term if both Zperms not = 0
              IF (zk(iczm) > 0.0 .AND. zk(ic) > 0.0) THEN
                 tz(mzm) = 2._kdp*zk(iczm)*zk(ic)*pi*(rsq(i+1)-rsq(i))  &
                      /(zk(iczm)*dz(k)+zk(ic)*dz(k-1))
              ELSE
                 tz(mzm) = 0._kdp
              END IF
              ! ...        compute harmonic mean of heat conduction term if both ZKC's ne. 0
              IF (zkc(iczm) > 0.0 .AND. zkc(ic) > 0.0) THEN
                 tzk(mzm) = 2._kdp*zkc(iczm)*zkc(ic)*  &
                      pi*(rsq(i+1)-rsq(i))/(zkc(iczm)*dz(k)+zkc(ic)*dz(k-1))
              ELSE
                 tzk(mzm) = 0._kdp
              END IF
           END DO
        END DO
     END IF
  END IF
  ! ... set transmissibility of boundary cells for a global rectangular region
  ! ...      to zero; no flow
  ! ...     Z boundary cella
  DO  j = 1,ny
     DO  i = 1,nx
        mzf = i + (j-1)*nzz*nx
        mzl = nz*nx + i + (j-1)*nzz*nx
        tz(mzf) = 0._kdp
        tz(mzl) = 0._kdp
        tzk(mzf) = 0._kdp
        tzk(mzl) = 0._kdp
     END DO
  END DO
  ! ...     X boundary cells
  DO  j = 1,ny
     DO  k = 1,nz
        mxf = (k-1)*nxx + 1 + (j-1)*nxx*nz
        mxl = (k-1)*nxx + nxx + (j-1)*nxx*nz
        tx(mxf) = 0._kdp
        tx(mxl) = 0._kdp
        txk(mxf) = 0._kdp
        txk(mxl) = 0._kdp
     END DO
  END DO
  ! ...     Y boundary cells
  IF (ny > 1) THEN
     DO  k = 1, nz
        DO  i = 1, nx
           myf = (k-1)*nx + i
           myl = myf + nx*nz*ny
           ty(myf) = 0._kdp
           ty(myl) = 0._kdp
           tyk(myf) = 0._kdp
           tyk(myl) = 0._kdp
        END DO
     END DO
  END IF
  ! ...   set transmissibility of cell faces forming external boundaries of region
  ! ...        to zero
  DO  ic=1,nxyz
     CALL ictoijk(ic,i,j,k,nx,nz)
     ik = (k-1)*nx + i
     !..                m=(k-1)*NXY + (j-1)*NX+ i
     IF (np(ik,j) == 0) THEN
        mxm = (k-1)*nxx + i + (j-1)*nxx*nz
        mxp = mxm + 1
        mzm = (k-1)*nx + i + (j-1)*nx*nzz
        mzp = mzm + nx
        mym = (k-1)*nx + i + (j-1)*nxz
        myp = mym + nxz
        tx(mxm) = 0._kdp
        tx(mxp) = 0._kdp
        txk(mxm) = 0._kdp
        txk(mxp) = 0._kdp
        tz(mzm) = 0._kdp
        tz(mzp) = 0._kdp
        tzk(mzm) = 0._kdp
        tzk(mzp) = 0._kdp
        IF (ny > 1) THEN
           ty(mym) = 0._kdp
           ty(myp) = 0._kdp
           tyk(mym) = 0._kdp
           tyk(myp) = 0._kdp
        END IF
     END IF
  END DO
END SUBROUTINE tcalc
