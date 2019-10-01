SUBROUTINE printar(ndim,array,lprnt,fu,cnv,jfmt,nnoppr)
  !.....Prints out the contents of 2-dimensional array with selected
  !.....     format
  !.....ARRAY - The full 3-dimensional array
  ! ...      Stored in x-z slice order of node numbering
  !.....CNV - Factor for units conversion; HT usage is the inverse, from
  !.....      input to internal units
  !.....Selects layers or slices, single or multiple
  !.....NNOPPR - Number of nodes or planes printed
  !.....     for NDIM=1
  !.....     number of nodes printed
  !.....     for NDIM=2
  !.....     -100 - NX-1 planes; -010 - NY-1 planes; -001 - NZ-1 planes
  !.....     100 - One X-plane; 010 - One Y-plane; 001 - One Z-plane
  !.....     200 - NX+1 planes; 020 - NY+1 planes; 002 - NZ+1 planes
  !.....     000 - NX, NY, and NZ planes are printed in their respective directions
  ! ... Presently prints only x-z planes of data for Hydrotherm
  USE machine_constants, ONLY: kdp
!  USE control
  USE mesh
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndim
  REAL(kind=kdp), DIMENSION(:), INTENT(IN) :: array
  INTEGER, DIMENSION(:,:), INTENT(IN) :: lprnt
  INTEGER, INTENT(IN) :: fu
  REAL(kind=kdp), INTENT(IN) :: cnv
  INTEGER, INTENT(IN) :: jfmt
  INTEGER, INTENT(IN) :: nnoppr
  !
  CHARACTER (LEN=1) :: ir3lbl
  CHARACTER (LEN=9) :: form
  CHARACTER (LEN=12) :: sp12='            '
  CHARACTER (LEN=9), DIMENSION(8) :: aform=(/'  (F12.0)','  (F12.1)','  (F12.2)',  &
       '  (F12.3)','  (F12.4)','  (F12.5)','(1PG12.4)','(1PG12.5)'/)
  CHARACTER (LEN=12), DIMENSION(10) :: caprnt
  INTEGER :: i, iifmt, iir2, ior, ir2, ir3, ir3p, mic, mik, n1, n2, n3, nnpr,  &
       npr2, npr3, nxpr, nypr, nzpr, nxypr, nxzpr
  INTEGER, PARAMETER :: orenpr=13
  !.....ORENPR: 12 - Areal layers; 13 - Vertical slices
  LOGICAL :: prwin
  REAL(kind=kdp) :: cnvi
  REAL(kind=kdp), DIMENSION(10) :: aprnt
  !.....Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2003/01/22 22:00:03 $'
  !-----------------------------------------------------------------------------
  cnvi = 1._kdp/cnv
  ior=MOD(ABS(orenpr),10)
  SELECT CASE (jfmt)
     CASE (10)
        iifmt=1
     CASE (11)
        iifmt=2
     CASE (12)
        iifmt=3
     CASE (13)
        iifmt=4
     CASE (14)
        iifmt=5
     CASE (15)
        iifmt=6
     CASE (24)
        iifmt=7
     CASE (25)
        iifmt=8
  END SELECT
  form=aform(iifmt)
  IF(ndim == 1) THEN
     nxpr=nnoppr
     npr3=1
  ELSE
     SELECT CASE (nnoppr/100)
        CASE (-1)
           nxpr = nx-1
        CASE (0)
           nxpr = nx
        CASE (1)
           nxpr = 1
        CASE (2)
           nxpr = nx+1
     END SELECT
     SELECT CASE (MOD(nnoppr,100)/10)
        CASE (-1)
           nypr = ny-1
        CASE (0)
           nypr = ny
        CASE (1)
           nypr = 1
        CASE (2)
           nypr = ny+1
     END SELECT
     SELECT CASE (MOD(nnoppr,10))
        CASE (-1)
           nzpr = nz-1
        CASE (0)
           nzpr = nz
        CASE (1)
           nzpr = 1
        CASE (2)
           nzpr = nz+1
     END SELECT
     nxypr = nxpr*nypr
     nxzpr = nxpr*nzpr
     IF(ior == 2) THEN
        !.....Areal layers
        ir3lbl='K'
        npr3=nzpr
        npr2=nypr
        IF(npr3 > 1) WRITE(fu,2001) 'Areal Layers'
2001    FORMAT(/tr20,a)
     ELSE IF(ior == 3) THEN
        !.....Vertical slices
        ir3lbl='J'
        npr3=nypr
        npr2=nzpr
        IF(npr3 > 1) WRITE(fu,2001) 'Vertical Slices'
     END IF
  END IF
  DO  ir3=1,npr3
     IF(ndim == 2) THEN
        ! ... Do not print this plane if all cells are excluded
        nnpr=0
        DO  ir2=1,npr2
           DO  i=1,nxpr
!              IF(ior == 2) m=(ir3-1)*nxy+(ir2-1)*nx+i
!              IF(ior == 3) m=(ir2-1)*nxy+(ir3-1)*nx+i
              mik=(ir2-1)*nxpr+i
              IF(lprnt(mik,ir3) /= -1) THEN
                 nnpr=1
                 EXIT
              END IF
           END DO
        END DO
        IF(nnpr == 0) CYCLE     ! ... No active cells in this plane
        ir3p=ir3
        IF(MOD(nnoppr,10) == 1) ir3p=nz
        IF(npr3 > 1) WRITE(fu,2003) ir3lbl,' =',ir3p
2003    FORMAT(/tr60,a1,a,i3)
     END IF
     n1=1
     n2=10
50   n2=MIN(n2,nxpr)
     WRITE(fu,2004) (i,i=n1,n2)
2004 FORMAT(/i12,9I12)
     n3=n2-n1+1
     IF(ndim == 1) THEN
        DO  i=1,n3
           aprnt(i)=cnvi*array(n1-1+i)
           WRITE(caprnt(i),form) aprnt(i)
        END DO
        WRITE(fu,2005) (caprnt(i),i=1,n3)
2005    FORMAT(tr3,10A12)
     ELSE IF(ndim == 2) THEN
        prwin=.false.
        DO  ir2=1,npr2
           iir2=npr2+1-ir2
           DO  i=1,n3
!              IF(ior == 2) m=(ir3-1)*nxy+(iir2-1)*nx+n1-1+i
!              IF(ior == 3) m=(iir2-1)*nxy+(ir3-1)*nx+n1-1+i
              mik=(iir2-1)*nxpr+n1-1+i
              IF(lprnt(mik,ir3) /= -1) THEN
                 prwin=.true.     ! ... At least one active cell in this plane
                 EXIT
              END IF
           END DO
        END DO
        IF(prwin) THEN
           DO ir2=1,npr2
              iir2=npr2+1-ir2
              DO i=1,n3
                 IF(ior == 2) mic=(iir2-1)*nxzpr+(ir3-1)*nxpr+n1-1+i
                 IF(ior == 3) mic=(ir3-1)*nxzpr+(iir2-1)*nxpr+n1-1+i
                 mik=(iir2-1)*nxpr+n1-1+i
                 IF(lprnt(mik,ir3) /= -1) THEN
                    aprnt(i)=cnvi*array(mic)
                    WRITE(caprnt(i),form) aprnt(i)
                 ELSE
                    WRITE(caprnt(i),3001) sp12
3001                FORMAT(a12)
                 END IF
              END DO
              IF(nnpr == 0) THEN
                 WRITE(fu,2006) iir2
2006             FORMAT(i3)
              ELSE
                 WRITE(fu,2007) iir2,(caprnt(i),i=1,n3)
2007             FORMAT(i3,10A12)
              END IF
           END DO
        END IF
     END IF
     IF(n2 == nxpr) CYCLE
     n1=n1+10
     n2=n2+10
     WRITE(fu,2008)
2008 FORMAT(/)
     GO TO 50
  END DO
  WRITE(fu,2008)
END SUBROUTINE printar
