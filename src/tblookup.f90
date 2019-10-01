SUBROUTINE tblookup(px, hx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp, errflag)
  ! ... Purpose:  to call the appropriate lookup table for water properties
  ! ...           in the range (excluding the two-phase field):
  ! ...               0.01e9 <= H <= 52.0e9 erg/g
  ! ...                0.5 <= P < 10,000 bar
  USE machine_constants, ONLY: kdp
  USE f_units
  USE bc, ONLY: unconfined
  USE control
  USE tables
  IMPLICIT NONE
  ! ...   HX - enthalpy [erg/g]
  ! ...   PX -  pressure [dyne/cm^2]
  ! ...   Da -  density of water [g/cm^3]
  ! ...      if Da<0 then set Da to -HERGS1(ilo)
  ! ...   PDaPH -  derivative of density w.r.t. enthalpy
  ! ...   PDaPP -  derivative of density w.r.t. pressure
  ! ...   Ta -  temperature of water
  ! ...   PTaPH -  derivative of temperature w.r.t. enthalpy
  ! ...   PTaPP -  derivative of temperature w.r.t. pressure
  ! ...   Va -  viscosity of water
  ! ...   PVaPH -  derivative of viscosity w.r.t. enthalpy
  ! ...   PVaPP -  derivative of viscosity w.r.t. pressure
  REAL(KIND=kdp), INTENT(IN) :: hx, px
  REAL(KIND=kdp), INTENT(OUT) :: da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp
  LOGICAL, INTENT(INOUT) :: errflag
  !
  REAL(KIND=kdp) :: pxx
  REAL(KIND=kdp) :: db, dum1, dum2, dum3, dum4, dum5, dum6, dum7,  &
       dum8, dum9, dumy1, dumy2, dumy3, hh, hl, hs, ht0,  &
       hw, pdbph, pdbpp, pssd1ph, pssd1pp, psst1ph,  &
       psst1pp, pssv1ph, pssv1pp, pswd1ph, pswd1pp,  &
       pswt1ph, pswt1pp, pswv1ph, pswv1pp, pt0d1ph,  &
       pt0d1pp, pt0dph1, pt0dpp1, pt0t1ph, pt0t1pp,  &
       pt0tph1, pt0tpp1, pt0v1ph, pt0v1pp, pt0vph1,  &
       pt0vpp1, ptbph, ptbpp, pvbph, pvbpp, ssd1, sst1,  &
       ssv1, swd1, swt1, swv1, t0d1, t0t1, t0v1, tb, vb
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.10 $//$Date: 2007/08/03 23:42:47 $'
  !     ------------------------------------------------------------------
  errflag = .FALSE.
  pxx = px
  ! ... Check that input P & H are within the limits of the table
  ! ...      You should not be here for air-water zone with pressure less than 0.5 bar
  ! ... ***** but at the present state of code, you can be
  ! ... Special patch to allow for low pressures; unconfined flow
  IF (px < 0.5e6_kdp) THEN
     IF(unconfined) THEN
        itmcntl(20) = itmcntl(20)+1
        IF(MOD(itmcntl(20),1000) == 1) THEN
           WRITE (fustdout,9005)  &
                ' ***** Pressure too low for limits of the table, for P= ',  &
                px, ' (dyne/cm^2);','Low limit will be used; occurrence ',  &
                itmcntl(20),' *****'
           WRITE (fuclog,9005)  &
                ' ***** Pressure too low for limits of the table, for P= ',  &
                px, ' (dyne/cm^2);','Low limit will be used; occurrence ',  &
                itmcntl(20),' *****'
9005       FORMAT (a,1pg14.7,a/tr10,a,i8,a)
        END IF
        pxx = 0.5e6_kdp
     ELSE
        WRITE (fustdout,9005)  &
             ' *****ERROR: Pressure too low for limits of the table, for P= ',  &
             px, ' (dyne/cm^2) *****'
        WRITE (fuclog,9005)  &
             ' *****ERROR: Pressure too low for limits of the table, for P= ',  &
             px, ' (dyne/cm^2) *****'
        ierr(135) = .TRUE.
        errflag = .TRUE.
        RETURN
     END IF
  ELSE IF (px > 1.0e10_kdp) THEN
     WRITE (fustdout, 9005)  &
          ' *****ERROR: Pressure too high for limits of the table, for P= ',px,'(dyne/cm^2) *****'
     WRITE (fuclog, 9005)  &
          ' *****ERROR: Pressure too high for limits of the table, for P= ',px,'(dyne/cm^2) *****'
     ierr(123) = .TRUE.
     errflag = .TRUE.
     RETURN
  END IF
  IF (hx > 52.0e9_kdp) THEN
     WRITE (fustdout, 9005)  &
          ' *****ERROR: Enthalpy too high for limits of the table, for H= ',hx,'(erg/g) *****'
     WRITE (fuclog, 9005)  &
          ' *****ERROR: Enthalpy too high for limits of the table, for H= ',hx,'(erg/g) *****'
     ierr(124) = .TRUE.
     errflag = .TRUE.
     RETURN
  END IF
  IF ((px < 240.e6_kdp .AND. hx < 0.01e9_kdp) .OR.  &
       (px >= 240.0e6_kdp .AND. hx < 0.5e9_kdp)) THEN
     WRITE (fustdout, 9006)  &
          ' *****ERROR: Enthalpy too low for limits of the table, for H= ' , hx, '(erg/g)',  &
          ' and P= ',px,'(dyne/cm^2) *****'
     WRITE (fuclog, 9006)  &
          ' *****ERROR: Enthalpy too low for limits of the table, for H= ' , hx, '(erg/g)',  &
          ' and P= ',px,'(dyne/cm^2) *****'
9006 FORMAT (a,1PG14.7,a/tr10,a,1pg14.7,a)
     ierr(132) = .TRUE.
     errflag = .TRUE.
     RETURN
  END IF
  ! ...     for Pressures < 240 bars
  IF (pxx < 240.0e6_kdp) THEN
     IF (hx < 16.0e9_kdp) THEN
        ! ... For P < 240 bar and H < 16.0e9 erg/g
        CALL table1(hx, pxx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
        IF (da <= 0._kdp) THEN           ! ... Point is between table and 2-phase curve
           db = da
           DO WHILE (db <= 0._kdp)
              hl = -db
              CALL table1(hl, pxx, db, pdbph, pdbpp, tb, ptbph, ptbpp, vb, pvbph, pvbpp)
           END DO
           CALL phsbndy2(pxx, hw, dumy1, dumy2, dumy3,  &
                swd1, pswd1ph, pswd1pp, swt1, pswt1ph, pswt1pp,  &
                swv1, pswv1ph, pswv1pp,  &
                dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9)
           CALL spline42(hw, swd1, pswd1ph, pswd1pp, swt1, pswt1ph, pswt1pp,  &
                swv1, pswv1ph, pswv1pp, hl, db, pdbph, pdbpp, tb, ptbph, ptbpp,  &
                vb, pvbph, pvbpp, hx, da, pdaph, pdapp, ta, ptaph, ptapp,  &
                va, pvaph, pvapp)
        END IF
     ELSE
        IF (hx < 26.0e9_kdp) THEN
           ! ...             for P < 240 bar and 16.0D9 < H < 26.0D9 ergs/g
           CALL table3(hx, pxx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
           IF (da <= 0._kdp) THEN         ! ... Point is between table and 2-phase curve
              IF (hx <= 20.86e9_kdp) THEN
                 ! ... for points to the left of the two-phase field
                 ! ...      and right of tabulated values
                 db = da
                 DO WHILE (db <= 0._kdp)
                    hl = -db
                    IF (hl < 16.0e9_kdp) THEN
                       CALL table1(hl, pxx, db, pdbph, pdbpp,  &
                            tb, ptbph, ptbpp, vb, pvbph, pvbpp)
                    ELSE
                       CALL table3(hl, pxx, db, pdbph, pdbpp,  &
                            tb, ptbph, ptbpp, vb, pvbph, pvbpp)
                    END IF
                 END DO
                 CALL phsbndy2(pxx, hw, dumy1, dumy2, dumy3,  &
                      swd1, pswd1ph, pswd1pp, swt1, pswt1ph, pswt1pp,  &
                      swv1, pswv1ph, pswv1pp,  &
                      dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9)
                 CALL spline42(hw, swd1, pswd1ph, pswd1pp, swt1, pswt1ph, pswt1pp,  &
                      swv1, pswv1ph, pswv1pp, hl, db, pdbph, pdbpp, tb, ptbph, ptbpp,  &
                      vb, pvbph, pvbpp, hx, da, pdaph, pdapp, ta, ptaph, ptapp,  &
                      va, pvaph, pvapp)
              ELSE
                 ! ... for points to the right of the two-phase field
                 ! ...      and left of tabulated values
                 db = da
                 DO WHILE (db <= 0._kdp)
                    hh = -db
                    IF (hh < 26.0e9_kdp) THEN
                       CALL table3(hh, pxx, db, pdbph, pdbpp,  &
                            tb, ptbph, ptbpp, vb, pvbph, pvbpp)
                    ELSE
                       CALL table2(hh, pxx, db, pdbph, pdbpp,  &
                            tb, ptbph, ptbpp, vb, pvbph, pvbpp)
                    END IF
                 END DO
                 CALL phsbndy2(pxx, dumy1, hs, dumy2, dumy3,  &
                      dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9,  &
                      ssd1, pssd1ph, pssd1pp, sst1, psst1ph, psst1pp,  &
                      ssv1, pssv1ph, pssv1pp)
                 CALL spline42(hh, db, pdbph, pdbpp, tb, ptbph, ptbpp,  &
                      vb, pvbph, pvbpp, hs, ssd1, pssd1ph, pssd1pp,  &
                      sst1, psst1ph, psst1pp, ssv1, pssv1ph, pssv1pp,  &
                      hx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
              END IF
           END IF
        ELSE
           ! ... for P < 240 bar and 26.0D9 < H < 52.0D9 ergs/g
           CALL table2(hx, pxx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
           IF (da <= 0._kdp) THEN         ! ... point is between table and 2-phase curve
              db = da
              DO WHILE (db <= 0._kdp)
                 hh = -db
                 CALL table2(hh, pxx, db, pdbph, pdbpp,  &
                      tb, ptbph, ptbpp, vb, pvbph, pvbpp)
              END DO
              CALL phsbndy2(pxx, dumy1, hs, dumy2, dumy3,  &
                   dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9,  &
                   ssd1, pssd1ph, pssd1pp, sst1, psst1ph, psst1pp,  &
                   ssv1, pssv1ph, pssv1pp)
              CALL spline42(hh, db, pdbph, pdbpp, tb, ptbph, ptbpp,  &
                   vb, pvbph, pvbpp, hs, ssd1, pssd1ph, pssd1pp,  &
                   sst1, psst1ph, psst1pp, ssv1, pssv1ph, pssv1pp,  &
                   hx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
           END IF
        END IF
     END IF
  ELSE
     ! ... for Pressures >= 240 bars
     CALL table4(hx, pxx, da, pdaph, pdapp, ta, ptaph, ptapp, va, pvaph, pvapp)
     IF (da <= 0._kdp) THEN
        ! ... point is between table and Temperature = 0 C curve
        ! ... for points to the right of the T=0 field
        ! ...      and left of tabulated values
        db = da
        DO WHILE (db <= 0._kdp)
           hh = -db
           CALL table4(hh, pxx, db, pdbph, pdbpp, tb, ptbph, ptbpp, vb, pvbph, pvbpp)
        END DO
        CALL bndyt0(pxx, ht0, t0d1, pt0dph1, pt0dpp1, t0t1, pt0tph1,  &
             pt0tpp1, t0v1, pt0vph1, pt0vpp1)
        pt0d1ph = 0._kdp
        pt0d1pp = 0._kdp
        pt0t1ph = 0._kdp
        pt0t1pp = 0._kdp
        pt0v1ph = 0._kdp
        pt0v1pp = 0._kdp
        CALL spline42(hh, db, pdbph, pdbpp, tb, ptbph, ptbpp, vb, pvbph, pvbpp,  &
             ht0, t0d1, pt0d1ph, pt0d1pp, t0t1, pt0t1ph, pt0t1pp,  &
             t0v1, pt0v1ph, pt0v1pp, hx, da, pdaph, pdapp, ta, ptaph, ptapp,  &
             va, pvaph, pvapp)
     END IF
  END IF
END SUBROUTINE tblookup
