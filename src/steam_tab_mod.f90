MODULE steam_table_funct
  ! ... Set of data and routines to determine enthalpy and density of
  ! ...      water and steam as a function of pressure and temperature
  USE machine_constants, ONLY: kdp
  IMPLICIT NONE
  PRIVATE     ! Makes everything private
  PUBLIC :: tempress     !  The only accessible routine in this module
  INTEGER, PARAMETER :: nc=36
!!$  INTEGER, DIMENSION(40), PARAMETER :: ii=(/4*0, 4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*8,   &
!!$       2*2, 0, 4, 3*2, 4/),  &
!!$       jj=(/2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5,  &
!!$       7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 1, 3*4, 0, 2, 0, 0/)
  INTEGER, DIMENSION(40), PARAMETER :: ii=(/0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4,  &
       5,5,5,5, 6,6,6,6, 8,8,8,8,   &
       2,2, 0, 4, 2,2,2, 4/),  &
       jj=(/2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5,  &
       7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 1, 4,4,4, 0, 2, 0, 0/)
  REAL(KIND=kdp) :: q0, q5
  REAL(KIND=kdp) :: ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, cjtt, cjth
  REAL(KIND=kdp), PARAMETER :: gascon=0.461522_kdp, tzren=647.073_kdp, aa=1._kdp,  &
       uref=-4328.4549774261_kdp, sref=7.6180720166752_kdp
  REAL(KIND=kdp) :: z, dzren, y
  REAL(KIND=kdp), DIMENSION(40), PARAMETER :: gren=(/  &
       -0.53062968529023E3_kdp, 0.22744901424408E4_kdp, 0.78779333020687E3_kdp,   &
       -0.69830527374994E2_kdp, 0.17863832875422E5_kdp, -0.39514731563338E5_kdp,  &
       0.33803884280753E5_kdp, -0.13855050202703E5_kdp, -0.25637436613260E6_kdp,  &
       0.48212575981415E6_kdp, -0.34183016969660E6_kdp, 0.12223156417448E6_kdp,   &
       0.11797433655832E7_kdp, -0.21734810110373E7_kdp, 0.10829952168620E7_kdp,   &
       -0.25441998064049E6_kdp, -0.31377774947767E7_kdp, 0.52911910757704E7_kdp,  &
       -0.13802577177877E7_kdp, -0.25109914369001E6_kdp, 0.46561826115608E7_kdp,  &
       -0.72752773275387E7_kdp, 0.41774246148294E6_kdp, 0.14016358244614E7_kdp,   &
       -0.31555231392127E7_kdp, 0.47929666384584E7_kdp, 0.40912664781209E6_kdp,   &
       -0.13626369388386E7_kdp, 0.69625220862664E6_kdp, -0.10834900096447E7_kdp,  &
       -0.22722827401688E6_kdp, 0.38365486000660E6_kdp, 0.68833257944332E4_kdp,   &
       0.21757245522644E5_kdp, -0.26627944829770E4_kdp, -0.70730418082074E5_kdp,  &
       -0.225_kdp, -1.68_kdp, 0.055_kdp, -93._kdp/)
  REAL(KIND=kdp), PARAMETER :: g1=11._kdp, g2=44.333333333333_kdp, gf=3.5_kdp
  REAL(KIND=kdp) :: b1, b2, b1t, b2t, b1tt, b2tt
  REAL(KIND=kdp), DIMENSION(10), PARAMETER :: bp=(/0.7478629_kdp, -0.3540782_kdp, &
       0._kdp, 0._kdp,  &
       0.7159876e-2_kdp, 0._kdp, -0.3528426e-2_kdp, 0._kdp, 0._kdp, 0._kdp/),  &
       bq=(/1.1278334_kdp, 0._kdp, -0.5944001_kdp, -5.010996_kdp, 0._kdp, 0.63684256_kdp,  &
       0._kdp, 0._kdp, 0._kdp, 0._kdp/)
  REAL(KIND=kdp), DIMENSION(4), PARAMETER :: adz=(/0.319_kdp, 0.319_kdp, 0.319_kdp, &
       1.55_kdp/),  &
       atz=(/64.e1_kdp, 64.e1_kdp, 641.6_kdp, 27.e1_kdp/),  &
       aad=(/34._kdp, 4.e1_kdp, 3.e1_kdp, 1.05e3_kdp/), aat=(/2.e4_kdp, 2.e4_kdp, &
       4.e4_kdp, 25._kdp/)
  REAL(KIND=kdp) :: ab, gb, sb, ub, hb, cvb, dpdtb
  REAL(KIND=kdp) :: ai, gi, si, ui, hi, cvi, cpi
  REAL(KIND=kdp) :: ar, sr, ur, cvr, dpdtr
  ! ... $Revision: 1.7 $//$Date: 2011/07/19 21:24:24 $

CONTAINS

  SUBROUTINE tempress(subt, subp, subd, subh)
    ! ... ********a poor subroutine name for what it does
    ! ... Purpose:  To calculate the enthaply and density of water or steam
    ! ...         as a function of temperature and pressure
    ! ... Method: This subroutine uses the program developed for the
    ! ...         Steam Tables.  The code in this subroutine is an
    ! ...         abreviated version of the one given by Haar et al., 1984.
    ! ... Reference:
    ! ...     Adapted from p. 302-312 of the NBS/NRC steam tables
    ! ...     by Haar, Gallagher, and Kell, and published by Hemisphere
    ! ...     Publishing Corp., 1984.
    ! ... Units used in this subroutine:
    ! ...     Temperature [K]
    ! ...     Pressure    [MPa]
    ! ...     Density     [g/cm^3]
    ! ...     Enthalpy    [kJ/kg]
    !
    ! ...   The values set in this module supply most of the fixed parameters
    ! ...  used in the rest of the routines. PX is the b(i) of table I,
    ! ...  Q is the b(i) of table I, G1,G2 and FG are the alpha, beta
    ! ...  and gamma of eq 2, and G,II,JJ are the g(i), k(i) and l(i) of eq 5.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: subt            !.. [Deg.C]
    REAL(KIND=kdp), INTENT(IN) :: subp            !.. [dyne/cm^2] 
    REAL(KIND=kdp), INTENT(OUT) :: subd, subh     !.. [g/cm^3]; [erg/g]
    !
    REAL(KIND=kdp) :: d, dgss, dll, dq, dvv, h, psat, px, pxx
    REAL(KIND=kdp) :: rt, t 
    REAL(KIND=kdp), PARAMETER :: tcrit=647.1260000001_kdp
    ! ... Set string for use with RCS ident command
    !CHARACTER(LEN=80) :: ident_string='$Revision: 1.7 $//$Date: 2011/07/19 21:24:24 $'
    !     ------------------------------------------------------------------
    ! ... special patch to handle low pressures; set to a minimum pressure
    ! ...      for the tables
    pxx = MAX(subp,5.e5_kdp)
    ! ...        convert input P&T to local units; MPa, K
    px = pxx*1.0e-7_kdp
    t = subt + 273.15_kdp
    rt = gascon*t
    CALL bb(t)
    ! ... Calculate enthalpy and density as a function of pressure and temperature
    dgss = px/(t*0.4_kdp)
    psat = 2.e4_kdp
    dll = 0._kdp
    dvv = 0._kdp
    IF (t < tcrit) CALL pcorr(t, psat, dll, dvv)
    IF (px > psat) dgss = dll
    CALL dfind(d, px, dgss, t, dq)
    CALL therm(d, t)        ! ... calculate dimensionless thermo params
    h = hd*rt
    ! ...        convert density and enthalpy to proper units for output
    subh = h*1.0e7_kdp
    subd = d
  END SUBROUTINE tempress

  SUBROUTINE base(d, t, zbase)
    ! ... This routine calculates z (=pbase/(drt) of eq. q) (called base),
    ! ... and also abase,gbase,sbase,ubase,hbase,cvbase, and 1/(drt) * dp/dt
    ! ... for the base fct (called dpdtb).
    ! ...  the ab,gb,sb,ub,hb and cvb are calculated in dimensionless units
    ! ...  ab/rt, gb/rt, sb/r, etc.
    ! ...  g1,g2 and gf are the alpha, beta and gamma of eq 2, which are
    ! ...  supplied in the initialization.  b1 and b2 are the "excluded
    ! ...  volume" and "2nd virial" (eqs 3 and 4) supplied by the subroutine
    ! ...  bb(t). which also supplies the 1st and 2nd derivatives with
    ! ...  respect to t (b1t,b2t,b1tt,b2tt).
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: d, t
    REAL(KIND=kdp), INTENT(OUT) :: zbase
    !
    REAL(KIND=kdp) :: bb2tt, dz0, x, z0
    !     ------------------------------------------------------------------
    !...
    y = 0.25_kdp*b1*d
    x = 1._kdp - y
    z0 = (1._kdp+g1*y+g2*y*y)/x**3
    z = z0 + 4._kdp*y*(b2/b1-gf)
    dz0 = (g1+2._kdp*g2*y)/x**3 + 3._kdp*(1._kdp+g1*y+g2*y*y)/x**4
    dzren = dz0 + 4._kdp*(b2/b1-gf)
    ab = -LOG(x) - (g2-1._kdp)/x + 28.16666667_kdp/(x*x) +  &
         4._kdp*y*(b2/b1-gf) + 15.166666667_kdp + LOG(d*t*gascon/0.101325_kdp)
    gb = ab + z
    zbase = z
    bb2tt = t*t*b2tt
    ub = -t*b1t*(z-1._kdp-d*b2)/b1 - d*t*b2t
    hb = z + ub
    cvb = 2._kdp*ub + (z0-1._kdp)*((t*b1t/b1)**2-t*t*b1tt/b1)  &
         - d*(bb2tt-gf*b1tt*t*t) - (t*b1t/b1)**2*y*dz0
    dpdtb = zbase/t + zbase*d/z*(dzren*b1t/4._kdp+b2t-b2/b1*b1t)
    sb = ub - ab
  END SUBROUTINE base

  SUBROUTINE bb(t)
    ! ...  This subroutine calculates the b's of eqs. 3,4 using coefficients
    ! ...  from initialization, calculating also the first and second derivative
    ! ...  with respect to tempeature. the b's calculated here are in cm^3/g.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: t
    !
    INTEGER :: i
    REAL(KIND=kdp), DIMENSION(10) :: v
    !     ------------------------------------------------------------------
    !...
    v(1) = 1._kdp
    DO  i = 2,10
       v(i) = v(i-1)*tzren/t
    END DO
    b1 = bp(1) + bp(2)*LOG(1._kdp/v(2))
    b2 = bq(1)
    b1t = bp(2)*v(2)/tzren
    b2t = 0._kdp
    b1tt = 0._kdp
    b2tt = 0._kdp
    DO  i = 3,10
       b1 = b1 + bp(i)*v(i-1)
       b2 = b2 + bq(i)*v(i-1)
       b1t = b1t - (i-2)*bp(i)*v(i-1)/t
       b2t = b2t - (i-2)*bq(i)*v(i-1)/t
       b1tt = b1tt + bp(i)*(i-2)**2*v(i-1)/t/t
       b2tt = b2tt + bq(i)*(i-2)**2*v(i-1)/t/t
    END DO
    b1tt = b1tt - b1t/t
    b2tt = b2tt - b2t/t
  END SUBROUTINE bb

  SUBROUTINE corr(t, px, dl, dv, delg)
    ! ...  subroutine corr will calculate, for an input t and px at or near
    ! ...  the  vapor pressure, the corresponding liquid and vapor densities
    ! ...  and also  delg = (gl-gv)/rt for use in calculating the correction
    ! ...  to the vapor pressure in pcorr or vapor temperature in tcorr for
    ! ...  delg = 0.  note that delg is defined as zero for all temperatures
    ! ...  above 646.3K.
    ! ... px and dl and dv do not have the right logical intents.
    !***  px is overwritten. Would be better to have an output px consistent 
    ! ***    with dl and dv.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: t
    REAL(KIND=kdp), INTENT(IN OUT) :: px
    REAL(KIND=kdp), INTENT(IN OUT) :: dl, dv
!!$    REAL(KIND=kdp), INTENT(IN OUT) :: delg
    REAL(KIND=kdp), INTENT(OUT) :: delg
    !
    REAL(KIND=kdp) :: dliq, dq, dvap, gl, gv, rt, tau, zdum
    REAL(KIND=kdp), PARAMETER :: pcrit=22.055_kdp, tcrit=647.1260000001_kdp
    !     ------------------------------------------------------------------
    !...
    rt = gascon*t
    CALL bb(t)
    IF (t <= 646.3_kdp) THEN
       dliq = dl
       IF (dl <= 0.0_kdp) dliq = 1.11_kdp - 0.0004_kdp*t
       CALL dfind(dl, px, dliq, t, dq)
       CALL therm(dl, t)
       gl = gd
       dvap = dv
       IF (dv <= 0.0_kdp) dvap = px/rt
       CALL dfind(dv, px, dvap, t, dq)
       IF (dv < 5.e-7_kdp) dv = 5.e-7_kdp
       CALL therm(dv, t)
       gv = gd
       delg = gl - gv
    ELSE
       px = pcrit
       delg = 0.0_kdp
       IF (t <= tcrit) THEN
          tau = 0.657128_kdp*(1.0_kdp-t/647.126_kdp)**0.325_kdp
          dl = 0.322_kdp + tau
          dv = 0.322_kdp - tau
          CALL base(dv, t, zdum)
!!$          CALL qq(t, dv)
          CALL qqsub(t, dv)     !!$*** for flint
          px = rt*dv*z + q0
       END IF
    END IF
  END SUBROUTINE corr

  SUBROUTINE dfind(dout, px, d, t, dpd)
    ! ...  Routine to find density corresponding to input pressure P(MPa), and
    ! ...  temperature T(K), using initial guess density D(g/cm^3). The output
    ! ...  density is in g/cm^3, also, as a byproduct, DP/DRHO is calculated
    ! ...  ("DPD", MPa cm^3/g)
!!$    USE machine_constants, ONLY: kdp
    USE f_units
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(OUT) :: dout, dpd
    REAL(KIND=kdp), INTENT(IN) :: px, d, t
    !
    INTEGER :: l
    REAL(KIND=kdp) :: dd, dp, dpdx, pp, rt, x, zdum
    !     ------------------------------------------------------------------
    !...
    dd = d
    rt = gascon*t
    IF (dd <= 0._kdp) dd = 1.e-8_kdp
    dd = MIN(dd,1.9_kdp) 
    l = 0
10  l = l + 1
    IF (dd <= 0._kdp) dd = 1.e-8_kdp
    dd = MIN(dd,1.9_kdp) 
!!$    CALL qq(t, dd)
    CALL qqsub(t, dd)     !!$*** for flint
    CALL base(dd, t, zdum)
    pp = rt*dd*zdum + q0
    dpd = rt*(z+y*dzren) + q5
    ! ...  check for negative DP/DRHO
    IF (dpd <= 0._kdp) THEN
       ! ...  assume guess to be in 2-phase region, adjust guess accordingly
       IF (d >= 0.2967_kdp) dd = dd*1.02_kdp
       IF (d < 0.2967_kdp) dd = dd*0.98_kdp
       IF (l <= 10) GO TO 10
    END IF
    dpdx = dpd*1.1_kdp
    dpdx = MAX(dpdx,0.1_kdp)
    dp = ABS(1._kdp-pp/px)
    IF (dp >= 1.e-8_kdp) THEN
       IF (d <= 0.3_kdp .OR. dp >= 1.e-7_kdp) THEN
          IF (d <= 0.7_kdp .OR. dp >= 1.e-6_kdp) THEN
             x = (px-pp)/dpdx
             IF (ABS(x) > 0.1_kdp) x = x*0.1_kdp/ABS(x)
             dd = dd + x
             IF (dd <= 0._kdp) dd = 1.e-8_kdp
             IF (l <= 30) GO TO 10
             WRITE (fustdout, 9005) px, t, d
             WRITE (fuclog, 9005) px, t, d
9005         FORMAT ('***Convergence failure in DFIND'/' PX =',1pg15.7,' T =',  &
                  1pg15.7,' D =',1pg15.7)
          END IF
       END IF
    END IF
    dout = dd
  END SUBROUTINE dfind

  SUBROUTINE ideal(t)
    ! ...  This subroutine calculates the thermodynamic properties for
    ! ...  water in the ideal gas state from function of H.W. Woolley
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: t
    !
    INTEGER :: i
    REAL(KIND=kdp) :: tl, tt
    REAL(KIND=kdp), DIMENSION(18) :: c=(/0.19730271018e2_kdp, 0.209662681977e2_kdp,  &
         -0.483429455355_kdp, 0.605743189245e1_kdp, 22.56023885_kdp, -9.87532442_kdp,  &
         -0.43135538513e1_kdp, 0.458155781_kdp, -0.47754901883e-1_kdp,  &
         0.41238460633e-2_kdp, -0.27929052852e-3_kdp, 0.14481695261e-4_kdp,  &
         -0.56473658748e-6_kdp, 0.16200446e-7_kdp, -0.3303822796e-9_kdp,  &
         0.451916067368e-11_kdp, -0.370734122708e-13_kdp, 0.137546068238e-15_kdp/)
    !     ------------------------------------------------------------------
    !...
    tt = t/100._kdp
    tl = LOG(tt)
    gi = -(c(1)/tt+c(2))*tl
    hi = (c(2)+c(1)*(1._kdp-tl)/tt)
    cpi = c(2) - c(1)/tt
    DO  i = 3,18
       gi = gi - c(i)*tt**(i-6)
       hi = hi + c(i)*(i-6)*tt**(i-6)
       cpi = cpi + c(i)*(i-6)*(i-5)*tt**(i-6)
    END DO
    ai = gi - 1._kdp
    ui = hi - 1._kdp
    cvi = cpi - 1._kdp
    si = ui - ai
  END SUBROUTINE ideal

  SUBROUTINE pcorr(t, px, dl, dv)
    ! ...  Calculates the vapor pressure p and the liquid
    ! ...  and vapor densities corresponding to the input t, corrected such
    ! ...  that gl-gv=0.  the function ps is required to give a
    ! ...  reasonably good approximation to the vapor pressure to be used as
    ! ...  the starting point for the iteration.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: t
    REAL(KIND=kdp), INTENT(OUT) :: px
    REAL(KIND=kdp), INTENT(INOUT) :: dl, dv       ! ... because initial guesses are input
    !
    REAL(KIND=kdp) :: delg, dp
    !     ------------------------------------------------------------------
    !...
    px = ps(t)
10  CALL corr(t, px, dl, dv, delg)
    dp = delg*gascon*t/((1._kdp/dv-1._kdp/dl)+1.0e-13_kdp)
    px = px + dp
    IF (ABS(delg) > 1.e-8_kdp) GOTO 10     ! ... no limit on the number of iterations
  END SUBROUTINE pcorr

  FUNCTION ps(t) RESULT(app_vp)
    ! ...  calculates an approximation to the vapor pressure,
    ! ...  app_vp, as a function of input temperature.  the vapor pressure
    ! ...  calculated agrees with the vapor pressure predicted by the surface
    ! ...  to within 0.02% to within a degree or so of the critical temperature
    ! ...  and can serve as an initial guess for further refinement by
    ! ...  imposing the condition that g1=gv.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp) :: app_vp
    REAL(KIND=kdp), INTENT(IN) :: t
    !
    INTEGER :: i
    REAL(KIND=kdp) :: b, pl, q, v, w, z
    REAL(KIND=kdp), DIMENSION(8) :: a=(/-7.8889166_kdp, 2.5514255_kdp, -6.716169_kdp,  &
         33.239495_kdp, -105.38479_kdp, 174.35319_kdp, -148.39348_kdp, 48.631602_kdp/)
    !     ------------------------------------------------------------------
    !...
    IF (t <= 314._kdp) THEN
       pl = 6.3573118_kdp - 8858.843_kdp/t + 607.56335_kdp*t**(-0.6_kdp)
       app_vp = 0.1_kdp*EXP(pl)
    ELSE IF (t <= 646.3_kdp) THEN
       v = t/647.25_kdp
       w = ABS(1._kdp-v)
       b = 0._kdp
       DO  i = 1,8
          z = i
          b = b + a(i)*w**((z+1._kdp)/2._kdp)
       END DO
       q = b/v
       app_vp = 22.093_kdp*EXP(q)
    ELSE
       app_vp = -147.7854939_kdp + 0.262453402_kdp*t
    END IF
  END FUNCTION ps

!!$  SUBROUTINE qq(t, d)
  SUBROUTINE qqsub(t, d)      !!$*** for flint
    ! ... Calculates, for input t(K) and d(g/cm^3), the
    ! ...  residual contributions to pressure (q), Helmholtz function (ar),
    ! ...  dp/drhp (q5), and also the the Gibbs function, entropy,internal
    ! ...  energy, enthalpy, isochoric heat capacity and dpdt. (eq 5).
    ! ...  terms 37 thru 39 are the additional terms affecting only the
    ! ...  immediate vicinity of the critical point, and term 40 is the
    ! ...  additional term to improve results in the low t, high p region.
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: t, d
    !
    INTEGER :: i, j, k, km, l
    REAL(KIND=kdp) :: att, d2f, dadt, dd, ddz, del, dex, dfdt, dpt, e, ex1, ex2, &
         fct, q10, q20, q2a, q5t, qm, qp, rt, tau, tex, tx, v, zz
    REAL(KIND=kdp), DIMENSION(11), TARGET :: qr
    REAL(KIND=kdp), DIMENSION(10), TARGET :: qt
    REAL(KIND=kdp), DIMENSION(:), POINTER :: qzr, qzt
    !     ------------------------------------------------------------------
    !...
    qzr => qr(3:11)
    qzt => qt(2:10)
    rt = gascon*t
    qr(1) = 0._kdp
    q5 = 0._kdp
    q0 = 0._kdp
    ar = 0._kdp
    dadt = 0._kdp
    cvr = 0._kdp
    dpdtr = 0._kdp
    e = EXP(-aa*d)
    q10 = d*d*e
    q20 = 1._kdp - e
    qr(2) = q10
    v = tzren/t
    qt(1) = t/tzren
    DO  i=2,10
       qr(i+1) = qr(i)*q20
       qt(i) = qt(i-1)*v
    END DO
    DO  i=1,nc
       k = ii(i) + 1
       l = jj(i)
       zz = k
!!$     QP=Gren(I)*AA*QZR(K-1)*QZT(L)
       qp = gren(i)*aa*qr(k+1)*qzt(l)     ! ... correction of previous line because qzr(k-1)
                                          ! ...     gives zero subscript error when K=1
       q0 = q0 + qp
       q5 = q5 + aa*(2._kdp/d-aa*(1._kdp-e*(k-1)/q20))*qp
       ar = ar + gren(i)*qzr(k)*qzt(l)/(q10*zz*rt)
       dfdt = (q20**k)*(1-l)*qzt(l+1)/(tzren*k)
       d2f = l*dfdt
       dpt = dfdt*q10*aa*k/q20
       dadt = dadt + gren(i)*dfdt
       dpdtr = dpdtr + gren(i)*dpt
       cvr = cvr + gren(i)*d2f/gascon
    END DO
    qp = 0._kdp
    q2a = 0._kdp
    DO  j=37,40
       IF (gren(j) == 0._kdp) CYCLE
       k = ii(j)
       km = jj(j)
       ddz = adz(j-36)
       del = d/ddz - 1._kdp
       IF(ABS(del) < 1.e-10_kdp) del = 1.e-10_kdp
!!$       dd = del*del
       ex1 = -aad(j-36)*del**k
       dex = EXP(ex1)*del**km
       att = aat(j-36)
       tx = atz(j-36)
       tau = t/tx - 1._kdp
       ex2 = -att*tau*tau
       tex = EXP(ex2)
       q10 = dex*tex
       qm = km/del - k*aad(j-36)*del**(k-1)
       fct = qm*(d**2)*q10/ddz
       q5t = fct*(2._kdp/d+qm/ddz) - ((d/ddz)**2)*q10*  &
            (km/del/del+k*(k-1)*aad(j-36)*(del**(k-2)))
       q5 = q5 + q5t*gren(j)
       qp = qp + gren(j)*fct
       dadt = dadt - 2._kdp*gren(j)*att*tau*q10/tx
       dpdtr = dpdtr - 2._kdp*gren(j)*att*tau*fct/tx
       q2a = q2a + t*gren(j)*(4._kdp*att*ex2+2._kdp*att)*q10/(tx*tx)
       ar = ar + q10*gren(j)/rt
    END DO
    sr = -dadt/gascon
    ur = ar + sr
    cvr = cvr + q2a/gascon
    q0 = q0 + qp
!!$  END SUBROUTINE qq
  END SUBROUTINE qqsub     !!$*** for flint

  SUBROUTINE therm(d, t)
    ! ... Calculates the thermodynamic functions in
    ! ...  dimensionless units (ad=a-rt, gd=g/rt, sd=s/r, ud=u/rt,
    ! ...  hd=h/rt, cvd=cv/r, and cpd=cp/r)
!!$    USE machine_constants, ONLY: kdp
    IMPLICIT NONE
    REAL(KIND=kdp), INTENT(IN) :: d, t
    !
    REAL(KIND=kdp) :: dpdd, rt, zz
    !     ------------------------------------------------------------------
    !...
    CALL ideal(t)
    rt = gascon*t
    zz = z + q0/(rt*d)
    dpdd = rt*(z+y*dzren) + q5
    ad = ab + ar + ai - uref/t + sref
    gd = ad + zz
    ud = ub + ur + ui - uref/t
    dpdt = rt*d*dpdtb + dpdtr
    cvd = cvb + cvr + cvi
    cpd = cvd + t*(dpdt**2)/(d*d*dpdd*gascon)
    hd = ud + zz
    sd = sb + sr + si - sref
    dvdt = dpdt/(dpdd*d*d)
    cjtt = 1._kdp/d - t*dvdt
    cjth = -cjtt/(cpd*gascon)
  END SUBROUTINE therm

END MODULE steam_table_funct

