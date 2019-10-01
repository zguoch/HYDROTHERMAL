SUBROUTINE getunits(array, dumstring, kodd)
  ! ... Purpose:  To determine the units of measure used for input for each of the
  ! ...           array variables
  USE machine_constants, ONLY: kdp
  USE f_units
  USE control
  USE units
  IMPLICIT NONE
  CHARACTER(LEN=4), INTENT(IN) :: array
  CHARACTER(LEN=76), INTENT(IN) :: dumstring
  INTEGER, INTENT(OUT) :: kodd
  ! ...   KODD -  array i.d. number
  ! ...        1= Porosity
  ! ...        2= X - Permeability
  ! ...        3= Y - Permeability
  ! ...        4= Z - Permeability
  ! ...        5= Thermal Conductivity
  ! ...        6= Specific Heat of Rock
  ! ...        7= Rock Density
  ! ...        8= Rock compressibility
  ! ...        9= Conductive Heat Flux at Base of region
  ! ...       10= Enthalpy
  ! ...       11= Pressure
  ! ...       12= X - spacing
  ! ...       13= Y - spacing
  ! ...       14= Z - spacing
  ! ...       15= Temperature
  ! ...       16= Source/sink 1 - mass input rate
  ! ...       17= Source/sink 2 - enthalpy or temperature of source
  ! ...       18= Precipitation flux at top of region
  ! ...       19= Associated Temperature for Precipitation flux at top of region
  !
!!$  CHARACTER(LEN=20), EXTERNAL :: uppercase
  CHARACTER(LEN=20) :: uunits
  CHARACTER(LEN=13) :: unitsavail(6), ht_cgs
  INTEGER :: ldum1, ldum2, lll
  LOGICAL :: lproblem
  INTERFACE
     FUNCTION uppercase(string) RESULT(outstring)
       CHARACTER(LEN=*), INTENT(IN) :: string
       CHARACTER(LEN=LEN(string)) :: outstring
     END FUNCTION uppercase
  END INTERFACE
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 1.13 $//$Date: 2007/12/10 21:57:45 $'
  !     ------------------------------------------------------------------
   IF (informt >= 3) THEN
     CALL getunits_3
   ELSE
     CALL getunits_2
   END IF
   RETURN
CONTAINS
   SUBROUTINE getunits_3
!**** Units allowed for version 3 and onward
     IMPLICIT NONE
     ! ... initialize
     kodd = 0
     unitsavail = '             '
     ! ... extract the portion of input line between parentheses containing units
     ldum1 = INDEX(dumstring,'(')
     ldum2 = INDEX(dumstring,')')
     ! ... Check to see that input line contains parentheses
     IF (ldum1 == 0 .OR. ldum2 == 0 .OR. ldum2 == ldum1+1) THEN
        lproblem = .TRUE.
        uunits = '???'
     ELSE
        lproblem = .FALSE.
        ldum1 = ldum1 + 1
        ldum2 = ldum2 - 1
        uunits = dumstring(ldum1:ldum2)
        WRITE(fupdef,'(tr13,2a)') 'User specified input units: ', uunits
        uunits = uppercase(uunits)
     END IF
     ! ... Porosity
     IF (array == 'PORO') THEN
        kodd = 1
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(-)'
        unitsavail(2) = '(%)'
        ht_cgs = unitsavail(1)
        ! ...        cgs units: unitless
        IF (uunits == 'UNITLESS' .OR. uunits == 'NONE' .OR. uunits == '-') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(-)'
           ! ...           units: %
        ELSE IF (uunits == '%' .OR. uunits == 'PERCENT') THEN
           unitfac(kodd) = 1.0e-2_kdp
           unitlabel(kodd) = '(%)'
        END IF
        ! ... X, Y, or Z Permeability
     ELSE IF (array == 'XPER' .OR. array == 'YPER' .OR. array == 'ZPER') THEN
        IF (array == 'XPER') kodd = 2
        IF (array == 'YPER') kodd = 3
        IF (array == 'ZPER') kodd = 4
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(cm^2)'
        unitsavail(2) = '(m^2)'
        unitsavail(3) = '(darcy)'
        ht_cgs = unitsavail(1)
        ! ...      cgs units: cm2
        IF (uunits == 'CM2' .OR. uunits == 'CM^2' .OR. uunits == 'CM*2') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(cm^2)'
           ! ...           units: m2
        ELSE IF (uunits == 'M2' .OR. uunits == 'M^2' .OR. uunits == 'M*2') THEN
           unitfac(kodd) = 1.0e4_kdp
           unitlabel(kodd) = '(m^2)'
           ! ...           units: darcy
        ELSE IF (uunits == 'DARCY') THEN
           unitfac(kodd) = 9.87e-9_kdp
           unitlabel(kodd) = '(darcy)'
        END IF
        ! ... Thermal Conductivity
     ELSE IF (array == 'THER') THEN
        kodd = 5
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(erg/s-cm-K)'
        unitsavail(2) = '(erg/s-cm-C)'
        unitsavail(3) = '(J/s-m-K)'
        unitsavail(4) = '(J/s-m-C)'
        unitsavail(5) = '(W/m-K)'
        unitsavail(6) = '(W/m-C)'
        ht_cgs = unitsavail(1)
        ! ...          cgs units: erg/s-cm-K
        IF (uunits == 'ERG/CM SEC K' .OR. uunits == 'ERG/SEC-CM-K' .OR.  &
             uunits == 'ERG/CM SEC C' .OR. uunits == 'ERG/SEC-CM-C' .OR.  &
             uunits == 'ERG/CM S K' .OR. uunits == 'ERG/S-CM-K' .OR.  &
             uunits == 'ERG/CM S C' .OR. uunits == 'ERG/S-CM-C') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(erg/s-cm-K)'
           ! ...           units: W/m K
        ELSE IF (uunits == 'W/M K' .OR. uunits == 'W/M-K' .OR.  & 
             uunits == 'W/M C' .OR. uunits == 'W/M-C') THEN 
           unitfac(kodd) = 1.0e5_kdp
           unitlabel(kodd) = '(W/m-K)'
           ! ...           units: J/s m K
        ELSE IF (uunits == 'J/SEC M K' .OR. uunits == 'J/SEC-M-K' .OR.  &
             uunits == 'J/SEC M C' .OR. uunits == 'J/SEC-M-C' .OR.  &
             uunits == 'J/S M K' .OR. uunits == 'J/S-M-K' .OR.  &
             uunits == 'J/S M C' .OR. uunits == 'J/S-M-C') THEN
           unitfac(kodd) = 1.0e5_kdp
           unitlabel(kodd) = '(J/s-m-K)'
        END IF
        ! ... Specific Heat (Heat Capacity)
     ELSE IF (array == 'SPEC') THEN
        kodd = 6
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(erg/g-K)'
        unitsavail(2) = '(J/g-K)'
        unitsavail(3) = '(J/kg-K)'
        unitsavail(4) = '(kJ/kg-K)'
        unitsavail(5) = '(erg/g-C)'
        unitsavail(6) = '(J/g-C)'
        ht_cgs = unitsavail(1)
        ! ...          cgs units: erg/g-K
        IF (uunits == 'ERG/G K' .OR. uunits == 'ERG/G-K' .OR.  &
             uunits == 'ERG/G C' .OR. uunits == 'ERG/G-C') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(erg/g-K)'
           ! ...           units: J/g-K
        ELSE IF (uunits == 'J/G K' .OR. uunits == 'J/G-K' .OR.  &
             uunits == 'J/G C' .OR. uunits == 'J/G-C') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(J/g-K)'
           ! ...           units: J/kg-K
        ELSE IF (uunits == 'J/KG K' .OR. uunits == 'J/KG-K' .OR.  & 
             uunits == 'J/KG C' .OR. uunits == 'J/KG-C') THEN
           unitfac(kodd) = 1.0e4_kdp
           unitlabel(kodd) = '(J/kg-K)'
           ! ...           units: kJ/kg-K
        ELSE IF (uunits == 'KJ/KG K' .OR. uunits == 'KJ/KG-K' .OR.  &
             uunits == 'KJ/KG C' .OR. uunits == 'KJ/KG-C') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(kJ/kg-K)'
        END IF
        ! ... Rock Density
     ELSE IF (array == 'RDEN') THEN
        kodd = 7
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(g/cm^3)'
        unitsavail(2) = '(kg/m^3)'
        ht_cgs = unitsavail(1)
        ! ...             cgs units: g/cm^3
        IF (uunits == 'G/CM3' .OR. uunits == 'G/CM^3') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(g/cm^3)'
           ! ...           units: kg/m3
        ELSE IF (uunits == 'KG/M3' .OR. uunits == 'KG/M^3') THEN
           unitfac(kodd) = 1.0e-3_kdp
           unitlabel(kodd) = '(kg/m^3)'
        END IF
        ! ... Porous Matrix Compressibility
     ELSE IF (array == 'COMP') THEN
        kodd = 8
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(cm^2/dyne)'
        unitsavail(2) = '(bar^-1)'
        unitsavail(3) = '(Pa^-1)'
        unitsavail(4) = '(MPa^-1)'
        ht_cgs = unitsavail(1)
        ! ...             cgs units: cm^2/dyne
        IF (uunits == 'CM2/DYNE' .OR. uunits == 'CM^2/DYNE') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(cm^2/dyne)'
           ! ...          units: bar^-1
        ELSE IF (uunits == '1/BARS' .OR. uunits == '1/BAR' .OR.  &
             uunits == 'BARS^-1' .OR. uunits == 'BAR^-1' .OR.  &
             uunits == 'BARS-1' .OR. uunits == 'BAR-1') THEN
           unitfac(kodd) = 1.0e-6_kdp
           unitlabel(kodd) = '(bar^-1)'
           ! ...          units: Pa^-1
        ELSE IF (uunits == '1/PA' .OR. uunits == 'PA^-1' .OR.  &
             uunits == 'PA-1') THEN
           unitfac(kodd) = 1.0e-1_kdp
           unitlabel(kodd) = '(Pa^-1)'
           ! ...          units: MPa^-1
        ELSE IF (uunits == '1/MPA' .OR. uunits == 'MPA^-1' .OR.  &
             uunits == 'MPA-1') THEN
           unitfac(kodd) = 1.0e-7_kdp
           unitlabel(kodd) = '(MPa^-1)'
        END IF
        ! ... Conductive Heat Flux along the Base
     ELSE IF (array == 'COND') THEN
        kodd = 9
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(erg/s-cm^2)'
        unitsavail(2) = '(mW/m^2)'
        unitsavail(3) = '(W/m^2)'
        unitsavail(4) = '(J/s-m^2)'
        ht_cgs = unitsavail(1)
        ! ...      cgs units: erg/s-cm^2
        IF (uunits == 'ERGS/SEC-CM^2' .OR.  &
             uunits == 'ERGS/SEC-CM2' .OR. uunits == 'ERGS/S-CM^2' .OR.  &
             uunits == 'ERGS/S-CM2' .OR. uunits == 'ERG/SEC-CM^2' .OR.  &
             uunits == 'ERG/SEC-CM2' .OR. uunits == 'ERG/S-CM^2' .OR.  &
             uunits == 'ERG/S-CM2' .OR. uunits == 'ERGS/SEC CM^2' .OR.  &
             uunits == 'ERGS/SEC CM2' .OR. uunits == 'ERGS/S CM^2' .OR.  &
             uunits == 'ERGS/S CM2' .OR. uunits == 'ERG/SEC CM^2' .OR.  &
             uunits == 'ERG/SEC CM2' .OR. uunits == 'ERG/S CM^2' .OR.  &
             uunits == 'ERG/S CM2') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(erg/s-cm^2)'
           ! ...          units:  mW/m^2
        ELSE IF (uunits == 'MW/M^2' .OR. uunits == 'MW/M2') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(mW/m^2)'
           ! ...          units:  W/m^2
        ELSE IF (uunits == 'W/M^2' .OR. uunits == 'W/M2') THEN 
           unitfac(kodd) = 1.0e3_kdp
           unitlabel(kodd) = '(W/m^2)'
           ! ...          units: J/s-m^2
        ELSE IF (uunits == 'J/SEC-M^2' .OR. uunits == 'J/SEC-M2' .OR.  &
             uunits == 'J/S-M^2' .OR. uunits == 'J/S-M2' .OR.  &
             uunits == 'J/SEC M^2' .OR. uunits == 'J/SEC M2' .OR.  &
             uunits == 'J/S M^2' .OR. uunits == 'J/S M2') THEN
           unitfac(kodd) = 1.0e3_kdp
           unitlabel(kodd) = '(J/s-m^2)'
        END IF
        ! ... Enthalpy
     ELSE IF (array == 'ENTH') THEN
        kodd = 10
        ! ...      set default units for temperature when enthalpy data input
        unitfac(15) = 1._kdp
        unitlabel(15) = '(C)'
        WRITE(fupdef,'(tr13,a,a)')  &
             'Default output units assigned for Temperature: ',unitlabel(15)
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(erg/g)'
        unitsavail(2) = '(J/g)'
        unitsavail(3) = '(kJ/kg)'
        unitsavail(4) = '(J/kg)'
        ht_cgs = unitsavail(1)
        ! ...         cgs units: erg/g
        IF (uunits == 'ERGS/G' .OR. uunits == 'ERG/G') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(erg/g)'
           ! ...          units: J/g
        ELSE IF (uunits == 'J/G') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(J/g)'
           ! ...          units: kJ/kg
        ELSE IF (uunits == 'KJ/KG') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(kJ/kg)'
           ! ...          units: J/kg
        ELSE IF (uunits == 'J/KG') THEN
           unitfac(kodd) = 1.0e4_kdp
           unitlabel(kodd) = '(J/kg)'
        END IF
        ! ... Pressure
     ELSE IF (array == 'PRES') THEN
        kodd = 11
        unitfac(kodd) = 0._kdp  
        ! ...          available input units
        unitsavail(1) = '(dyne/cm^2)'
        unitsavail(2) = '(MPa)'
        unitsavail(3) = '(bar)'
        unitsavail(4) = '(Pa)'
        ht_cgs = unitsavail(1)
        !..        unitsavail(5) = '(atm)'
        ! ...         cgs units: dyne/cm^2
        IF (uunits == 'DYNES/CM2' .OR. uunits == 'DYNES/CM^2' .OR.  &
             uunits == 'DYNE/CM2' .OR. uunits == 'DYNE/CM^2') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(dyne/cm^2)'
           ! ...           units: MPa
        ELSE IF (uunits == 'MPA') THEN
           unitfac(kodd) =  1.e7_kdp
           unitlabel(kodd) = '(MPa)'
           ! ...           units: bars
        ELSE IF (uunits == 'BARS' .OR. uunits == 'BAR') THEN
           unitfac(kodd) = 1.0e6_kdp
           unitlabel(kodd) = '(bar)'
           ! ...           units: Pa
        ELSE IF (uunits == 'PA') THEN
           unitfac(kodd) = 10._kdp
           unitlabel(kodd) = '(Pa)'
           ! ...           units: atm
           ! ...    ELSEIF (uunits.EQ.'atm') THEN
           ! ...      Unitfac(KODD) = 1.01325D6
           ! ...      Unitlabel(KODD) = '(atm)'
        END IF
        ! ... X, Y, or Z  Spacing
     ELSE IF (array == 'XSPA' .OR. array == 'YSPA' .OR. array == 'ZSPA') THEN
        IF (array == 'XSPA') kodd = 12
        IF (array == 'YSPA') kodd = 13
        IF (array == 'ZSPA') kodd = 14
        ! ...          available input units
        unitsavail(1) = '(cm)'
        unitsavail(2) = '(m)'
        unitsavail(3) = '(km)'
        ht_cgs = unitsavail(1)
        ! ...           cgs units: cm
        IF (uunits == 'CM') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(cm)'
           ! ...           units: m
        ELSE IF (uunits == 'M') THEN
           unitfac(kodd) = 1.0e2_kdp
           unitlabel(kodd) = '(m)'
           ! ...           units: km
        ELSE IF (uunits == 'KM') THEN
           unitfac(kodd) = 1.0e5_kdp
           unitlabel(kodd) = '(km)'
        END IF
        ! ... Temperature
     ELSE IF (array == 'TEMP') THEN
        kodd = 15
        ! ...      set default units for enthalpy when temperature data input
        unitfac(10) = 1.0e7_kdp
        unitlabel(10) = '(kJ/kg)'
        WRITE(fupdef,'(tr13,a,a)')  &
             'Default output units assigned for Enthalpy: ',unitlabel(10)
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(Deg.C)'
        ht_cgs = unitsavail(1)
        ! ...           cgs units: C
        IF (uunits == 'C' .OR. uunits == 'DEG C' .OR.  uunits == 'DEG.C' .OR.  &
             uunits == 'DEG-C' .OR. uunits == 'DEGREES CELSIUS' .OR.  &
             uunits == 'DEGREES-CELSIUS' .OR. uunits == 'DEGREE C' .OR.  &
             uunits == 'DEGREE-C' .OR. uunits == 'DEGREES C' .OR.  &
             uunits == 'DEGREES-C') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(Deg.C)'
        END IF
        ! ... Source/sink input 1 - Mass input rate
     ELSE IF (array == 'SRC1') THEN
        kodd = 16
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(g/s)'
        unitsavail(2) = '(kg/s)'
        unitsavail(3) = '(g/yr)'
        unitsavail(4) = '(kg/yr)'
        ht_cgs = unitsavail(1)
        ! ...         cgs units: g/s
        IF (uunits == 'G/S' .OR. uunits == 'G/SEC') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(g/s)'
           ! ...          units: kg/s
        ELSE IF (uunits == 'KG/S' .OR. uunits == 'KG/SEC') THEN
           unitfac(kodd) = 1.0e3_kdp
           unitlabel(kodd) = '(kg/s)'
           ! ...          units: g/yr
        ELSE IF (uunits == 'G/YR' .OR. uunits == 'G/YEAR') THEN
           unitfac(kodd) = 3.17e-8_kdp
           unitlabel(kodd) = '(g/yr)'
           ! ...          units: kg/yr
        ELSE IF (uunits == 'KG/YR' .OR. uunits == 'KG/YEAR') THEN
           unitfac(kodd) = 3.17e-5_kdp
           unitlabel(kodd) = '(kg/yr)'
        END IF
        ! ... Source/sink input 2 - Enthalpy or Temperature of Source
     ELSE IF (array == 'SRC2') THEN
        kodd = 17
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(erg/g)'
        unitsavail(2) = '(J/g)'
        unitsavail(3) = '(kJ/kg)'
        unitsavail(4) = '(J/kg)'
        unitsavail(5) = '(Deg.C)'
        ht_cgs = unitsavail(1)
        ! ...         cgs units: erg/g
        IF (uunits == 'ERGS/G' .OR. uunits == 'ERG/G') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(erg/g)'
           ! ...          units: J/g
        ELSE IF (uunits == 'J/G') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(J/g)'
           ! ...          units: kJ/kg
        ELSE IF (uunits == 'KJ/KG') THEN
           unitfac(kodd) = 1.0e7_kdp
           unitlabel(kodd) = '(kJ/kg)'
           ! ...          units: J/kg
        ELSE IF (uunits == 'J/KG') THEN
           unitfac(kodd) = 1.0e4_kdp
           unitlabel(kodd) = '(J/kg)'
           ! ...          units: C
        ELSE IF (uunits == 'C' .OR. uunits == 'DEG C' .OR. uunits == 'DEG.C' .OR.  &
             uunits == 'DEG-C' .OR. uunits == 'DEGREES CELSIUS' .OR.  &
             uunits == 'DEGREES-CELSIUS' .OR.  &
             uunits == 'DEGREE C' .OR. uunits == 'DEGREE-C' .OR.  &
             uunits == 'DEGREES C' .OR. uunits == 'DEGREES-C') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(Deg.C)'
           ht_cgs = '(Deg.C)'
        END IF
        ! ... Precipitation flux at top of region
     ELSE IF (array == 'PREC') THEN
        kodd = 18
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(cm^3/s-cm^2)'
        unitsavail(2) = '(cm/s)'
        unitsavail(3) = '(m^3/s-m^2)'
        unitsavail(4) = '(m/s)'
        unitsavail(5) = '(mm/yr)'
        unitsavail(6) = '(m/yr)'
        ht_cgs = unitsavail(1)
        ! ...      cgs units: cm^3/s-cm^2
        IF (uunits == 'CM^3/SEC-CM^2' .OR.  &
             uunits == 'CM3/SEC-CM2' .OR. uunits == 'CM^3/S-CM2' .OR.  &
             uunits == 'CM3/S-CM^2' .OR. uunits == 'CM^3/SEC-CM^2' .OR.  &
             uunits == 'CM3/SEC-CM2' .OR. uunits == 'CM^3/S-CM^2' .OR.  &
             uunits == 'CM3/S-CM2' .OR. uunits == 'CM^3/SEC CM^2' .OR.  &
             uunits == 'CM3/SEC CM2' .OR. uunits == 'CM^3/S CM^2' .OR.  &
             uunits == 'CM^3/S CM2' .OR. uunits == 'CM3/SEC CM^2' .OR.  &
             uunits == 'CM3/SEC CM2' .OR. uunits == 'CM^3/S CM^2' .OR.  &
             uunits == 'CM3/S CM2') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(cm^3/s-cm^2)'
           ! ...          units:  cm/s
        ELSE IF (uunits == 'CM/S' .OR. uunits == 'CM/SEC') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(cm/s)'
           ! ...          units: m^3/s-m^2
        ELSE IF (uunits == 'M^3/SEC-M^2' .OR.  &
             uunits == 'M3/SEC-M2' .OR. uunits == 'M^3/S-M2' .OR.  &
             uunits == 'M3/S-M^2' .OR. uunits == 'M^3/SEC-M^2' .OR.  &
             uunits == 'M3/SEC-M2' .OR. uunits == 'M^3/S-M^2' .OR.  &
             uunits == 'M3/S-M2' .OR. uunits == 'M^3/SEC M^2' .OR.  &
             uunits == 'M3/SEC M2' .OR. uunits == 'M^3/S M^2' .OR.  &
             uunits == 'M^3/S M2' .OR. uunits == 'M3/SEC M^2' .OR.  &
             uunits == 'M3/SEC M2' .OR. uunits == 'M^3/S M^2' .OR.  &
             uunits == 'M3/S M2') THEN
           unitfac(kodd) = 1.0e2_kdp
           unitlabel(kodd) = '(m^3/s-m^2)'
           ! ...          units:  m/s
        ELSE IF (uunits == 'M/S' .OR. uunits == 'M/SEC') THEN
           unitfac(kodd) = 1.0e2_kdp
           unitlabel(kodd) = '(m/s)'
           ! ...          units:  mm/yr
        ELSE IF (uunits == 'MM/YR') THEN
           unitfac(kodd) = 3.1688e-9_kdp
           unitlabel(kodd) = '(mm/yr)'
           ! ...          units:  m/yr
        ELSE IF (uunits == 'M/YR') THEN
           unitfac(kodd) = 3.1688e-6_kdp
           unitlabel(kodd) = '(m/yr)'
        END IF
        ! ... Associated Temperature for precipitation flux
     ELSE IF (array == 'TFLU') THEN
        kodd = 19
        unitfac(kodd) = 0._kdp
        ! ...          available input units
        unitsavail(1) = '(Deg.C)'
        ht_cgs = unitsavail(1)
        ! ...           cgs units: Deg.C
        IF (uunits == 'C' .OR. uunits == 'DEG C' .OR. uunits == 'DEG.C' .OR.  &
             uunits == 'DEG-C' .OR. uunits == 'DEGREES CELSIUS' .OR.  &
             uunits == 'DEGREES-CELSIUS' .OR. uunits == 'DEGREE C' .OR.  &
             uunits == 'DEGREE-C' .OR. uunits == 'DEGREES C' .OR.  &
             uunits == 'DEGREES-C') THEN
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = '(Deg.C)'
        END IF
     END IF               ! ... end of array types
     ! ... check for errors
     IF (lproblem) THEN       ! ... array line did not contain units enclosed in parentheses
        IF (informt == 1) THEN          ! ... old input format, write warning and use cgs units
           WRITE (fustdout, 9005) 'WARNING - Problem', array
           WRITE (fustdout, 9025) unitsavail(1)
           WRITE (fustdout, 9035) (unitsavail(lll), lll=2, 6)
           WRITE (fustdout, 9010)
           WRITE (fuclog, 9005) 'WARNING - Problem', array
   9005    FORMAT (/, 5X, ' **** ', a, ' determining units for array ', a,  &
                ' *****')
           WRITE (fuclog, 9025) unitsavail(1)
   9025    FORMAT (13X, 'Using HYDROTHERM cgs units ', a)
           WRITE (fuclog, 9035) (unitsavail(lll), lll=2, 6)
           WRITE (fuclog, 9010)
   9010    FORMAT (13X, 'Check that units are enclosed in parentheses and',  &
                /, 13X, 'that they are on the same line as the array name')
           unitfac(kodd) = 1._kdp
           unitlabel(kodd) = unitsavail(1)
           ierr(11) = .TRUE.
        ELSE          ! ... version 2 or version 3 input format
           WRITE (fustdout, 9005) '**** ERROR', array
           WRITE (fustdout, 9010)
           WRITE (fustdout, 9030) array, unitsavail(1)
           WRITE (fustdout, 9035) (unitsavail(lll), lll=2, 6)
           WRITE (fuclog, 9005) '**** ERROR', array
           WRITE (fuclog, 9010)
           WRITE (fuclog, 9030) array, unitsavail(1)
           WRITE (fuclog, 9035) (unitsavail(lll), lll=2, 6)
           ierr(12) = .TRUE.
           RETURN
        END IF
     END IF
     IF (kodd == 0) THEN
        WRITE (fustdout, 9015) array
        WRITE (fuclog, 9015) array
   9015 FORMAT (/, ' **** Problem with input array named: ', a,  &
             ' *****', /, 6X, 'Check array name')
        ierr(13) = .TRUE.
        RETURN
     END IF
     IF (unitfac(kodd) <= 0._kdp) THEN
        WRITE (fustdout, 9020) array, uunits
        WRITE (fustdout, 9030) array, unitsavail(1)
        WRITE (fustdout, 9035) (unitsavail(lll), lll=2, 6)
        WRITE (fuclog, 9020) array, uunits
   9020 FORMAT (/, ' **** ERROR determining units for array ', a,  &
             ' *****', /, 6X, 'Input units - ', a, ' not appropriate')
        WRITE (fuclog, 9030) array, unitsavail(1)
   9030 FORMAT (13X, 'HYDROTHERM cgs units for array ', a, ' are ', a)
        WRITE (fuclog, 9035) (unitsavail(lll), lll=2, 6)
   9035 FORMAT (13X, 'Other available input units are:', /, 11X, 5(2X,a))
        ierr(14) = .TRUE.
        RETURN
     END IF
     ! ...     if cgs units are used, check that conversion factor =1
     IF (unitlabel(kodd) == ht_cgs) THEN
        IF (ABS(unitfac(kodd)-1._kdp) > 1.0e-8_kdp) THEN
           WRITE (fustdout, 9040) unitfac(kodd), kodd
           WRITE (fuclog, 9040) unitfac(kodd), kodd
   9040    FORMAT (/, ' **** Conversion factor for cgs units ',  &
                'does not equal 1.0 ****', /, 10X,  &
                'Value for UNITFAC(KODD) = ', 1P, g12.4, ' KODD = ', 0P,  &
                i3, /, 10X, 'Check value in file getunits.f90')
           ierr(15) = .TRUE.
        END IF
     END IF
     ! ... write input units and conversion factors for array
     WRITE(fupdef,'(13X,3A,1P,G10.3)') 'Output label for units: ',  &
          unitlabel(kodd), ' Conversion factor: ',  &
          unitfac(kodd)
     IF (unitlabel(kodd) == ht_cgs) THEN
        WRITE(fupdef, '(13X,2A)') 'No conversion necessary, ',  &
             'HYDROTHERM cgs units used for input'
     ELSE
        WRITE(fupdef, '(13X,2A,/,20X,A)')  &
             'Converting input values to HYDROTHERM cgs units ' , unitsavail(1),  &
             '(Value in input units * factor = value in cgs units)'
     END IF
   END SUBROUTINE getunits_3
!
   SUBROUTINE getunits_2
   ! **********
   ! Units allowed for Version 2 and earlier
     IMPLICIT NONE
!      ---set initial values
          KODD = 0
          ldum1 = 0
          ldum2 = 0
! 
!             initialize character strings for available units
          unitsavail(1) = ' '
          unitsavail(2) = ' '
          unitsavail(3) = ' '
          unitsavail(4) = ' '
          unitsavail(5) = ' '
! 
!      ---separate portion of input line between parentheses () containing units
          ldum1 = INDEX(DUMSTRING, '(')
          ldum2 = INDEX(DUMSTRING, ')')
! 
!      ---check to see that input line containes parentheses
          IF (ldum1.EQ.0 .OR. ldum2.EQ.0 .OR. ldum2.EQ.ldum1+1) THEN
            lproblem = .TRUE.
            uunits = '???'
          ELSE
            lproblem = .FALSE.
            ldum1 = ldum1 + 1
            ldum2 = ldum2 - 1
            uunits = DUMSTRING(ldum1:ldum2)
            WRITE (fupdef, '(13X,2A)') 'User specified input units: ', uunits
          ENDIF
! 
! 
! 
!      ---Porosity
! 
          IF (ARRAY.EQ.'PORO') THEN
            KODD = 1
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(unitless)'
            unitsavail(2) = '(%)'
            ht_cgs = '(-)'
! 
!                base units: unitless
            IF (uunits.EQ.'unitless' .OR. uunits.EQ.'none') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(-)'
! 
!                   units: %
            ELSEIF (uunits.EQ.'%' .OR. uunits.EQ.'per cent' .OR.     &
                   uunits.EQ.'Per cent') THEN
              Unitfac(KODD) = 1.0D-2
              Unitlabel(KODD) = '(%)'
            ENDIF
! 
!      ---X, Y, or Z Permeability
! 
          ELSEIF (ARRAY.EQ.'XPER' .OR. ARRAY.EQ.'YPER' .OR. ARRAY.EQ.'ZPER')    &
                 THEN
            IF (ARRAY.EQ.'XPER') KODD = 2
            IF (ARRAY.EQ.'YPER') KODD = 3
            IF (ARRAY.EQ.'ZPER') KODD = 4
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(cm**2)'
            unitsavail(2) = '(m**2)'
            unitsavail(3) = '(darcy)'
            ht_cgs = '(cm^2)'
! 
!              base units: cm2
            IF (uunits.EQ.'cm2' .OR. uunits.EQ.'cm**2' .OR.     &
               uunits.EQ.'cm*2') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(cm^2)'
! 
!                   units: m2
            ELSEIF (uunits.EQ.'m2' .OR. uunits.EQ.'m**2' .OR.     &
                   uunits.EQ.'m*2') THEN
              Unitfac(KODD) = 1.0D4
              Unitlabel(KODD) = '(m^2)'
! 
!                   units: darcy
            ELSEIF (uunits.EQ.'darcy' .OR. uunits.EQ.'Darcy') THEN
              Unitfac(KODD) = 9.87D-9
              Unitlabel(KODD) = '(darcy)'
            ENDIF
! 
!      ---Thermal Conductivity
! 
          ELSEIF (ARRAY.EQ.'THER') THEN
            KODD = 5
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(erg/cm_s_K)'
            unitsavail(2) = '(J/m_s_K)'
            ht_cgs = '(erg/s-cm-K)'
! 
!                  base units: erg/cm_s_K
            IF (uunits.EQ.'erg/cm sec K' .OR. uunits.EQ.'erg/cm_sec_K' .OR.     &
               uunits.EQ.'erg/cm sec C' .OR. uunits.EQ.'erg/cm_sec_C' .OR.     &
               uunits.EQ.'erg/cm s K' .OR. uunits.EQ.'erg/cm_s_K' .OR.     &
               uunits.EQ.'erg/cm s C' .OR. uunits.EQ.'erg/cm_s_C') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(erg/s-cm-K)'
! 
!                   units: W/m K
            ELSEIF (uunits.EQ.'W/m K' .OR. uunits.EQ.'W/m_K' .OR.     &
                   uunits.EQ.'W/m C' .OR. uunits.EQ.'W/m_C') THEN
              Unitfac(KODD) = 1.0D5
              Unitlabel(KODD) = '(W/m-K)'
! 
!                   units: J/m s K
            ELSEIF (uunits.EQ.'J/m sec K' .OR. uunits.EQ.'J/m_sec_K' .OR.     &
                   uunits.EQ.'J/m sec C' .OR. uunits.EQ.'J/m_sec_C' .OR.     &
                   uunits.EQ.'J/m s K' .OR. uunits.EQ.'J/m_s_K' .OR.     &
                   uunits.EQ.'J/m s C' .OR. uunits.EQ.'J/m_s_C') THEN
              Unitfac(KODD) = 1.0D5
              Unitlabel(KODD) = '(J/s-m-K)'
            ENDIF
! 
!      ---Specific Heat (Heat Capacity)
! 
          ELSEIF (ARRAY.EQ.'SPEC') THEN
            KODD = 6
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(erg/g_K)'
            unitsavail(2) = '(J/g_K)'
            unitsavail(3) = '(J/kg_K)'
            unitsavail(4) = '(kJ/kg_K)'
            ht_cgs = '(erg/g-K)'
! 
!                  base units: erg/g_K
            IF (uunits.EQ.'erg/g K' .OR. uunits.EQ.'erg/g_K' .OR.     &
               uunits.EQ.'erg/g C' .OR. uunits.EQ.'erg/g_C') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(erg/g-K)'
! 
!                   units: J/g_K
            ELSEIF (uunits.EQ.'J/g K' .OR. uunits.EQ.'J/g_K' .OR.     &
                   uunits.EQ.'J/g C' .OR. uunits.EQ.'J/g_C') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(J/g-K)'
! 
!                   units: J/kg_K
            ELSEIF (uunits.EQ.'J/kg K' .OR. uunits.EQ.'J/kg_K' .OR.     &
                   uunits.EQ.'J/kg C' .OR. uunits.EQ.'J/kg_C') THEN
              Unitfac(KODD) = 1.0D4
              Unitlabel(KODD) = '(J/kg-K)'
! 
!                   units: kJ/kg_K
            ELSEIF (uunits.EQ.'kJ/kg K' .OR. uunits.EQ.'kJ/kg_K' .OR.     &
                   uunits.EQ.'kJ/kg C' .OR. uunits.EQ.'kJ/kg_C') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(kJ/kg-K)'
            ENDIF
! 
!      ---Rock Density
! 
          ELSEIF (ARRAY.EQ.'RDEN') THEN
            KODD = 7
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(g/cm**3)'
            unitsavail(2) = '(kg/m**3)'
            ht_cgs = '(g/cm^3)'
! 
!                     base units: g/cm3
            IF (uunits.EQ.'g/cm3' .OR. uunits.EQ.'g/cm**3') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(g/cm^3)'
! 
!                   units: kg/m3
            ELSEIF (uunits.EQ.'kg/m3' .OR. uunits.EQ.'kg/m**3') THEN
                      Unitfac(KODD) = 1.0D-3
              Unitlabel(KODD) = '(kg/m^3)'
            ENDIF
! 
!      ---Rock Compressibility
! 
          ELSEIF (ARRAY.EQ.'COMP') THEN
            KODD = 8
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(cm**2/dyne)'
            unitsavail(2) = '(bars**-1)'
            unitsavail(3) = '(Pa**-1)'
            unitsavail(4) = '(MPa**-1)'
            ht_cgs = '(cm^2/dyne)'
! 
!                     base units: cm**2/dyne
            IF (uunits.EQ.'cm2/dyne' .OR. uunits.EQ.'cm**2/dyne') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(cm^2/dyne)'
! 
!                  units: bars**-1
            ELSEIF (uunits.EQ.'1/bars' .OR. uunits.EQ.'1/bar' .OR.     &
                   uunits.EQ.'bars**-1' .OR. uunits.EQ.'bar**-1' .OR.     &
                   uunits.EQ.'bars-1' .OR. uunits.EQ.'bar-1') THEN
              Unitfac(KODD) = 1.0D-6
              Unitlabel(KODD) = '(bar^-1)'
! 
!                  units: Pa**-1
            ELSEIF (uunits.EQ.'1/Pa' .OR. uunits.EQ.'Pa**-1' .OR.     &
                   uunits.EQ.'Pa-1') THEN
              Unitfac(KODD) = 1.0D-1
              Unitlabel(KODD) = '(Pa^-1)'
! 
!                  units: MPa**-1
            ELSEIF (uunits.EQ.'1/MPa' .OR. uunits.EQ.'MPa**-1' .OR.     &
                   uunits.EQ.'MPa-1') THEN
              Unitfac(KODD) = 1.0D-7
              Unitlabel(KODD) = '((MPa^-1)'
            ENDIF
! 
!      ---Conductive Heat Flux along the Base
! 
          ELSEIF (ARRAY.EQ.'COND') THEN
            KODD = 9
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(erg/s_cm**2)'
            unitsavail(2) = '(mW/m**2)'
            unitsavail(3) = '(W/m**2)'
            unitsavail(4) = '(J/s_m**2)'
            ht_cgs = '(erg/s-cm^2)'
! 
!              base units: erg/s_cm**2
            IF (uunits.EQ.'ergs/sec_cm**2' .OR.     &
               uunits.EQ.'ergs/sec_cm2' .OR. uunits.EQ.'ergs/s_cm**2' .OR.     &
               uunits.EQ.'ergs/s_cm2' .OR. uunits.EQ.'erg/sec_cm**2' .OR.     &
               uunits.EQ.'erg/sec_cm2' .OR. uunits.EQ.'erg/s_cm**2' .OR.     &
               uunits.EQ.'erg/s_cm2' .OR. uunits.EQ.'ergs/sec cm**2' .OR.     &
               uunits.EQ.'ergs/sec cm2' .OR. uunits.EQ.'ergs/s cm**2' .OR.     &
               uunits.EQ.'ergs/s cm2' .OR. uunits.EQ.'erg/sec cm**2' .OR.     &
               uunits.EQ.'erg/sec cm2' .OR. uunits.EQ.'erg/s cm**2' .OR.     &
               uunits.EQ.'erg/s cm2') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(erg/s-cm^2)'
! 
!                  units:  mW/m**2
            ELSEIF (uunits.EQ.'mW/m**2' .OR. uunits.EQ.'mW/m2') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(mW/m^2)'
! 
!                  units:  W/m**2
            ELSEIF (uunits.EQ.'W/m**2' .OR. uunits.EQ.'W/m2') THEN
              Unitfac(KODD) = 1.0D3
              Unitlabel(KODD) = '(W/m^2)'
! 
!                  units: J/s_m**2
            ELSEIF (uunits.EQ.'J/sec_m**2' .OR. uunits.EQ.'J/sec_m2' .OR.     &
                   uunits.EQ.'J/s_m**2' .OR. uunits.EQ.'J/s_m2' .OR.     &
                   uunits.EQ.'J/sec m**2' .OR. uunits.EQ.'J/sec m2' .OR.     &
                   uunits.EQ.'J/s m**2' .OR. uunits.EQ.'J/s m2') THEN
              Unitfac(KODD) = 1.0D3
              Unitlabel(KODD) = '(J/s-m^2)'
            ENDIF
! 
!      ---Enthalpy
! 
          ELSEIF (ARRAY.EQ.'ENTH') THEN
            KODD = 10
!              set default units for temperature when enthalpy data input
            Unitfac(15) = 1._kdp
            Unitlabel(15) = '(C)'
            WRITE (fupdef, '(13X,A,A)')     &
                         'Default output units assigned for Temperature: '    &
                         , Unitlabel(15)
! 
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(erg/g)'
            unitsavail(2) = '(J/g)'
            unitsavail(3) = '(kJ/kg)'
            unitsavail(4) = '(J/kg)'
            ht_cgs = unitsavail(1)
! 
!                 base units: erg/g
            IF (uunits.EQ.'ergs/g' .OR. uunits.EQ.'erg/g') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(erg/g)'
! 
!                  units: J/g
            ELSEIF (uunits.EQ.'J/g') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(J/g)'
! 
!                  units: kJ/kg
            ELSEIF (uunits.EQ.'kJ/kg') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(kJ/kg)'
! 
!                  units: J/kg
            ELSEIF (uunits.EQ.'J/kg') THEN
              Unitfac(KODD) = 1.0D4
              Unitlabel(KODD) = '(J/kg)'
            ENDIF
! 
!      ---Pressure
! 
          ELSEIF (ARRAY.EQ.'PRES') THEN
            KODD = 11
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(dyne/cm**2)'
            unitsavail(2) = '(MPa)'
            unitsavail(3) = '(bars)'
            unitsavail(4) = '(Pa)'
            unitsavail(5) = '(atm)'
            ht_cgs = '(dyne/cm^2)'
! 
!                 base units: dyne/cm**2
            IF (uunits.EQ.'dynes/cm2' .OR. uunits.EQ.'dynes/cm**2' .OR.     &
               uunits.EQ.'dyne/cm2' .OR. uunits.EQ.'dyne/cm**2') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(dyne/cm^2)'
! 
!                   units: MPa
            ELSEIF (uunits.EQ.'MPa') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(MPa)'
! 
!                   units: bars
            ELSEIF (uunits.EQ.'bars' .OR. uunits.EQ.'bar') THEN
              Unitfac(KODD) = 1.0D6
              Unitlabel(KODD) = '(bar)'
! 
!                   units: Pa
            ELSEIF (uunits.EQ.'Pa') THEN
              Unitfac(KODD) = 1.0D1
              Unitlabel(KODD) = '(Pa)'
! 
!                   units: atm
            ELSEIF (uunits.EQ.'atm') THEN
              Unitfac(KODD) = 1.01325_kdp
              Unitlabel(KODD) = '(atm)'
            ENDIF
! 
!      ---X, Y, or Z  Spacing
! 
          ELSEIF (ARRAY.EQ.'XSPA' .OR. ARRAY.EQ.'YSPA' .OR. ARRAY.EQ.'ZSPA')    &
                 THEN
            IF (ARRAY.EQ.'XSPA') KODD = 12
            IF (ARRAY.EQ.'YSPA') KODD = 13
            IF (ARRAY.EQ.'ZSPA') KODD = 14
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(cm)'
            unitsavail(2) = '(m)'
            unitsavail(3) = '(km)'
            ht_cgs = unitsavail(1)
! 
!                   base units: cm
            IF (uunits.EQ.'cm') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(cm)'
! 
!                   units: m
            ELSEIF (uunits.EQ.'m') THEN
              Unitfac(KODD) = 1.0D2
              Unitlabel(KODD) = '(m)'
! 
!                   units: km
            ELSEIF (uunits.EQ.'km') THEN
              Unitfac(KODD) = 1.0D5
              Unitlabel(KODD) = '(km)'
            ENDIF
! 
!      ---Temperature
! 
          ELSEIF (ARRAY.EQ.'TEMP') THEN
            KODD = 15
!              set default units for enthalpy when temperature data input
            Unitfac(10) = 1.0e7_kdp
            Unitlabel(10) = '(kJ/kg)'
            WRITE (fupdef, '(13X,A,A)')     &
                            'Default output units assigned for Enthalpy: '    &
                            , Unitlabel(10)
! 
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(C)'
            ht_cgs = '(Deg.C)'
! 
!                   base units: C
            IF (uunits.EQ.'C' .OR. uunits.EQ.'deg C' .OR.     &
               uunits.EQ.'deg_C' .OR. uunits.EQ.'degrees Celsius' .OR.     &
               uunits.EQ.'degrees_Celsius' .OR. uunits.EQ.'degree C' .OR.     &
               uunits.EQ.'degree_C' .OR. uunits.EQ.'degrees C' .OR.     &
               uunits.EQ.'degrees_C') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(Deg.C)'
            ENDIF
! 
!      ---Source/sink input 1 - Mass input rate
! 
          ELSEIF (ARRAY.EQ.'SRC1') THEN
            KODD = 16
! 
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(g/s)'
            unitsavail(2) = '(kg/s)'
            unitsavail(3) = '(g/yr)'
            unitsavail(4) = '(kg/yr)'
            ht_cgs = unitsavail(1)
! 
!                 base units: g/s
            IF (uunits.EQ.'g/s' .OR. uunits.EQ.'g/sec') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(g/s)'
! 
!                  units: kg/s
            ELSEIF (uunits.EQ.'kg/s' .OR. uunits.EQ.'kg/sec') THEN
              Unitfac(KODD) = 1.0D3
              Unitlabel(KODD) = '(kg/s)'
! 
!                  units: g/yr
            ELSEIF (uunits.EQ.'g/yr' .OR. uunits.EQ.'g/year') THEN
              Unitfac(KODD) = 3.17D-8
              Unitlabel(KODD) = '(g/yr)'
! 
!                  units: kg/yr
            ELSEIF (uunits.EQ.'kg/yr' .OR. uunits.EQ.'kg/year') THEN
              Unitfac(KODD) = 3.17D-5
              Unitlabel(KODD) = '(kg/yr)'
            ENDIF
! 
!      ---Source/sink input 2 - Enthalpy / Temperature of Source
! 
          ELSEIF (ARRAY.EQ.'SRC2') THEN
            KODD = 17
! 
!                  set initalial value for conversion factor
            Unitfac(KODD) = 0._kdp
!                  available input units
            unitsavail(1) = '(erg/g)'
            unitsavail(2) = '(J/g)'
            unitsavail(3) = '(kJ/kg)'
            unitsavail(4) = '(J/kg)'
            unitsavail(5) = '(C)'
            ht_cgs = unitsavail(1)
! 
! 
!                 base units: erg/g
            IF (uunits.EQ.'ergs/g' .OR. uunits.EQ.'erg/g') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(erg/g)'
! 
!                  units: J/g
            ELSEIF (uunits.EQ.'J/g') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(J/g)'
! 
!                  units: kJ/kg
            ELSEIF (uunits.EQ.'kJ/kg') THEN
              Unitfac(KODD) = 1.0e7_kdp
              Unitlabel(KODD) = '(kJ/kg)'
! 
!                  units: J/kg
            ELSEIF (uunits.EQ.'J/kg') THEN
              Unitfac(KODD) = 1.0D4
              Unitlabel(KODD) = '(J/kg)'
! 
!                  units: C
            ELSEIF (uunits.EQ.'C' .OR. uunits.EQ.'deg C' .OR.     &
                   uunits.EQ.'deg_C' .OR. uunits.EQ.'degrees Celsius' .OR.     &
                   uunits.EQ.'degrees_Celsius' .OR.     &
                   uunits.EQ.'degree C' .OR. uunits.EQ.'degree_C' .OR.     &
                   uunits.EQ.'degrees C' .OR. uunits.EQ.'degrees_C') THEN
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = '(Deg.C)'
              ht_cgs = '(Deg.C)'
            ENDIF
! 
!      ---end of array types
          ENDIF
! 
!      ---check for errors
! 
!               if array line did not contain units enclosed in parentheses
          IF (lproblem) THEN
!              if old input format is used, write warning and use cgs units,
            IF (Informt.EQ.1) THEN
              write (fuclog, 9005) 'WARNING - Problem', ARRAY
              write (fuclog, 9025) unitsavail(1)
              write (fuclog, 9035) (unitsavail(lll), lll=2, 5)
              write (fuclog, 9010)
              Unitfac(KODD) = 1._kdp
              Unitlabel(KODD) = unitsavail(1)
 
              write (fustdout, 9005) 'WARNING - Problem', ARRAY
              write (fustdout, 9025) unitsavail(1)
              write (fustdout, 9035) (unitsavail(lll), lll=2, 5)
              write (fustdout, 9010)
              ierr(11) = .TRUE.
            ELSE
              write (fustdout, 9005) '**** ERROR', ARRAY
              write (fustdout, 9010)
              write (fustdout, 9030) ARRAY, unitsavail(1)
              write (fustdout, 9035) (unitsavail(lll), lll=2, 5)
              write (fuclog, 9005) '**** ERROR', ARRAY
              write (fuclog, 9010)
              write (fuclog, 9030) ARRAY, unitsavail(1)
              write (fuclog, 9035) (unitsavail(lll), lll=2, 5)
             ierr(12) = .TRUE.
             RETURN
            ENDIF
          ENDIF
! 
      IF (KODD.EQ.0) THEN
      write (fuclog, 9015) ARRAY
      write (fustdout, 9015) ARRAY
        ierr(13) = .TRUE.
        RETURN
      ENDIF
      IF (Unitfac(KODD).EQ.0._kdp) THEN
         write (fustdout, 9020) ARRAY, uunits
         write (fustdout, 9030) ARRAY, unitsavail(1)
         write (fustdout, 9035) (unitsavail(lll), lll=2, 5)
         write (fuclog, 9020) ARRAY, uunits
         write (fuclog, 9030) ARRAY, unitsavail(1)
         write (fuclog, 9035) (unitsavail(lll), lll=2, 5)
        ierr(14) = .TRUE.
        RETURN
      ENDIF
     ! ...     if cgs units are used, check that conversion factor =1
      IF (Unitlabel(KODD).EQ.ht_cgs) THEN
         IF (ABS(Unitfac(KODD)-1._kdp).GT.1.0e-8_kdp) THEN
            write (fustdout, 9040) Unitfac(KODD), KODD
            write (fuclog, 9040) Unitfac(KODD), KODD
           ierr(15) = .TRUE.
         ENDIF
      ENDIF
!      ---write input units and conversion factors for array
      WRITE (fupdef, '(13X,3A,1P,G10.3)') 'Output label for units: ',     &
                                    Unitlabel(KODD),     &
                                    ' Conversion factor: ',     &
                                    Unitfac(KODD)
      IF (Unitlabel(KODD).EQ.ht_cgs) THEN
        WRITE (fupdef, '(13X,2A)') 'No conversion necessary, ',     &
                              'HYDROTHERM cgs units used for input'
      ELSE
        WRITE (fupdef, '(13X,2A,/,20X,A)')     &
                    'Converting input values to HYDROTHERM cgs units '    &
                    , ht_cgs,     &
              '(Value in input units * factor = value in cgs units)'
      ENDIF
   9005    FORMAT (/, 5X, ' **** ', a, ' determining units for array ', a,  &
                ' *****')
   9025    FORMAT (13X, 'Using HYDROTHERM cgs units ', a)
   9010    FORMAT (13X, 'Check that units are enclosed in parentheses and',  &
                /, 13X, 'that they are on the same line as the array name')
   9015 FORMAT (/, ' **** Problem with input array named: ', a,  &
             ' *****', /, 6X, 'Check array name')
   9020 FORMAT (/, ' **** ERROR determining units for array ', a,  &
             ' *****', /, 6X, 'Input units - ', a, ' not appropriate')
   9030 FORMAT (13X, 'HYDROTHERM cgs units for array ', a, ' are ', a)
   9035 FORMAT (13X, 'Other available input units are:', /, 11X, 5(2X,a))
   9040    FORMAT (/, ' **** Conversion factor for cgs units ',  &
                'does not equal 1.0 ****', /, 10X,  &
                'Value for UNITFAC(KODD) = ', 1P, g12.4, ' KODD = ', 0P,  &
                i3, /, 10X, 'Check value in file getunits.f90')
   END SUBROUTINE getunits_2
END SUBROUTINE getunits
