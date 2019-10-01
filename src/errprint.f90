SUBROUTINE errprint(ie1,ie2)
  ! ... Prints the index numbers of errors encountered and a brief message
  ! ... *****DO NOT RENUMBER THIS ROUTINE****
  USE f_units
  USE control
  USE units, ONLY: phifile
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ie1, ie2
  !
  INTEGER :: ie
  CHARACTER (LEN=120) :: errline, errline2
  ! ... Set string for use with RCS ident command
  !CHARACTER(LEN=80) :: ident_string='$Revision: 2.2 $//$Date: 2007/10/16 21:17:22 $'
  !------------------------------------------------------------------------
  DO  ie=ie1,ie2
     IF (ierr(ie)) GO TO 20
  END DO
  RETURN
20 WRITE(fuclog,2001) '**** The Following Errors Were Detected ****'
2001 FORMAT(//tr40,a)
  ! ... The error messages
  IF(ierr(1)) THEN
     errline = '1  - GDATA - Problem with node sequence numbers'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
9001 FORMAT(tr5,a)
  END IF
  IF(ierr(2)) THEN
     errline = '2 - GFILES - Input file does not exist'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(3)) THEN
     errline = '3 - GFILES - Input file for i.c. pressure and '//  &
          'porosity does not exist'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(4)) THEN
     errline = '4 - GFILES - Invalid value for NPLACE '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(5)) THEN
     errline = '5 - GFILES - Unable to locate hydrotherm tables'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(6)) THEN
     errline = '6 - GDATA - First line of input data file does '//  &
          'not begin with TITLE'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(7)) THEN
     errline = '7 - GDATA - Unable to recognize format flag for '//  &
          'input data file'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(8)) THEN
     errline = '8 - GDATA - Problem grid dimensions larger than '//  &
          'compiled limits'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(9)) THEN
     errline = '9 - GDATA - MPNXYYZ Dimension must be greater '// 'than one'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(10)) THEN
     errline = '10 - GDATA - Initial condition P & PHI should '//  &
          'not be input from file '
     WRITE(fuclog,9010) errline,TRIM(phifile),  &
          'because either BETA = 0 or start time = 0'
     WRITE(fustdout,9010) errline,TRIM(phifile),  &
          'because either BETA = 0 or start time = 0'
9010 FORMAT(tr5,a/tr20,a/tr10,a)
  END IF
  IF(ierr(11)) THEN
     errline = '11 - GETUNITS - Units problem with array input '//  &
          'Will default to base cgs units'
     WRITE(fuclog,9011) errline
     WRITE(fustdout,9011) errline
9011 FORMAT(tr5,a,tr2,a)
  END IF
  IF(ierr(12)) THEN
     errline = '12 - GETUNITS - Units problem with array input '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(13)) THEN
     errline = '13 - GETUNITS - Array name not recognized'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(14)) THEN
     errline = '14 - GETUNITS - Array units not recognized '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(15)) THEN
     errline = '15 - GETUNITS - Conversion factor error for '// 'array units'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(16)) THEN
     errline = '16 - INPUTPH - Initial conditions for '//  &
          'temperature or enthalpy are incomplete'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(17)) THEN
     errline = '17 - INPUTPH - Pressure must not be redefined '//  &
          'for temperature calculation option 33'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(19)) THEN
     errline = '19 - GDATA - Factor for maximum time step '//  &
          'increase must be greater than one'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(20)) THEN
     errline = '20 - GDATA - Invalid index for permeability '// 'function'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(21)) THEN
     errline = '21 - RD - A parameter is using ROCK type input '//  &
          'which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(22)) THEN
     errline = '22 - RD - A parameter is using RTFTEMP type input '//  &
          'which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(23)) THEN
     errline = '23 - RD - Parameter values not in order of '//  &
          'increasing temperature'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(24)) THEN
     errline = '24 - RD - Rock type index is undefined for a parameter'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(25)) THEN
     errline = '25 - RDSPACE - A parameter is using a type '//  &
          'of input which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(26)) THEN
     errline = '26 - RDSPACE - Array dimensions of zero are invalid'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(27)) THEN
     errline = '27 - READARRY - Heat conduction data not at '//  &
          'bottom boundary of mesh'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(28)) THEN
     errline = '28 - READARRY - Temperature for precipitation '//  &
          ' flux data not at '// 'top boundary of mesh'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(29)) THEN
     errline = '29 - READARRY - An Index of Rock type is '//  &
          'undefined for the simulation region'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(30)) THEN
     errline = '30 - READARRY - Input temperature is outside '//  &
          'of table limits'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(31)) THEN
     errline = '31 - READARRY - Input enthalpy is outside '// 'of table limits'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(32)) THEN
     errline = '32 - READARRY - Input temperature for '//  &
          'precipitation flux is outside '// 'of table limits'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(33)) THEN
     errline = '33 - READARRY - Insufficient number of rock '//  &
          'types defined for RTFTEMP input'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(34)) THEN
     errline = '34 - READARRY - Missing TOP or BOTTOM data record'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(35)) THEN
     errline = '35 - GDATA - Incorrect number of nodes in an '// 'x-z slice'
     errline2 = 'Sum of nodes of each type not equal to total '//  &
          ' number in slice'
     WRITE(fuclog,9035) errline, errline2
     WRITE(fustdout,9035) errline, errline2
9035 FORMAT(tr5,a/tr10,a)
  END IF
  IF(ierr(36)) THEN
     errline = '36 - GDATA - NEWFMT must be integer not logical'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(37)) THEN
     errline = '37 - GDATA - Value for bandwidth must not be '// 'specified'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(38)) THEN
     errline = '38 - SINK - Time step of zero specified'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(39)) THEN
     errline = '39 - SINK - Ending time of simulation period is less than '//  &
          'starting time'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(40)) THEN
     errline = '40 - SINK - Time step of -1 specified for first simulation period'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(41)) THEN
     errline = '41 - SINK - NSRCE specified as -1 for first simulation period'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(42)) THEN
     errline = '42 - SINK - Input type keyword not recognized'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(43)) THEN
     errline = '43 - SINK - Attempt to redefine cell dimensions'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(44)) THEN
     errline = '44 - INPUTPH - Hydrostatic pressure gradient is '//  &
          'required for boiling point with depth'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(45)) THEN
     errline = '45 - INPUTPH - Pressure must not be redefined '//  &
          'when assigning nodes to boiling temperature'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(46)) THEN
     errline = '46 - READARRY - Y-dimensions of cells must not '//  &
          'be defined for cylindrical coordinates'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(47)) THEN
     errline = '47 - READARRY - Pressure initial/boundary '//  &
          'conditions is outside range of table (too high)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(48)) THEN
     errline = '48 - READARRY - Precipitation flux must be '//  &
          'specified at the top boundary of the mesh'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(49)) THEN
     errline = '49 - READARRY - Array name not recognized'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(50)) THEN
     errline = '50 - READARRY - Spatial variation in conductive '//  &
          'heat flux not allowed for one cell in x-direction'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(51)) THEN
     errline = '51 - RELPERM - Invalid index for permeability '// 'function'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(52)) THEN
     errline = '52 - READARRY - Input surface temperature '//  &
          'outside the areal range of the mesh'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(53)) THEN
     errline = '53 - READARRY - More input surface temperature '//  &
          'values than the number of surface cells'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(54)) THEN
     errline = '54 - READARRY - Pressure for boiling point with '//  &
          'depth profile is outside the table range (too high)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(55)) THEN
     errline = '55 - READARRY - Input surface temperature '//  &
          'outside the areal range of the mesh for boiling '//  &
          'point with depth profile'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(56)) THEN
     errline = '56 - READARRY - Pressure at top of mesh '//  &
          'is outside the table range (too high)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(57)) THEN
     errline = '57 - READARRY - Spatial variation in '// 'precipitation '//  &
          'not allowed for one cell in x-direction'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(58)) THEN
     errline = '58 - READARRY - Spatial variation in '//  &
          'associated temperature for precipitation '//  &
          'not allowed for one cell in x-direction'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(59)) THEN
     errline = '59 - READARRY - Invalid input format for array '//  &
          'type being processed'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(60)) THEN
     errline = '60 - READARRY - Invalid option flag for '//  &
          'hydrostatic pressure distribution'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(61)) THEN
     errline = '61 - RD - Invalid option flag for '//  &
          'RTFTEMP distribution of a parameter'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(62)) THEN
     errline = '62 - RD - Number of points defining RTFTEMP '//  &
          'function exceeds limit of 4'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(63)) THEN
     errline = '63 - RD - A parameter is using a type '//  &
          'of input which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(64)) THEN
     errline = '64 - RD - A parameter is using a type '//  &
          'of input format which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(65)) THEN
     errline = '65 - RD - Insufficient number of rock '//  &
          'types have been defined'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(66)) THEN
     errline = '66 - GDATA - Missing TOP or BOTTOM data record'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(67)) THEN
     errline = '67 - READARRY - Only 5 linear segments of '//  &
          'temperature or pressure with depth '//'are allowed'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(68)) THEN
     errline = '68 - READARRY - Precipitation flux too large for permeability'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(70)) THEN
     errline = '70 - bndyT0 - Pressure outside of cubic spline '// 'table limits'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(71)) THEN
     errline = '71 - PHREG - Pressure or enthalpy is outside '//  &
          'table limits for phase diagram'
     WRITE(fuclog,9071) errline
     WRITE(fustdout,9071) errline
9071 FORMAT(tr5,a,tr2,a)
  END IF
  IF(ierr(72)) THEN
     errline = '72 - PHSBNKLO - Pressure is outside table '//  &
          'limits for spline function'
     WRITE(fuclog,9072) errline
     WRITE(fustdout,9072) errline
9072 FORMAT(tr5,a,tr2,a)
  END IF
  IF(ierr(73)) THEN
     errline = '73 - PRESBOIL -Error computing boiling point '//  &
          'with depth profile.'
     errline2 = 'Pressure at top greater than critical point.'
     WRITE(fuclog,9073) errline, errline2
     WRITE(fustdout,9073) errline, errline2
9073 FORMAT(tr5,a/tr10,a)
  END IF
  IF(ierr(74)) THEN
     errline = '74 - TBLOOKUP -Error computing average water '// 'properties'
     WRITE(fuclog,9073) errline
     WRITE(fustdout,9073) errline
  END IF
  IF(ierr(75)) THEN
     errline = '75 - TBLOOKUP -Error computing average water '//  &
          'and steam properties '
     WRITE(fuclog,9073) errline
     WRITE(fustdout,9073) errline
  END IF
  IF(ierr(76)) THEN
     errline = '76 - RD - A parameter is using RTFTIME type input '//  &
          'which is not available'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(77)) THEN
     errline = '77 - RD - Invalid option flag for '//  &
          'RTFTIME function of a parameter'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(78)) THEN
     errline = '78 - RD - Number of points defining RTFTIME '//  &
          'function exceeds limit of 4 or is less than 2'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(79)) THEN
     errline = '79 - RD - Parameter values for RTFTIME not in order of '//  &
          'increasing time'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(81)) THEN
     errline = '81 - GMRES_ASSEMBLE_SOLVE - Linear equation scaling error; '//  &
          'Row or column number ierr is all zeros'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(82)) THEN
     errline = '82 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'Input matrix wrong. Length of row in L or U greater than n'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(83)) THEN
     errline = '83 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'L matrix overflows ALU storage'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(84)) THEN
     errline = '84 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'U matrix overflows ALU storage'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(85)) THEN
     errline = '85 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'Illegal value for lfil'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(86)) THEN
     errline = '86 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'Row of zeros encountered during factorization'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(87)) THEN
     errline = '87 - GMRES_ASSEMBLE_SOLVE - Preconditioning failure; '//  &
          'Zero pivot encountered during factorization step ier'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(88)) THEN
     errline = '88 - GMRES_ASSEMBLE_SOLVE - Solver failure; '//  &
          'Iteration limit reached without convergence'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(89)) THEN
     errline = '89 - GMRES_ASSEMBLE_SOLVE - Solver failure; '//  &
          'Initial solution gives residual of zero'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(90)) THEN
     errline = '90 - FORMEQ - Linear equation scaling failure; '//  &
          'A zero scale factor encountered'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(91)) THEN
     errline = '91 - SINK - Time step specified less than minimum allowed'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(101)) THEN
     errline = '101 - PRPTY - Porosity is out of range (0-1)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(102)) THEN
     errline = '102 - PRPTY - Invalid index of thermodynamic region'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(103)) THEN
     errline = '103 - SINK - Error in ichgpr parameter'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(104)) THEN
     errline = '104 - SINKNODE - Input enthalpy of fluid must be '//  &
          '0 for sink nodes'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(105)) THEN
     errline = '105 - SINKNODE - One or more nodal indices '//  &
          'is not in the region '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(106)) THEN
     errline = '106 - SINKNODE - Error in NODE input for point source/sink'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(107)) THEN
     errline = '107 - SINKWELL - Input enthalpy of fluid must be '//  &
          '0 for sink wells'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(108)) THEN
     errline = '108 - SINKWELL - Well location '//  &
          'indices are not in the region '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(109)) THEN
     errline = '109 - SINKWELL - Error in well screen elevation '// 'index'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(110)) THEN
     errline = '110 - TIME_STEP_HT - Time step cut more than '//  &
          'maximum number of times allowed'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(111)) THEN
     errline = '111 - TIME_STEP_HT - Time step is less than '// 'minimum allowed'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(112)) THEN
     errline = '112 - SINKNODE - Node source not allowed in cell column with well '//  &
          'source '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(113)) THEN
     errline = '113 - SINKWELL - Well source not allowed in cell column with node '//  &
          'source '
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(114)) THEN
     errline = '114 - SINKWELL - Error in WELL input for source/sink'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(121)) THEN
     errline = '121 - AVGVALUE - Unable to calculate average '//  &
          'pressure and enthalpy at cell face'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(122)) THEN
     errline = '122 - WRNEW - Conversion factor of zero'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(123)) THEN
     errline = '123 - TBLOOKUP - Pressure outside of table (high)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(124)) THEN
     errline = '124 - TBLOOKUP - Enthalpy outside of table (high)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(125)) THEN
     errline = '125 - TIME_STEP_HT - Time step criterion flag is '// 'invalid: 0'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(126)) THEN
     errline = '126 - TIME_STEP_HT - N-R convergence failure'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(127)) THEN
     errline = '127 - PRINTOPT - Invalid indices of location '//  &
          'for writing of temporal plot data'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(128)) THEN
     errline = '128 - GDATA - Invalid option number for linear equation '//  &
          'solver'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(129)) THEN
     errline = '129 - GDATA - SSOR (Direct) solver must be selected '//  &
          'for a one-cell problem'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(130)) THEN
     errline = '130 - READARRY - Input surface pressure '//  &
          'outside the areal range of the mesh'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(131)) THEN
     errline = '131 - READARRY - More input surface pressure '//  &
          'values than the number of surface cells'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(132)) THEN
     errline = '132 - TBLOOKUP - Enthalpy outside of table (low)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(133)) THEN
     errline = '133 - WATER_PROP_LOP - Pressure too high for functions'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(134)) THEN
     errline = '134 - WATER_PROP_LOP - Enthalpy too high for functions'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(135)) THEN
     errline = '135 - TBLOOKUP - Pressure outside of table (low)'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(172)) THEN
     errline = '172 - READARRY - Function no. 3 for temperature-depencence of '//  &
                 'specific heat is currently disabled. Please use an alternative function.'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(173)) THEN
     errline = '173 - READARRY - Option no. 22 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(174)) THEN
     errline = '174 - READARRY - Option no. 22 for calculating temperature '//  &
                 'distribution is currently disabled. Please use an alternative option.'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(175)) THEN
     errline = '175 - WT_ELEV - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(176)) THEN
     errline = '176 - WRNEW - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(177)) THEN
     errline = '177 - WELLALLO - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(178)) THEN
     errline = '178 - SINKWELL - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(179)) THEN
     errline = '179 - SINK - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(180)) THEN
     errline = '180 - READARRY - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(181)) THEN
     errline = '181 - RDSPACE - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(182)) THEN
     errline = '182 - RD - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(183)) THEN
     errline = '183 - PLOTXYZ - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(184)) THEN
     errline = '184 - PLOTIDL - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(185)) THEN
     errline = '185 - PLOTGNU - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(186)) THEN
     errline = '186 - PLOTEXPL - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(187)) THEN
     errline = '187 - PLOTB2 - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF

  IF(ierr(188)) THEN
     errline = '188 - PDATA - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(189)) THEN
     errline = '189 - SSOR - Deallocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(190)) THEN
     errline = '190 - HYDROTHERM - Maximum number of time steps reached'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(191)) THEN
     errline = '191 - GMRES - Deallocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(192)) THEN
     errline = '192 - GMRES_ASEM_SOLVE - Deallocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(193)) THEN
     errline = '193 - ILUPC - Deallocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(194)) THEN
     errline = '194 - ILUPC - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(195)) THEN
     errline = '195 - GMRES_ASEM_SOLVE - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(196)) THEN
     errline = '196 - GMRES - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(197)) THEN
     errline = '197 - SSOR - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(198)) THEN
     errline = '198 - READARRY - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(199)) THEN
     errline = '199 - GDATA - Allocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
  IF(ierr(200)) THEN
     errline = '200 - DEALLOC_ARR - Deallocation error'
     WRITE(fuclog,9001) errline
     WRITE(fustdout,9001) errline
  END IF
END SUBROUTINE errprint
