PROGRAM hydrotherm
  ! ...                 ************************
  ! ...                 HYDROTHERM, version  3.2
  ! ...                 ************************
  ! ... Three-dimensional finite-difference simulator of 
  ! ...       ground-water flow and heat transport in the temperature
  ! ...       range of 0 to 1,200 degrees Celsius
  ! ... Purpose:  To simulate vapor- or liquid-dominated hydrothermal
  ! ...           systems in three-dimensions using finite difference techniques.
  ! ...           Direct matrix solution is employed for X-Z, 2-d regions and
  ! ...           one-dimensional regions.
  ! ...           Iterative matrix solution is employed for 3-d regions using 
  ! ...           slice successive overrelaxation methd (SSOR) or 
  ! ...           generalized minimum-residual method (GMRES).
  ! ...           Newton-Raphson iterations are used
  ! ...           to solve non-linear mass and energy governing
  ! ...           equations simultaneously for pressure and enthalpy.
  !
  ! ... Originally programed by Charles R. Faust and James W. Mercer
  ! ...               Geotrans, Inc., 1980
  ! ... Modified by Steven E. Ingebritsen
  ! ...               AES Dept., Stanford University, 1983
  ! ... Rewritten by Daniel O. Hayba
  ! ...               U.S. Geological Survey, 1990
  ! ... Extensively modified, extending the temperature range to 1200 C
  ! ...       and the pressure range to 10 kb
  ! ...               Daniel O. Hayba and Steven E. Ingebritsen
  ! ...               U.S. Geological Survey, 1993
  ! ... Modified to handle unconfined flow with a stagnant air phase
  ! ... Restructured main program to allow replacement with GUI driver
  ! ... Converted to Fortran 90
  ! ... Modified by addition of seepage surface boundary condition and
  ! ...     precipitation recharge boundary condition
  ! ...               Kenneth L. Kipp
  ! ...               U.S. Geological Survey, 2008
  ! ...                            **********
  ! ...                        DISCLAIMER STATEMENT
  ! ... Although this program has been used by the USGS, no warranty, expressed 
  ! ... or implied, is made by the USGS or the United States Government as to 
  ! ... the accuracy and functioning of the program and related program material 
  ! ... nor shall the fact of distribution constitute any such warranty, and no 
  ! ... responsibility is assumed by the USGS in connection therewith.
  ! ...                            **********
  USE control
  USE units
  IMPLICIT NONE
  INTEGER :: flag, i1, i2
  ! ... Set string for use with RCS ident command
  CHARACTER(LEN=80) ::  &
       ident_string='$Name: v3_1_1 $//$Revision: 1.16 $//$Date: 2009/03/31 21:06:19 $'
  !   ---------------------------------------------------------------------------------
  ! ... Extract the version name for the header
  CALL version
  ! ... Set a flag to indicate this is a standalone run without GUI
  gui = .false.
  ! ... Initialize variables and parameters
  flag = 0
  CALL init_ht(flag)
  ! ... Time loop 
  DO WHILE(flag == 0)
     CALL time_step_ht(flag)
  END DO
  ! ... Termination procedure
  CALL terminate_ht(flag)
  CALL dealloc_arr
  STOP 'Simulation Completed'
END PROGRAM hydrotherm
