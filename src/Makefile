# compiler
FC = gfortran
BIN=../
# compile flags
FCFLAGS = -g -c -fdefault-real-8 -freal-8-real-16 #-fbacktrace -fno-align-commons
# link flags
FLFLAGS = -g -fbacktrace

SRC_MOD=*_mod.f90 
 # program name
PROGRAM = HYDROTHERMAL
# follow the following steps make one by one
modules_mach: mach_mod.f90 modules.f90
	${FC} ${FCFLAGS} mach_mod.f90 modules.f90
mod: *_mod.f90 
	${FC} ${FCFLAGS} ${SRC_MOD}
fmod: *.f 
	${FC} ${FCFLAGS} *.f
obj: *.f90
	${FC} ${FCFLAGS} *.f90
program: *.o
	$(FC) -o ${PROGRAM} $(FLFLAGS) *.o
	mv ${PROGRAM} ${BIN}
clean:
	rm -f *.o *.mod
