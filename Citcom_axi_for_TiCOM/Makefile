#
#  Makefile for "citcom" allowing considerable variations
#

#
#  NOTES:
#  1) Macro substitution is different for different version of make, for SGI, CRAY, HP
#  It seems necessary to go through and replace the CEXT ... FOBJFLAG macros manually.
#
#  2) Put FLAGS=$(YOUR_MACHINE_FLAGS) chosen from the list that follows, add machines
#  where necessary, try and find a better way than just defining "const" to be null. 
#
#  3) Which flags give best performance ?  I'm doing my best with the architectures	
#  as I try them.

##############################################################################
#	Individual Machine Variations
##############################################################################

#FFTLIB= /home/arthur/d5/zhong/lib/fft.a /home/arthur/d5/zhong/lib/realfft.a
COMPRESS= /usr/bin/compress
LIB_PATH=
LIB_LIST=
LIB= $(LIB_PATH) $(LIB_LIST) -lm


###################################################################
#	Operating System Variables
###################################################################

CC=cc #/usr/lib/cmplrs/cc/gemc_cc
F77=f77
CPP=

CEXT=c
FEXT=F   # which implies further action of cpp
OBJEXT=o
FOBJEXT=o
OBJFLAG=-c
FOBJFLAG=-c


###################################
# Choose your machine from here.
###################################

####################################
# Dec Alpha, OSF 1
AXPFLAGS= -unsigned  -non_shared  \
	-math_library=fast -prefix=all -reentrancy=none -assume=noaccuracy_sensitive \
	-unsigned_char -extern=strict_refdef -trapuv #  -D_INTRINSICS 
AXPLDFLAGS= -unsigned -assume=noaccuracy_sensitive -non_shared -D_INTRINSICS 
AXPOPTIM= -O -Ublas -float -Olimit 1000 # -cord  -feedback citcom.feedback
####################################

####################################
# CRAY Unicos systems
CRAYFLAGS=
CRAYLDFLAGS=
CRAYOPTIM=
####################################

####################################
#IBM AIX systems
AIXFLAGS= -D__aix__
AIXLDFLAGS=
AIXOPTIM= -O3 -qarch=pwr -qopt=pwr -Ublas
####################################

####################################
# SUNOS systems
SUNFLAGS= -D__sunos__ -Dconst=""
SUNLDFLAGS=
SUNOPTIM=-O3 -fsingle
####################################

####################################
# Solaris systems
SOLARISFLAGS= -D__solaris -Dconst="" -I/opt/mpi/include
SOLARISLDFLAGS=-fast -lsocket -lnsl -lthread
SOLARISOPTIM=-fast -xO4 -dalign -xtarget=ultra -xarch=v8plus -fsingle
####################################

####################################
#HP running HPUX
HPUXFLAGS=-Dconst=""
HPUXLDFLAGS=
HPUXOPTIM=+O3
####################################

####################################
# SGI with IRIX 
SGIFLAGS=
SGILDFLAGS=
SGIOPTIM=-O -fsingle
####################################

####################################
#PARAGON 
PARAGONFLAGS=-nx
PARAGONLDFLAGS=
PARAGONOPTIM=-O4
####################################

####################################
LinuxFLAGS=-I/usr/local/mpi/include -I/usr/include -I/usr/include/malloc
LinuxLDFLAGS=-no-pie
LinuxOPTIM=-O2



SILENCE_FLAG +=
FLAGS= $(LinuxFLAGS) -DCOMPRESS_BINARY=\"$(COMPRESS)\" $(SILENCE_FLAG) -g
LDFLAGS= $(LinuxLDFLAGS)
OPTIM= $(LinuxOPTIM)

#FLAGS= $(SOLARISFLAGS) -DCOMPRESS_BINARY=\"$(COMPRESS)\"
#LDFLAGS= $(SOLARISLDFLAGS)
#OPTIM= $(SOLARISOPTIM)

####################################################################
#	Citcom's files .....
####################################################################

CFILES= Advection_diffusion.c\
	Boundary_conditions.c\
	Citcom.c\
	Composition_adv.c\
	Construct_arrays.c\
	Convection.c\
	Drive_solvers.c\
	Element_calculations.c\
	Geometry_cartesian.c\
	General_matrix_functions.c\
	Stokes_flow_Incomp.c\
	Instructions.c\
	Nodal_mesh.c\
	Output.c\
	Pan_problem_misc_functions.c\
	Parsing.c\
	Phase_change.c\
	Process_buoyancy.c\
	Process_velocity.c\
	Profiling.c\
	Shape_functions.c\
	Size_does_matter.c\
	Solver_multigrid.c\
	Solver_conj_grad.c\
	Topo_gravity.c\
	Viscosity_structures.c

HEADER = element_definitions.h\
	 global_defs.h\
	 viscosity_descriptions.h\
	 advection.h


FFILES=#Blas_lapack_interfaces.F

# Note, if make on your machine does not support nesting of
# Macro definitions (CRAY), you'll have to do the substitution yourself
# at this point...

OBJFILES=$(CFILES:.c=.o) $(FFILES:.f=.o)
#
global_defs.h: viscosity_descriptions.h advection.h\
		Convection_variables.h
		
# The following entries can probably be automated from $CFILES etc

Advection_diffusion.o: $(HEADER) Advection_diffusion.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Advection_diffusion.c
#
Boundary_conditions.o: $(HEADER) Boundary_conditions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Boundary_conditions.c
#	
Citcom.o: $(HEADER) Citcom.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Citcom.c
#	
Composition_adv.o: $(HEADER) Composition_adv.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Composition_adv.c
#	
Construct_arrays.o: $(HEADER) Construct_arrays.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Construct_arrays.c
#	
Convection.o: $(HEADER) Convection.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Convection.c
#	
Drive_solvers.o: $(HEADER) Drive_solvers.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Drive_solvers.c
#	
Element_calculations.o: $(HEADER) Element_calculations.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Element_calculations.c
	
General_matrix_functions.o: $(HEADER) General_matrix_functions.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   General_matrix_functions.c
		
Geometry_cartesian.o: $(HEADER) Geometry_cartesian.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Geometry_cartesian.c
	
Instructions.o: $(HEADER) Instructions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Instructions.c
	
Nodal_mesh.o: $(HEADER) Nodal_mesh.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Nodal_mesh.c

Output.o: $(HEADER) Output.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Output.c
	
Pan_problem_misc_functions.o: $(HEADER)  Pan_problem_misc_functions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Pan_problem_misc_functions.c
		
Parsing.o: $(HEADER) Parsing.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Parsing.c

Phase_change.o: $(HEADER) Phase_change.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Phase_change.c

Process_velocity.o: $(HEADER) Process_velocity.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Process_velocity.c

Process_buoyancy.o: $(HEADER) Process_buoyancy.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Process_buoyancy.c
			
Profiling.o: $(HEADER) Profiling.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Profiling.c
	
Shape_functions.o: $(HEADER) Shape_functions.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Shape_functions.c
	
Size_does_matter.o: $(HEADER) Size_does_matter.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Size_does_matter.c
	
Solver_conj_grad.o: $(HEADER) Solver_conj_grad.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Solver_conj_grad.c
	
Solver_multigrid.o: $(HEADER) Solver_multigrid.c
	$(CC) $(OPTIM) $(FLAGS) $(OBJFLAG)   Solver_multigrid.c

Stokes_flow_Incomp.o: $(HEADER) Stokes_flow_Incomp.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Stokes_flow_Incomp.c
	
Topo_gravity.o: $(HEADER) Topo_gravity.c 
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Topo_gravity.c
	
Viscosity_structures.o: $(HEADER) Viscosity_structures.c
	$(CC) $(OPTIM) $(FLAGS)  $(OBJFLAG)  Viscosity_structures.c
	

Blas_lapack_interfaces.f: Blas_lapack_interfaces.F
	$(CPP) $(OPTIM) -P Blas_lapack_interfaces.F Blas_lapack_interfaces.f
	$(F77) $(OPTIM) $(FOPTIM) -c Blas_lapack_interfaces.f
	

smaller: 
	compress $(CFILES)

larger:
	uncompress $(CFILES)


citcom.x: $(OBJFILES)  $(HEADER) Makefile
	$(CC) $(OPTIM) $(FLAGS) $(LDFLAGS)  -o citcom.x $(OBJFILES)  $(FFTLIB)  $(LIB)

clean:
	rm -f *.o *.x *pixie citcom2 *.f
	@echo "All clean!"

# Profiling/testing versions 

citcom2: $(OBJFILES)  $(HEADER) Makefile
	$(CC) -p $(OPTIM) $(LDFLAGS)  -o citcom2  $(OBJFILES)  $(FFTLIB)  $(LIB)
	

citcom2.pixie: citcom2 # $(OBJFILES)  $(HEADER) Makefile
	pixie -nopids -o citcom2.pixie citcom2

pixie: citcom2.pixie
	echo "Pixie profiling code for DEC alpha in citcom2.pixie"
