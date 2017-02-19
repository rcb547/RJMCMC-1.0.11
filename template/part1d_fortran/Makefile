
#
# If you installed rjmcmc to a different location you may have
# to modify these paths
#
INCLUDES = -I/usr/local/include
LIBPATH = -L/usr/local/lib

#
# If mpi was enabled, you may need an extra include path for mpi.h
#
INCLUDES += -I/usr/include/openmpi

#
# If you're using GSL, you'll need to include GSL include paths
#
#INCLUDES += $(shell gsl-config --cflags)

FC = gfortran
FFLAGS = -c -g -Wall $(INCLUDES)

LIBS = -lrjmcmc

#
# If mpi was enabled, you will need to include the mpi libraries
#
#LIBS += -lmpi

#
# If you're using GSL, you'll need to include GSL libraries
#
#LIBS += $(shell gsl-config --libs)

OBJS = part1d_fortran.o

TARGETS = part1d_fortran

all : $(TARGETS)

part1d_fortran : $(OBJS)
	$(FC) -o part1d_fortran $(OBJS) $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -o $*.o $*.f90

clean :
	rm -f $(TARGETS) $(OBJS) 
