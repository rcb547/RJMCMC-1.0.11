
bin_PROGRAMS = MPI_part1d_forwardmodel_natural_c

BASEDIR=../../../../../..

MPI_part1d_forwardmodel_natural_c_SOURCES = MPI_part1d_forwardmodel_natural_c.c
MPI_part1d_forwardmodel_natural_c_CPPFLAGS = -I$(BASEDIR)/include
MPI_part1d_forwardmodel_natural_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
MPI_part1d_forwardmodel_natural_c_LDADD += -lmpi
endif
