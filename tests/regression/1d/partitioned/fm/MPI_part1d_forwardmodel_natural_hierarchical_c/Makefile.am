
bin_PROGRAMS = MPI_part1d_forwardmodel_natural_hierarchical_c

BASEDIR=../../../../../..

MPI_part1d_forwardmodel_natural_hierarchical_c_SOURCES = MPI_part1d_forwardmodel_natural_hierarchical_c.c
MPI_part1d_forwardmodel_natural_hierarchical_c_CPPFLAGS = -I$(BASEDIR)/include
MPI_part1d_forwardmodel_natural_hierarchical_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
MPI_part1d_forwardmodel_natural_hierarchical_c_LDADD += -lmpi
endif
