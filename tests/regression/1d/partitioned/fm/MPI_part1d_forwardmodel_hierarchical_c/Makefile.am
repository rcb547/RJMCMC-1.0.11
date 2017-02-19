
bin_PROGRAMS = MPI_part1d_forwardmodel_hierarchical_c

BASEDIR=../../../../../..

MPI_part1d_forwardmodel_hierarchical_c_SOURCES = MPI_part1d_forwardmodel_hierarchical_c.c
MPI_part1d_forwardmodel_hierarchical_c_hierarchical_cPPFLAGS = -I$(BASEDIR)/include
MPI_part1d_forwardmodel_hierarchical_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
MPI_part1d_forwardmodel_hierarchical_c_LDADD += -lmpi
endif
