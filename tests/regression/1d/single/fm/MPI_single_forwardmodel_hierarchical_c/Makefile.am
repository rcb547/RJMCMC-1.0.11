
bin_PROGRAMS = MPI_single_forwardmodel_hierarchical_c

BASEDIR=../../../../../..

MPI_single_forwardmodel_hierarchical_c_SOURCES = MPI_single_forwardmodel_hierarchical_c.c
MPI_single_forwardmodel_hierarchical_c_CPPFLAGS = -I$(BASEDIR)/include
MPI_single_forwardmodel_hierarchical_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
MPI_single_forwardmodel_hierarchical_c_LDADD += -lmpi
endif
