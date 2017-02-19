
bin_PROGRAMS = MPI_part2d_forwardmodel_c

BASEDIR=../../../../../..

MPI_part2d_forwardmodel_c_SOURCES = MPI_part2d_forwardmodel_c.c
MPI_part2d_forwardmodel_c_CPPFLAGS = -I$(BASEDIR)/include
MPI_part2d_forwardmodel_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
MPI_part2d_forwardmodel_c_LDADD += -lmpi
endif
