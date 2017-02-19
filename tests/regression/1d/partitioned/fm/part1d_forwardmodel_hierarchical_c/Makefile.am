
bin_PROGRAMS = part1d_forwardmodel_hierarchical_c

BASEDIR=../../../../../..

part1d_forwardmodel_hierarchical_c_SOURCES = part1d_forwardmodel_hierarchical_c.c
part1d_forwardmodel_hierarchical_c_hierarchical_cPPFLAGS = -I$(BASEDIR)/include
part1d_forwardmodel_hierarchical_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
part1d_forwardmodel_hierarchical_c_LDADD += -lmpi
endif
