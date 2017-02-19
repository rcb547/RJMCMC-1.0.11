
bin_PROGRAMS = single_forwardmodel_hierarchical_c

BASEDIR=../../../../../..

single_forwardmodel_hierarchical_c_SOURCES = single_forwardmodel_hierarchical_c.c
single_forwardmodel_hierarchical_c_CPPFLAGS = -I$(BASEDIR)/include
single_forwardmodel_hierarchical_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
single_forwardmodel_hierarchical_c_LDADD += -lmpi
endif
