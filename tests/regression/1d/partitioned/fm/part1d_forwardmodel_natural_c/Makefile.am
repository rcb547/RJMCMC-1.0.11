
bin_PROGRAMS = part1d_forwardmodel_natural_c

BASEDIR=../../../../../..

part1d_forwardmodel_natural_c_SOURCES = part1d_forwardmodel_natural_c.c
part1d_forwardmodel_natural_c_CPPFLAGS = -I$(BASEDIR)/include
part1d_forwardmodel_natural_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
part1d_forwardmodel_natural_c_LDADD += -lmpi
endif
