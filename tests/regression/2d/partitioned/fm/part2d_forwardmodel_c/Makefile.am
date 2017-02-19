
bin_PROGRAMS = part2d_forwardmodel_c

BASEDIR=../../../../../..

part2d_forwardmodel_c_SOURCES = part2d_forwardmodel_c.c
part2d_forwardmodel_c_CPPFLAGS = -I$(BASEDIR)/include
part2d_forwardmodel_c_LDADD = $(BASEDIR)/src/.libs/librjmcmc.a

if WITHOPENMPI
part2d_forwardmodel_c_LDADD += -lmpi
endif
