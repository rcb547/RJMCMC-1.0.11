
bin_PROGRAMS = spherefit

spherefit_SOURCES = spherefit.c
spherefit_CPPFLAGS = -I../../../../../include
spherefit_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
spherefit_LDADD += -lmpi
endif

EXTRADIST = README \
	plot_misfit.py \
	plot_r.py \
	plot_x.py \
	plot_y.py \
	plot_z.py \
	spherefit.py

