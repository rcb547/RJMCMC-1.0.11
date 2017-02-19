
bin_PROGRAMS = simplef

simplef_SOURCES = simplef.f90
simplef_FCFLAGS = -I../../../../../include
simplef_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
simplef_LDADD += -lmpi
endif

EXTRA_DIST = README \
	plot_misfit.py \
	plot_history.py
