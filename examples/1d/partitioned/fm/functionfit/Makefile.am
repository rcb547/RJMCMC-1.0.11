
bin_PROGRAMS = functionfit

functionfit_SOURCES = functionfit.c
functionfit_CPPFLAGS = -I../../../../../include
functionfit_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
functionfit_LDADD += -lmpi
endif

EXTRA_DIST = README \
	functionfit.py \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py
