
bin_PROGRAMS = step

step_SOURCES = step.c
step_CPPFLAGS = -I../../../../../include
step_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
step_LDADD += -lmpi
endif

EXTRA_DIST = README \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py \
	step.py

