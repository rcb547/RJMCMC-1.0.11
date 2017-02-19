
bin_PROGRAMS = regression

regression_SOURCES = regression.c
regression_CPPFLAGS = -I../../../../../include
regression_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
regression_LDADD += -lmpi
endif

EXTRA_DIST = regression.py \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py