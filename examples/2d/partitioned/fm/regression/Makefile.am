
bin_PROGRAMS = regression

regression_SOURCES = regression.c
regression_CPPFLAGS = -I../../../../../include
regression_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
regression_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	plot.py \
	plot_mean.py \
	plot_median.py \
	plot_partition_count_histogram.py \
	data.txt
