
bin_PROGRAMS = gaussian

gaussian_SOURCES = gaussian.c
gaussian_CPPFLAGS = -I../../../../../include
gaussian_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
gaussian_LDADD += -lmpi
endif

EXTRA_DIST = README \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py \
	gaussian.py

