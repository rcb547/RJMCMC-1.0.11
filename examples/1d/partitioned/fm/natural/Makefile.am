
bin_PROGRAMS = natural

natural_SOURCES = natural.c
natural_CPPFLAGS = -I../../../../../include
natural_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
natural_LDADD += -lmpi
endif

EXTRA_DIST = data.txt \
	mkdata.py \
	natural.py \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py
