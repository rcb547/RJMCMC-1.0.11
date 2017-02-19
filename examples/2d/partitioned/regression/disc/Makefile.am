
bin_PROGRAMS = disc

disc_SOURCES = disc.c
disc_CPPFLAGS = -I../../../../../include
disc_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
disc_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	data.txt \
	plot_mean.py \
	plot_median.py \
	plot_partition_count_histogram.py \
	plot.py


