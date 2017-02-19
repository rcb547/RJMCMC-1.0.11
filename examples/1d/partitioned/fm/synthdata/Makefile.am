
bin_PROGRAMS = synthdata

synthdata_SOURCES = synthdata.c
synthdata_CPPFLAGS = -I../../../../../include
synthdata_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
synthdata_LDADD += -lmpi
endif

EXTRA_DIST = data.txt \
	README \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py \
	synthdata.py
