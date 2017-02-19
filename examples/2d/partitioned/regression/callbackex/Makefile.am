
bin_PROGRAMS = callbackex

callbackex_SOURCES = callbackex.c
callbackex_CPPFLAGS = -I../../../../../include
callbackex_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
callbackex_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	data.txt \
	plot_mean.py \
	plot_median.py \
	plot_partition_count_histogram.py \
	plot.py

