
bin_PROGRAMS = ausmoho

ausmoho_SOURCES = ausmoho.c
ausmoho_CPPFLAGS = -I../../../../../include
ausmoho_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
ausmoho_LDADD += -lmpi
endif

EXTRA_DIST = README \
	data.txt \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot.py


