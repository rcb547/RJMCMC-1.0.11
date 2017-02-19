
bin_PROGRAMS = multistep

multistep_SOURCES = multistep.c
multistep_CPPFLAGS = -I../../../../../include
multistep_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
multistep_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	data.txt \
	multistep.py \
	plot_credible.py \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_mode.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py
