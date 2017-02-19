
bin_PROGRAMS = multiquad

multiquad_SOURCES = multiquad.c
multiquad_CPPFLAGS = -I../../../../../include
multiquad_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
multiquad_LDADD += -lmpi
endif

EXTRA_DIST = README \
	mkdata.py \
	data.txt \
	multiquad.py \
	plot_credible.py \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_mode.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py
