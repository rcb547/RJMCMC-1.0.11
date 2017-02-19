
bin_PROGRAMS = zeromultistep

zeromultistep_SOURCES = zeromultistep.c
zeromultistep_CPPFLAGS = -I../../../../../include
zeromultistep_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
zeromultistep_LDADD += -lmpi
endif

EXTRA_DIST = mkdata.py \
	data.txt \
	README \
	plot_credible.py \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_mode.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py \
	zeromultistep.py


