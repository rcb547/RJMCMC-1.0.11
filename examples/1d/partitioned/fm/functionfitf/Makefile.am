
bin_PROGRAMS = functionfitf

functionfitf_SOURCES = functionfitf.f90
functionfitf_FCFLAGS = -I../../../../../include
functionfitf_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
functionfitf_LDADD += -lmpi
endif

EXTRA_DIST = README \
	functionfitf.py \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py

