
bin_PROGRAMS = velocityf

velocityf_SOURCES = velocityf.f90
velocityf_FCFLAGS = -I../../../../../include
velocityf_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
velocityf_LDADD += -lmpi
endif

EXTRA_DIST = README \
	plot_credible.py \
	plot_mean.py \
	plot_misfit.py \
	plot_partition_count_histogram.py \
	plot_partition_x_hist.py \
	velocityf.py