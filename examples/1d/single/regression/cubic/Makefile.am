
bin_PROGRAMS = cubic

cubic_SOURCES = cubic.c
cubic_CPPFLAGS = -I../../../../../include
cubic_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
cubic_LDADD += -lmpi
endif

EXTRA_DIST = README \
	cubic.py \
	data.txt \
	mkdata.py \
	plot_credible.py \
	plot_lambda.py \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_mode.py \
	plot_order.py

