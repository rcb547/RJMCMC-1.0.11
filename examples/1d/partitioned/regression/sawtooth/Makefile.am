
bin_PROGRAMS = sawtooth

sawtooth_SOURCES = sawtooth.c
sawtooth_CPPFLAGS = -I../../../../../include
sawtooth_LDADD = ../../../../../src/.libs/librjmcmc.a 

if WITHOPENMPI
sawtooth_LDADD += -lmpi
endif

EXTRA_DIST = mkdata.py \
	data.txt \
	README \
	plot_credible.py \
	plot_mean.py \
	plot_median.py \
	plot_misfit.py \
	plot_mode.py

