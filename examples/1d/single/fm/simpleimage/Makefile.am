
bin_PROGRAMS = simpleimage \
	mkimage

simpleimage_SOURCES = simpleimage.c image.c image.h
simpleimage_CPPFLAGS = -I../../../../../include
simpleimage_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
simpleimage_LDADD += -lmpi
endif

mkimage_SOURCES = mkimage.c image.c image.h
mkimage_CPPFLAGS = -I../../../../../include
mkimage_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
mkimage_LDADD += -lmpi
endif

EXTRADIST = README \
	plot_amplitude.py \
	plot_background.py \
	plot_misfit.py \
	plot_noise.py \
	plot_sigma.py \
	plot_sigmax.py \
	plot_sigmay.py \
	plot_theta.py \
	plot_x0.py \
	plot_y0.py \
	simpleimage.py

