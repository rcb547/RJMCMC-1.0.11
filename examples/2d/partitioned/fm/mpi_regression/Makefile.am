
if WITHOPENMPI

bin_PROGRAMS = mpi_regression \
	mpi_regression_checkpoint \
	mpi_regression_restart

mpi_regression_SOURCES = mpi_regression.c
mpi_regression_CPPFLAGS = -I../../../../../include
mpi_regression_LDADD = ../../../../../src/.libs/librjmcmc.a

mpi_regression_LDADD += -lmpi

mpi_regression_checkpoint_SOURCES = mpi_regression_checkpoint.c
mpi_regression_checkpoint_CPPFLAGS = -I../../../../../include
mpi_regression_checkpoint_LDADD = ../../../../../src/.libs/librjmcmc.a

mpi_regression_checkpoint_LDADD += -lmpi

mpi_regression_restart_SOURCES = mpi_regression_restart.c
mpi_regression_restart_CPPFLAGS = -I../../../../../include
mpi_regression_restart_LDADD = ../../../../../src/.libs/librjmcmc.a

mpi_regression_restart_LDADD += -lmpi

endif

EXTRA_DIST = README \
	mkdata.py \
	plot.py \
	plot_mean.py \
	plot_median.py \
	plot_partition_count_histogram.py \
	data.txt
