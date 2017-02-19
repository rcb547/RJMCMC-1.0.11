
lib_LTLIBRARIES = librjmcmc.la

librjmcmc_la_SOURCES = \
	bbox2d.c \
	curvefit.c \
	dataset1d.c \
	dataset2d.c \
	delaunay2d.c \
	engine.c \
	forwardmodel_f.c \
	forwardmodel_part1d.c \
	forwardmodel_part1d_hierarchical.c \
	forwardmodel_part2d.c \
	forwardmodel_part2d_hierarchical.c \
	forwardmodel_single.c \
	forwardmodel_util.c \
	forwardmodelparameter.c \
	part1d_forwardmodel.c \
	part1d_natural_rj.c \
	part1d_regression_rj.c \
	part1d_zero.c \
	part2d_forwardmodel.c \
	part2d_regression_rj.c \
	position_map1d.c \
	position_map2d.c \
	position_map2d_delaunay.c \
	position_map2d_linear.c \
	position_map2d_quadtree.c \
	quadtree.c \
	regression.c \
	regression_part1d.c \
	regression_part1d_natural.c \
	regression_part1d_zero.c \
	regression_part2d.c \
	resultset1d.c \
	resultset1dfm.c \
	resultset2d.c \
	resultset2dfm.c \
	resultsetfm.c \
	rjmcmc_debug.c \
	rjmcmc_random.c \
	rjmcmc_util.c \
	single1d_regression.c \
	wellrng.c

librjmcmc_la_CPPFLAGS = $(AM_CPPFLAGS) -I../include


if WITHOPENMPI
librjmcmc_la_LIBADD = -lmpi
endif


