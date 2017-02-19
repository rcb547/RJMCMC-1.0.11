
TESTS_CPPFLAGS = -I../include
TESTS_LDADD = ../src/.libs/librjmcmc.a \
	-lcunit

if WITHOPENMPI
TESTS_LDADD += -lmpi
endif

bin_PROGRAMS = 	\
	curvefit_tests \
	dataset1d_tests \
	delaunay2d_tests \
	part1d_forwardmodel_tests \
	part1d_regression_tests \
	part2d_forwardmodel_tests \
	position_map_tests \
	position_map1d_thrash  \
	position_map2d_tests \
	profile_position_map2d \
	quadtree_tests \
	random_tests \
	resultset2dfm_tests \
	util_tests


curvefit_tests_SOURCES = curvefit_tests.c
curvefit_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
curvefit_tests_LDADD = $(TESTS_LDADD)

dataset1d_tests_SOURCES = dataset1d_tests.c
dataset1d_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
dataset1d_tests_LDADD = $(TESTS_LDADD)

delaunay2d_tests_SOURCES = delaunay2d_tests.c
delaunay2d_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
delaunay2d_tests_LDADD = $(TESTS_LDADD)

part1d_forwardmodel_tests_SOURCES = part1d_forwardmodel_tests.c
part1d_forwardmodel_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
part1d_forwardmodel_tests_LDADD = $(TESTS_LDADD)

part1d_regression_tests_SOURCES = part1d_regression_tests.c
part1d_regression_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
part1d_regression_tests_LDADD = $(TESTS_LDADD)

part2d_forwardmodel_tests_SOURCES = part2d_forwardmodel_tests.c
part2d_forwardmodel_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
part2d_forwardmodel_tests_LDADD = $(TESTS_LDADD)

position_map_tests_SOURCES = position_map_tests.c
position_map_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
position_map_tests_LDADD = $(TESTS_LDADD)

position_map1d_thrash = position_map1d_thrash.c
position_map1d_thrash_CPPFLAGS = $(TESTS_CPPFLAGS)
position_map1d_thrash_LDADD = $(TESTS_LDADD)

position_map2d_tests_SOURCES = position_map2d_tests.c
position_map2d_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
position_map2d_tests_LDADD = $(TESTS_LDADD)

profile_position_map2d_SOURCES = profile_position_map2d.c
profile_position_map2d_CPPFLAGS = $(TESTS_CPPFLAGS)
profile_position_map2d_LDADD = $(TESTS_LDADD)

quadtree_tests_SOURCES = quadtree_tests.c
quadtree_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
quadtree_tests_LDADD = $(TESTS_LDADD)

random_tests_SOURCES = random_tests.c
random_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
random_tests_LDADD = $(TESTS_LDADD)

resultset2dfm_tests_SOURCES = resultset2dfm_tests.c
resultset2dfm_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
resultset2dfm_tests_LDADD = $(TESTS_LDADD)

util_tests_SOURCES = util_tests.c
util_tests_CPPFLAGS = $(TESTS_CPPFLAGS)
util_tests_LDADD = $(TESTS_LDADD)


