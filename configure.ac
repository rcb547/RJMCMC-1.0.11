AC_PREREQ([2.13])
AC_INIT([rjmcmc], [1.0.11], [http://www.iearth.org.au/codes/RJMCMC], [RJMCMC])
AM_INIT_AUTOMAKE([1.10 no-define tar-ustar])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([include/rjmcmc/rjmcmc_config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_PROG_FC
AC_PROG_LIBTOOL

AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB(m, sqrt)

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([assert.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt log])

#
# Debugging/Devel compile flags
#
AC_ARG_ENABLE(
	debug,
	AS_HELP_STRING(--enable-debug, enable debug and warning gcc flags),
	[if test x$enableval = xyes; then
	    CFLAGS="-Wall -g"
	 fi])

#
# Open MPI
#

AC_ARG_WITH([openmpi-include-path],
  [AS_HELP_STRING([--with-openmpi-include-path],
    [location of the Openmpi headers, defaults to /usr/include/openmpi])],
  [OPENMPI_CPPFLAGS="-I$withval"],
  [OPENMPI_CPPFLAGS='-I/usr/include/openmpi'])

AC_SUBST([OPENMPI_CPPFLAGS])
CPPFLAGS="$CPPFLAGS $OPENMPI_CPPFLAGS"

AC_ARG_WITH([openmpi-lib-path],
  [AS_HELP_STRING([--with-openmpi-lib-path], [location of the Openmpi libraries, defaults to /usr/lib/openmpi])],
  [OPENMPI_LDFLAGS="-L$withval"],
  [OPENMPI_LDFLAGS='-L/usr/lib/openmpi'])

AC_SUBST([OPENMPI_LDFLAGS])
LDFLAGS="$LDFLAGS $OPENMPI_LDFLAGS"

AC_ARG_WITH(
	[openmpi],
	AS_HELP_STRING(--with-openmpi, enable openmpi support),
	[],
	[with_openmpi=check])

AS_IF([test "x$with_openmpi" != xno],
      [

        AC_CHECK_HEADERS([mpi.h],
                         [],
                         [if test "x$with_openmpi" != xcheck; then 
                            AC_MSG_FAILURE([--with-openmpi was given but mpi.h unavailable. Try setting CPPFLAGS="-I<path to mpi.h>".])
                          else
                            withopenmpi=no
                          fi])

        AC_CHECK_LIB([mpi], 
                     [MPI_Init],
	             [AC_SUBST([LIBMPI], ["-lmpi"])
                      AC_DEFINE([HAVE_LIBMPI], [1],
                      [Define if you have libmpi])],
 	             [if test "x$with_openmpi" != xcheck; then
                        AC_MSG_FAILURE([--with-openmpi was given, but test for openmpi failed])
                      else
                        withopenmpi=no
                      fi])
      ])
	    
AM_CONDITIONAL([WITHOPENMPI], [test "x$withopenmpi" != xno])

AC_CONFIG_FILES([Makefile
		 rjmcmc.pc
                 src/Makefile
		 include/Makefile
                 examples/1d/partitioned/fm/functionfit/Makefile
		 examples/1d/partitioned/fm/functionfitf/Makefile
		 examples/1d/partitioned/fm/gaussian/Makefile
		 examples/1d/partitioned/fm/natural/Makefile
                 examples/1d/partitioned/fm/regression/Makefile
		 examples/1d/partitioned/fm/step/Makefile
		 examples/1d/partitioned/fm/synthdata/Makefile
		 examples/1d/partitioned/fm/velocityf/Makefile
                 examples/1d/partitioned/regression/multiquad/Makefile
		 examples/1d/partitioned/regression/multistep/Makefile
		 examples/1d/partitioned/regression/sawtooth/Makefile
		 examples/1d/partitioned/regression/zeromultistep/Makefile
		 examples/1d/single/fm/gaussian/Makefile
		 examples/1d/single/fm/simplef/Makefile
		 examples/1d/single/fm/simpleimage/Makefile
		 examples/1d/single/fm/spherefit/Makefile
		 examples/1d/single/regression/cubic/Makefile
		 examples/2d/partitioned/fm/mpi_regression/Makefile
		 examples/2d/partitioned/fm/regression/Makefile
		 examples/2d/partitioned/fm/gaussian/Makefile
		 examples/2d/partitioned/regression/ausmoho/Makefile
		 examples/2d/partitioned/regression/callbackex/Makefile
		 examples/2d/partitioned/regression/disc/Makefile
		 examples/2d/partitioned/regression/gaussian/Makefile
		 examples/2d/partitioned/regression/square/Makefile
		 tests/Makefile
		 tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_c/Makefile
		 tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_hierarchical_c/Makefile
		 tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_natural_c/Makefile
		 tests/regression/1d/partitioned/fm/MPI_part1d_forwardmodel_natural_hierarchical_c/Makefile
		 tests/regression/1d/partitioned/fm/part1d_forwardmodel_c/Makefile
		 tests/regression/1d/partitioned/fm/part1d_forwardmodel_hierarchical_c/Makefile
		 tests/regression/1d/partitioned/fm/part1d_forwardmodel_natural_c/Makefile
		 tests/regression/1d/partitioned/fm/part1d_forwardmodel_natural_hierarchical_c/Makefile
		 tests/regression/1d/partitioned/fm/part1d_forwardmodel_zero_cubic_c/Makefile
		 tests/regression/1d/single/fm/MPI_single_forwardmodel_c/Makefile
		 tests/regression/1d/single/fm/MPI_single_forwardmodel_hierarchical_c/Makefile
		 tests/regression/1d/single/fm/single_forwardmodel_c/Makefile
		 tests/regression/1d/single/fm/single_forwardmodel_hierarchical_c/Makefile
		 tests/regression/2d/partitioned/fm/MPI_part2d_forwardmodel_c/Makefile
		 tests/regression/2d/partitioned/fm/MPI_part2d_forwardmodel_hierarchical_c/Makefile
		 tests/regression/2d/partitioned/fm/part2d_forwardmodel_c/Makefile
		 tests/regression/2d/partitioned/fm/part2d_forwardmodel_hierarchical_c/Makefile])

#
# Old/Unmaintained examples
#
#                 examples/fm/part1d/bore/Makefile
#		 examples/mpi/fm/spherefit/Makefile
#		 examples/mpi/single/cubic/Makefile
#		 examples/part1d/giss/Makefile
#		 examples/part1d/msl/Makefile
#		 examples/part1d/knot/fixed/Makefile

AC_OUTPUT
