#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([celerium], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_RANLIB

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Set default flags
: ${CXXFLAGS="-O3"}
: ${CFLAGS="-O3"}

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MKDIR_P

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
		 tests/ArithmeticVector/Makefile
                 tests/WFSolver/Makefile
		 tests/Potential/Makefile
		 tests/PotentialAndWFSolver/Makefile
		 tests/Interpolator/Makefile
		 tests/SphericalHarmonics/Makefile
		 tests/OrbitalClass/Makefile
		 tests/Element/Makefile
		 tests/ElementaryCell/Makefile
		 tests/CubaWrapper/Makefile
		 tests/PeriodicLowdin/Makefile])

AC_OUTPUT
