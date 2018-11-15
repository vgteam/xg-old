#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([xg], [1.0.0], [https://github.com/ekg], [xg], [https://github.com/vgteam/xg])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/main.cpp])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_LANG([C++])
AC_PROG_RANLIB
AM_PROG_AR

AX_CXX_COMPILE_STDCXX_11([ext], [mandatory])
AC_OPENMP([], [AC_MSG_ERROR([Unable to find OpenMP support])])
AX_PTHREAD([], [AC_MSG_ERROR([Unable to find POSIX thread support])])

# Checks for libraries.
PKG_CHECK_MODULES([ZLIB], [zlib >= 1.2.8], [], [AC_MSG_ERROR([Unable to find zlib])])
AC_CHECK_LIB([z], [compress])
PKG_CHECK_MODULES([PROTOBUF], [protobuf >= 3.1.0], [], [AC_MSG_ERROR([Unable to find protobuf])])
PKG_CHECK_MODULES([SDSL], [sdsl-lite >= 2.1.0], [],
                  [AC_MSG_ERROR([Unable to find sdsl-lite headers or library])])
AC_CHECK_LIB([sdsl], [main])
PKG_CHECK_MODULES([DIVSUFSORT], [libdivsufsort >= 2.0.1], [], [AC_MSG_ERROR([Unable to find libdivsufsort])])
PKG_CHECK_MODULES([DIVSUFSORT64], [libdivsufsort64 >= 2.0.1], [],
             [AC_MSG_ERROR([Unable to find libdivsufsort64])])
AC_CHECK_LIB([divsufsort], [divsufsort])
AC_CHECK_LIB([divsufsort64], [divsufsort])
PKG_CHECK_MODULES([LIBSPARSEHASH], [libsparsehash >= 2.0.2], [],
                  [AC_MSG_ERROR([Unable to find gcsa2 headers or library])])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile xg.pc])
AC_OUTPUT
