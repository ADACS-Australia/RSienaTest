# AC_ARG_WITH([mpi-lib-path],
#   AS_HELP_STRING([--with-mpi-lib-path=DIR],
#                  [location of MPI library files]),
#   [
#     if test "${withval}" != "no" ; then
#       MPI_LIB_PATH="${withval}"
#     fi
#   ])

# After changing this file run `autoconf` to produce the 'configure' script.
#
# MPI discovery options for the build process:
#
# - Run with automatic detection: `R CMD INSTALL {pkg}`
#
# - Run with a prefix folder.  Subfolders will be searched for binarys.
#   `R CMD INSTALL {pkg} --configure-args=--with-mpi=/path/to/location`
#
# - Fully specify paths and type.
#   `R CMD INSTALL {pkg} --configure-args="        \\
#    --with-mpi-include=/path/to/mpi_include_dir  \\
#    --with-mpi-libpath=/path/to/mpi_lib_dir      \\
#    --with-mpi-type={OPENMPI,MPICH,MPICH2,LAM}"`

AC_PREREQ([2.69])
AC_INIT([RSienaTest], 1.1-289)
AC_CONFIG_SRCDIR([src])

# Find R home and set CC/CFLAGS.  This is required for the custom openmp test.
R_HOME=`R RHOME`
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi
RBIN="${R_HOME}/bin/R"
# CC=`"${RBIN}" CMD config CC`;
# CFLAGS=`"${RBIN}" CMD config CFLAGS`
# LIBS="${PKG_LIBS}"

AC_PROG_CXX
## commented out since was not structured anyway
## (is a full check of depencies neccessary?)
# AC_PROG_CC
# AC_CHECK_HEADERS([sys/time.h unistd.h])
# AC_CHECK_HEADER_STDBOOL
# AC_C_INLINE
# AC_C_RESTRICT
# AC_TYPE_SIZE_T
# AC_CHECK_TYPES([ptrdiff_t])
# AC_FUNC_ERROR_AT_LINE
# AC_FUNC_MALLOC
# AC_HEADER_MAJOR
# AC_FUNC_REALLOC
# AC_CHECK_FUNCS([clock_gettime gettimeofday memmove memset pow select sqrt])


## MPI
AC_CHECK_PROG(HAS_MPICC,mpicc,yes)
if test "$HAS_MPICC" = yes ; then
  PKG_CPPFLAGS="$MPI_DEFS `mpic++ -showme:compile`"
  PKG_LIBS="`mpic++ -showme:link`"
  AC_LANG_CONFTEST(AC_LANG_SOURCE([
    #include <mpi.h>
    #include <R_ext/Rdynload.h>
    void
    R_init_conftest(DllInfo *info)
    { MPI_Init(NULL, NULL); }
  ]))
  acr_check_dynlib_result=no
  AC_MSG_CHECKING([check if MPI can be loaded by R])
  PKG_CFLAGS="$PKG_CPPFLAGS" PKG_LIBS="$PKG_LIBS" \
      "$RBIN" CMD SHLIB conftest.c \
      1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
  "$RBIN" --slave -e \
      "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" \
      1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD \
    && acr_check_dynlib_result=yes
  AC_MSG_RESULT($acr_check_dynlib_result)
  if test "$acr_check_dynlib_result" = no ; then
    AC_MSG_RESULT([will compile without MPI support])
    PKG_CPPFLAGS=
    PKG_LIBS=
  else
    PKG_CPPFLAGS="-DMPI2 $PKG_CPPFLAGS"
  fi
fi

if test "`uname`" = Darwin ; then
  PKG_CPPFLAGS="-DMAC $PKG_CPPFLAGS"
fi
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

## OpenMP
AC_OPENMP
AC_LANG_CONFTEST(AC_LANG_SOURCE([
  #include <omp.h>
  #include <R_ext/Rdynload.h>
  void
  R_init_conftest(DllInfo *info)
  { omp_get_num_threads(); }
]))
acr_check_dynlib_result=no
AC_MSG_CHECKING([check if OpenMP can be loaded by R])
PKG_CFLAGS="$PKG_CPPFLAGS" PKG_LIBS="$PKG_LIBS" \
    "$RBIN" CMD SHLIB conftest.c \
    1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
"$RBIN" --slave -e \
    "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" \
    1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD \
  && acr_check_dynlib_result=yes
AC_MSG_RESULT($acr_check_dynlib_result)
if test "$acr_check_dynlib_result" = no; then
  AC_MSG_WARN([OpenMP does NOT appear to be supported on this system.])
else
  OPENMP_CFLAGS="$OPENMP_CFLAGS -DHAVE_OPENMP"
  AC_SUBST(OPENMP_CFLAGS)
fi

## substitute externalized source list
PKG_SOURCES=$(cat src/sources.list)
AC_SUBST(PKG_SOURCES)

AC_OUTPUT(src/Makevars)
