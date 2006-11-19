dnl Process this file with autoconf to produce a configure script.
AC_INIT(sgEvolver, 0.5.0)

dnl Get the target and build system types and add appropriate options
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([no-define])

PKG_CHECK_MODULES(DEPS, libMems-1.3 >= 1.4.0)
AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_LIBS)

AC_PREFIX_DEFAULT(/usr/local)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S 

dnl Checks for header files.
AC_HEADER_STDC

dnl Check for getopt_long
AC_CHECK_FUNC(getopt_long,
[GETOPT_LONG_SYSTEM=""],
[GETOPT_LONG_SYSTEM='getopt.$(OBJEXT) getopt1.$(OBJEXT)']
)
AC_SUBST([GETOPT_LONG_SYSTEM])

dnl Check for getopt
AC_CHECK_FUNC(getopt,
[HAVE_GETOPT="true"],
AC_DEFINE(GETOPT_UNDEFINED,[],"Define this if the system does not provide getopt
()")
)

dnl Mac OS X won't allow static compilation...
STATIC_FLAG="-static"
if ( test "x$target_vendor" = "xapple") then
        STATIC_FLAG=""
fi
AC_SUBST(STATIC_FLAG)


dnl Allow debugging compilation
AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

AM_CONFIG_HEADER(config.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
dnl AC_C_BIGENDIAN
AC_HEADER_TIME

dnl Checks for library functions.
AC_PROG_GCC_TRADITIONAL

dnl SAVE_LIBRARY_VERSION
AC_SUBST(LIBTOOL_VERSION_INFO)


AC_OUTPUT(Makefile src/Makefile )
