#-*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([pbc], [0.5.12], [blynn@cs.stanford.edu])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([./])
AC_PROG_LIBTOOL
#AC_CANONICAL_HOST

CFLAGS=
default_fink_path=/sw
case $host_os in
     darwin*) 
     	      dnl nested function support
     	      AC_MSG_NOTICE([Adding -fnested-functions to CFLAGS])
	      CFLAGS="-fnested-functions"

	      dnl fink installation
     	      AC_MSG_CHECKING([for a fink installation at $default_fink_path])
     	      if test -d $default_fink_path; then 
	      	 AC_MSG_RESULT([found it!])
		 AC_MSG_NOTICE([Adding -I$default_fink_path/include to CPPFLAGS])
	      	 CPPFLAGS="-I$default_fink_path/include $CPPFLAGS"
		 AC_MSG_NOTICE([Adding -L$default_fink_path/lib to LDFLAGS])
		 LDFLAGS="-L$default_fink_path/lib $LDFLAGS"
	      else
		AC_MSG_RESULT(none)
		AC_MSG_NOTICE([You may need to add set CPPFLAGS and LDFLAGS for gmp, etc.])
	      fi
	      ;;
esac

############################
# Configs for Windows DLLs.
# Framework for the below was extracted and
# modeled after the libgmp configure script.

AC_LIBTOOL_WIN32_DLL
AC_SUBST(LIBPBC_DLL,0)

case $host in
  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
    if test -z "$enable_shared"; then
      enable_shared=no
    fi
    # Don't allow both static and DLL.
    if test "$enable_shared" != no && test "$enable_static" != no; then
      AC_MSG_ERROR([cannot build both static and DLL, since gmp.h is different for each.
Use "--disable-static --enable-shared" to build just a DLL.])
    fi

    # "-no-undefined" is required when building a DLL, see documentation on
    # AC_LIBTOOL_WIN32_DLL.  Also, -no-undefined needs a version number
    # or it will complain about not having a nonnegative integer.
    if test "$enable_shared" = yes; then
      PBC_LDFLAGS="$PBC_LDFLAGS -no-undefined 0 -Wl,--export-all-symbols"
      LIBPBC_LDFLAGS="$LIBPBC_LDFLAGS -Wl,--output-def,.libs/libpbc.dll.def"
      LIBPBC_DLL=1
    fi
    ;;
esac
case $host in
  *-*-mingw*)
    gcc_cflags_optlist="$gcc_cflags_optlist nocygwin"
    gcc_cflags_nocygwin="-mno-cygwin"
    ;;
esac

AC_SUBST(PBC_LDFLAGS)
AC_SUBST(LIBPBC_LDFLAGS)
############################

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
lib_err_msg="add its path to LDFLAGS\nsee ./configure --help"
AC_CHECK_LIB( [m], [pow], [],[
	      echo "************************"
	      echo "m library not found"
	      echo -e $lib_err_msg
	      echo "************************"
	      exit -1
	      ])
AC_CHECK_LIB( [gmp], [__gmpz_init], [],[
	      echo "************************"
	      echo "gmp library not found"
	      echo -e $lib_err_msg
	      echo "************************"
	      exit -1
	      ])
dnl Reset libs because most programs do not need to link against all of these libs.
LIBS=

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

dnl setup CFLAGS
with_enable_optimized="no"
AC_ARG_ENABLE( optimized, 
	       [AC_HELP_STRING([--enable-optimized],
			[Enable optimized build])],
         [with_enable_optimized="$withval"],
	 [with_enable_optimized="no"])

CFLAGS="$CFLAGS -Wall -W -Wfloat-equal -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wredundant-decls \
-Wendif-labels -Wshadow -pipe -ffast-math -U__STRICT_ANSI__ -std=gnu99"
if test "x$with_enable_optimized" != "no"; then
   	CFLAGS="$CFLAGS -g -O2"
else
	CFLAGS="$CFLAGS -fomit-frame-pointer -O3"
fi

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gettimeofday memmove memset pow sqrt strchr strdup])

AC_CONFIG_FILES([Makefile example/Makefile gen/Makefile])
AC_OUTPUT

echo -ne "\n"
echo "global build variables"
echo "-----------------------------------------"
echo `date`
echo "host info:        $host"
echo "optimized build:  $with_enable_optimized"
echo "LDFLAGS:          $LDFLAGS"
echo "CPPFLAGS:         $CPPFLAGS"
echo "CFLAGS:           $CFLAGS"
echo "-----------------------------------------"
echo -ne "\n"

syscmd(bison -d -b pbc/parser pbc/parser.y)
syscmd(flex -o pbc/lex.yy.c --header-file=pbc/lex.yy.h pbc/parser.lex)
