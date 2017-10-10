AC_PREREQ([2.68])
AC_INIT([HashClash], [1.0.0], [marc@marc-stevens.nl],[hashclash],[https://github.com/cr-marcstevens/hashclash])
AC_CONFIG_SRCDIR([README.md])
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H([lib/config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])


DEFAULT_CXXFLAGS="-march=native -O3 -Wall"

CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

LT_INIT


# The version of the libtool library is of the form current:revision:age
#
# See http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
#
# When doing a release, they should be updated like this:
# 1. If no interfaces changed, only implementations: just increment
#    revision.
# 2. If interfaces were added, none removed: increment current, set
#    revision to zero and increment age.
# 3. If interfaces were removed (breaks backward compatibility): increment
#    current, and set both revision and age to zero.
HASHCLASH_LT_CURRENT=1
HASHCLASH_LT_REVISION=0
HASHCLASH_LT_AGE=0

AC_SUBST(HASHCLASH_LT_CURRENT)
AC_SUBST(HASHCLASH_LT_REVISION)
AC_SUBST(HASHCLASH_LT_AGE)

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S
AC_PROG_MKDIR_P
AC_PROG_CC
AC_PROG_CPP

AC_LANG([C++])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AX_BOOST_BASE([1.54],,[AC_MSG_ERROR([HashClash needs Boost])])

AX_PTHREAD()
LIBS="$PTHREAD_LIBS $LIBS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"

AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

AC_CONFIG_FILES([
 Makefile
])

AC_OUTPUT
