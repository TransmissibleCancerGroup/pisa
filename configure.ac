# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(pisa, 0.1, kg8@sanger)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_MAKE_SET

# Checks for header files.
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])
AC_CHECK_HEADER([lzma.h])

# Checks for libraries.
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([strdup])

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")
AC_SUBST(LDFLAGS, "$LDFLAGS -pthread -std=c++11")

# Check for hash table headers
AC_CHECK_HEADERS([ \
	google/sparse_hash_set google/sparse_hash_map \
	unordered_map tr1/unordered_map ext/hash_map \
	unordered_set tr1/unordered_set ext/hash_set \
])

AC_CONFIG_FILES([
	Makefile
	src/Makefile])

AC_CONFIG_SUBDIRS([deps/SeqLib])

AC_OUTPUT
