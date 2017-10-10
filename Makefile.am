EXTRA_DIST       = README.md LICENSE.TXT
ACLOCAL_AMFLAGS  = -I m4

noinst_LTLIBRARIES=lib/libhashclash.la

lib_libhashclash_la_SOURCES=\
	lib/hashclash/booleanfunction.cpp lib/hashclash/booleanfunction.hpp \
	lib/hashclash/conditions.cpp lib/hashclash/conditions.hpp \
	lib/hashclash/cpuperformance.hpp \
	lib/hashclash/differentialpath.cpp lib/hashclash/differentialpath.hpp \
	lib/hashclash/md5detail.cpp lib/hashclash/md5detail.hpp \
	lib/hashclash/progress_display.hpp \
	lib/hashclash/rng.cpp lib/hashclash/rng.hpp \
	lib/hashclash/saveload_bz2.hpp lib/hashclash/saveload_gz.hpp lib/hashclash/saveload.hpp \
	lib/hashclash/sdr.cpp lib/hashclash/sdr.hpp \
	lib/hashclash/sha1detail.cpp lib/hashclash/sha1detail.hpp \
	lib/hashclash/sha1differentialpath.cpp lib/hashclash/sha1differentialpath.hpp \
	lib/hashclash/sha1messagespace.cpp lib/hashclash/sha1messagespace.hpp \
	lib/hashclash/timer.cpp lib/hashclash/timer.hpp \
	lib/hashclash/types.hpp

lib_libhashclash_la_LIBADD=-lboost_thread -lboost_serialization -lboost_program_options -lboost_filesystem -lboost_iostreams -lboost_system -lz -lbz2
lib_libhashclash_la_LDFLAGS=-no-undefined -version-info @HASHCLASH_LT_CURRENT@:@HASHCLASH_LT_REVISION@:@HASHCLASH_LT_AGE@

bin_PROGRAMS=\
	bin/md5_diffpathforward \
	bin/md5_diffpathbackward \
	bin/md5_diffpathconnect \
	bin/md5_diffpathhelper

AM_CPPFLAGS=-I$(top_srcdir)/lib
LDADD=lib/libhashclash.la

bin_md5_diffpathforward_SOURCES=\
	src/md5forward/dostep.cpp \
	src/md5forward/forward.cpp \
	src/md5forward/main.cpp \
	src/md5forward/main.hpp

bin_md5_diffpathbackward_SOURCES=\
	src/md5backward/backward.cpp \
	src/md5backward/dostep.cpp \
	src/md5backward/main.cpp \
	src/md5backward/main.hpp

bin_md5_diffpathconnect_SOURCES=\
	src/md5connect/connect.cpp \
	src/md5connect/dostep.cpp \
	src/md5connect/main.cpp \
	src/md5connect/main.hpp

bin_md5_diffpathhelper_SOURCES=\
	src/md5helper/collisionfinding.cpp \
	src/md5helper/convert.cpp \
	src/md5helper/main.cpp \
	src/md5helper/main.hpp \
	src/md5helper/startnearcollision.cpp
