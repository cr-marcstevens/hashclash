SUFFIXES=.cu

EXTRA_DIST       = README.md LICENSE.TXT
ACLOCAL_AMFLAGS  = -I m4

noinst_LTLIBRARIES=lib/libhashclash.la

lib_libhashclash_la_SOURCES=\
	lib/hashclash/bestof.hpp \
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

lib_libhashclash_la_LDFLAGS=-no-undefined

bin_PROGRAMS=\
	bin/md5_fastcoll \
	bin/md5_diffpathforward \
	bin/md5_diffpathbackward \
	bin/md5_diffpathconnect \
	bin/md5_diffpathhelper \
	bin/md5_birthdaysearch \
	bin/sha1_diffpathforward \
	bin/sha1_diffpathbackward \
	bin/sha1_diffpathconnect \
	bin/sha1_diffpathhelper \
	bin/sha1_attackgenerator \
	bin/sha1_nearcollisionattack

AM_CPPFLAGS=-I$(top_srcdir)/lib
LDADD=lib/libhashclash.la

bin_md5_fastcoll_SOURCES=\
	src/md5fastcoll/main.cpp \
	src/md5fastcoll/md5.cpp \
	src/md5fastcoll/block0.cpp \
	src/md5fastcoll/block1.cpp \
	src/md5fastcoll/block1wang.cpp \
	src/md5fastcoll/block1stevens00.cpp \
	src/md5fastcoll/block1stevens01.cpp \
	src/md5fastcoll/block1stevens10.cpp \
	src/md5fastcoll/block1stevens11.cpp

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

bin_md5_birthdaysearch_SOURCES=\
	src/md5birthdaysearch/birthday.cpp \
	src/md5birthdaysearch/birthday_types.hpp \
	src/md5birthdaysearch/config.h \
	src/md5birthdaysearch/distribution.hpp \
	src/md5birthdaysearch/dostep.cpp \
	src/md5birthdaysearch/main.cpp \
	src/md5birthdaysearch/main.hpp \
	src/md5birthdaysearch/storage.cpp \
	src/md5birthdaysearch/storage.hpp \
	src/md5birthdaysearch/simd_avx256.cpp

if HAVE_CUDA
bin_md5_birthdaysearch_SOURCES+=\
	src/md5birthdaysearch/cuda_md5.cu
bin_md5_birthdaysearch_LDADD=lib/libhashclash.la $(CUDA_LIBS)
else
bin_md5_birthdaysearch_LDADD=lib/libhashclash.la
endif

bin_sha1_diffpathforward_SOURCES=\
	src/sha1forward/dostep.cpp \
	src/sha1forward/forward.cpp \
	src/sha1forward/main.cpp \
	src/sha1forward/main.hpp

bin_sha1_diffpathbackward_SOURCES=\
	src/sha1backward/backward.cpp \
	src/sha1backward/dostep.cpp \
	src/sha1backward/main.cpp \
	src/sha1backward/main.hpp

bin_sha1_diffpathconnect_SOURCES=\
	src/sha1connect/connect.cpp \
	src/sha1connect/dostep.cpp \
	src/sha1connect/main.cpp \
	src/sha1connect/main.hpp

bin_sha1_diffpathhelper_SOURCES=\
	src/sha1helper/convert.cpp \
	src/sha1helper/filter.cpp \
	src/sha1helper/main.cpp \
	src/sha1helper/main.hpp

bin_sha1_attackgenerator_SOURCES=\
	src/sha1attackgenerator/checkokpaths.cpp \
	src/sha1attackgenerator/collisionfinding.cpp \
	src/sha1attackgenerator/main.cpp \
	src/sha1attackgenerator/tunnel_analysis.cpp \
	src/sha1attackgenerator/filters.cpp \
	src/sha1attackgenerator/main.hpp \
	src/sha1attackgenerator/program_generator.cpp

bin_sha1_nearcollisionattack_SOURCES=\
	src/sha1attackgenerator/collfind.cpp

NVCCFLAGS=-ccbin $(CXX) -O2 -v --maxrregcount=64 --ptxas-options=-v -Xcompiler -mcmodel=medium

CUDA_SMS=50 52 60 61 70 75
GENCODE_FLAGS=
$(foreach sm,$(CUDA_SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
HIGHEST_SM:=$(lastword $(sort $(CUDA_SMS)))
GENCODE_FLAGS+= -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)

.cu.o:
	$(NVCC) $(GENCODE_FLAGS) $(NVCCFLAGS) $(CUDA_CFLAGS) $(addprefix -Xcompiler ,$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)) -o $@ -c $<
