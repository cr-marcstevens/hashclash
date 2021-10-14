#!/bin/bash

# check for special configuration for macosx

if [[ "$OSTYPE" == "darwin"* ]]; then
	SDK_PATH=$(xcrun --sdk macosx --show-sdk-path)
	# github actions macos-latest 1.57.0 fails, so try 1.70.0 instead
	: ${BOOST_VERSION:=1.72.0}
	echo "[*] Detected macosx: SDK_PATH=${SDK_PATH} BOOST_VERSION=${BOOST_VERSION}"
fi
if [ ! -z $(command -v glibtoolize) ]; then
	LIBTOOLIZE=glibtoolize
	echo "[*] Detected glibtoolize instead of libtoolize"
else
	LIBTOOLIZE=libtoolize
fi

# default values for boost library version and install location
# (if not already defined in environment)

: ${BOOST_VERSION:=1.57.0}
: ${BOOST_INSTALL_PREFIX:=$(pwd)/boost-$BOOST_VERSION}
: ${INCLUDE_DIRS:="/usr/include /usr/local/include $SDK_PATH/usr/include"}


# helper functions

function check_for_tool
{
	echo -n "[*] checking for system tool: $1: "
	if [ -z $(command -v $2) ]; then
		echo "*** error missing ***"
		exit 1
	fi
	echo "found"
}

function check_for_library
{
	echo -n "[*] checking for system library: $1: "
	have_lib=no
	for d in ${INCLUDE_DIRS}; do
		if [ -f $d/$2 ]; then
			have_lib=yes
			break
		fi
	done
	if [ "${have_lib}" = "no" ]; then
		echo "*** error missing ***"
		exit 1
	fi
	echo "found"
}


## Check prerequisites

check_for_tool autoconf autoreconf
check_for_tool automake automake
check_for_tool libtool $LIBTOOLIZE

check_for_library zlib1g-dev zlib.h
check_for_library libbz2-dev bzlib.h


# Check for local boost libraries and compile it if necessary

echo -n "[*] Checking for local boost (version ${BOOST_VERSION}): "
if [ ! -f ${BOOST_INSTALL_PREFIX}/include/boost/version.hpp ]; then
	echo "need to compile"
	echo "[*] Run: ./install_boost.sh"
	export BOOST_VERSION
	export BOOST_INSTALL_PREFIX
	(./install_boost.sh &> boost.log) || (echo "Boost build failed! See boost.log"; exit 1)
else
	echo "found"
fi
if [ ! -f ${BOOST_INSTALL_PREFIX}/include/boost/version.hpp ]; then
	if [ -f boost.log ]; then cat boost.log; fi
	echo "[*] *** error *** cannot find local boost (version ${BOOST_VERSION})"
	exit 1
fi


# Compile HashClash

echo "[*] Run: autoreconf --install"
autoreconf --install || exit 1
echo "[*] Run: ./configure --with-boost=${BOOST_INSTALL_PREFIX}"
./configure --with-boost=${BOOST_INSTALL_PREFIX} || exit 1
echo "[*] Run: make clean"
make clean
echo "[*] Run: make -j 4"
make -j 4 || exit 1
echo "[*] Finished!"
