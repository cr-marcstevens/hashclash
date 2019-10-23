#!/bin/sh

: ${BOOST_VERSION:=1.57.0}
: ${BOOST_INSTALL_PREFIX:=$(pwd)/boost-$BOOST_VERSION}

## Check prerequisites

echo -n "[*] Checking for system packages: autoconf, automake & libtool: "
if [ -z $(command -v autoreconf) ] || [ -z $(command -v automake) ] || [ -z $(command -v libtoolize) ]; then
	echo "*** error missing ***"
	exit 1
fi
echo "found"

echo -n "[*] Checking for libraries: zlib and bz2: "
have_zlib=no
have_bz2=no
for d in /usr/include /usr/local/include; do
	if [ -f $d/zlib.h ]; then have_zlib=yes; fi
	if [ -f $d/bzlib.h ]; then have_bz2=yes; fi
done
if [ "${have_zlib}" = "no" ] || [ "${have_bz2}" = "no" ]; then
	echo "*** error missing ***"
	exit 1
fi
echo "found"

echo -n "[*] Checking for local boost (version ${BOOST_VERSION}): "
if [ ! -f ${BOOST_INSTALL_PREFIX}/include/boost/version.hpp ]; then
	echo "need to compile"
	export BOOST_VERSION
	export BOOST_INSTALL_PREFIX
	./install_boost.sh || exit 1
else
	echo "found"
fi
if [ ! -f ${BOOST_INSTALL_PREFIX}/include/boost/version.hpp ]; then
	echo "[*] Cannot find local boost (version ${BOOST_VERSION})"
fi

echo "[*] Run: autoreconf --install"
autoreconf --install || exit 1
echo "[*] Run: ./configure --with-boost=${BOOST_INSTALL_PREFIX}"
./configure --with-boost=${BOOST_INSTALL_PREFIX} || exit 1
echo "[*] Run: make clean"
make clean
echo "[*] Run: make -j 4"
make -j 4
