#!/bin/bash

: ${BOOST_VERSION:=1.57.0}
: ${BOOST_BUILD_OPTIONS:="-j4 --with-filesystem --with-iostreams --with-program_options --with-serialization --with-system --with-thread"}
: ${BOOST_BUILD_CXXFLAGS:="-O2"}
: ${BOOST_INSTALL_PREFIX:=$(pwd)/boost-$BOOST_VERSION}

if [ ! -z $BOOST_TOOLSET ]; then
	BOOST_BUILD_OPTIONS="$BUILD_OPTIONS toolset=$BOOST_TOOLSET"
fi
echo "Building boost version : $BOOST_VERSION"
echo "Installing in directory: $BOOST_INSTALL_PREFIX"
echo "Build options          : $BOOST_BUILD_OPTIONS"
echo "CXXFLAGS               : $BOOST_BUILD_CXXFLAGS"
if [ -f config.status ]; then
	CXXSTD1=`grep -o "\-std=c++[0-9a-zA-Z]*" config.status 2>/dev/null | sort -u`
	CXXSTD2=`echo $BOOST_BUILD_CXXFLAGS | grep -o "\-std=c++[0-9a-zA-Z]*"`
	if [ "$CXXSTD1" != "$CXXSTD2" ]; then
		echo "*** WARNING: detected use of $CXXSTD1 in config.status, add this flag to CXXFLAGS:"
		echo "***    BOOST_BUILD_CXXFLAGS=\"$BOOST_BUILD_CXXFLAGS $CXXSTD1\" $0"
	fi
fi
echo "Starting in 3 seconds..."
sleep 3

BOOST_DIRVERSION=`echo $BOOST_VERSION | tr . _`
BOOST_DIR=boost_$BOOST_DIRVERSION
BOOST_FILE=$BOOST_DIR.tar.gz

function test_tar_gz_file
{
	if [ ! -f $TGZ_FILE ]; then
		return 1
	fi
	TGZ_FILE=$1
	gzip -t $TGZ_FILE \
		|| ( echo "$TGZ_FILE corrupted" ; rm $TGZ_FILE ; return 2 )
	return 0
}

function retrieve_tar_gz_file
{
	TGZ_FILE=$1
	TGZ_URL=$2
	
	if test_tar_gz_file $TGZ_FILE ; then 
		return 0
	fi
	wget $TGZ_URL -O $TGZ_FILE.tmp \
		|| ( echo "wget $TGZ_FILE failed" ; rm $TGZ_FILE.tmp ; return 3 )
	mv $TGZ_FILE.tmp $TGZ_FILE
	if test_tar_gz_file $TGZ_FILE ; then 
		return 0
	fi
	return 1	
}

# retrieve file if necessary
retrieve_tar_gz_file $BOOST_FILE https://downloads.sourceforge.net/project/boost/boost/$BOOST_VERSION/$BOOST_FILE
retrieve_tar_gz_file $BOOST_FILE https://downloads.sourceforge.net/project/boost/boost/$BOOST_VERSION/$BOOST_FILE

# test if we have file
test_tar_gz_file $BOOST_FILE \
	|| ( echo "Cannot find or download $BOOST_FILE" ; exit 1 )

# remove existing boost directory if present
if [ -d $BOOST_DIR ]; then
	echo -n "Removing old directory $BOOST_DIR..."
	rm -rf "$BOOST_DIR"
	echo "done."
fi

# extract boost tar file
tar -xzf $BOOST_FILE \
	|| ( echo "Cannot extract $BOOST_FILE" ; exit 1 )

# go into boost directory
if [ ! -d $BOOST_DIR ]; then
	echo "Expected directory $BOOST_DIR after extracting $BOOST_FILE"
	exit 1
fi
cd $BOOST_DIR


# bootstrap bjam build system for boost
./bootstrap.sh \
	cxxflags="$BOOST_BUILD_CXXFLAGS" \
	linkflags="$BOOST_BUILD_CXXFLAGS" \
	--prefix=$BOOST_INSTALL_PREFIX

if [ ! -e ./b2 ]; then
	cat bootstrap.log
fi

# compile boost
./b2 -q \
	cxxflags="$BOOST_BUILD_CXXFLAGS" \
	linkflags="$BOOST_BUILD_CXXFLAGS" \
	--prefix=$BOOST_INSTALL_PREFIX \
	$BOOST_BUILD_OPTIONS \
	--build-type=minimal link=static \
	install \
	|| ( echo "Building boost failed!" ; exit 1 )

# exit boost directory
cd ..


# cleanup boost directory
echo -n "Removing old directory $BOOST_DIR..."
rm -rf $BOOST_DIR
echo "done."
