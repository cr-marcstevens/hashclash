#!/bin/bash

BOOST_VERSION=1.57.0

BOOST_BUILD_OPTIONS="-j4 --without-graph --without-graph_parallel --without-python --without-regex --without-wave "
BOOST_BUILD_CXXFLAGS="-std=c++11"

BOOST_DIRVERSION=`echo $BOOST_VERSION | tr . _`
BOOST_DIR=boost_$BOOST_DIRVERSION
BOOST_FILE=$BOOST_DIR.tar.gz

PWD=`pwd`

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

# test if we have file
test_tar_gz_file $BOOST_FILE \
	|| ( echo "Cannot find or download $BOOST_FILE" ; exit 1 )

# remove existing boost directory if present
if [ -d $BOOST_DIR ]; then
	echo -n "Removing old directory $BOOST_DIR..."
	rm -rf $BOOST_DIR
	echo "done."
fi

# extract boost tar file
tar -xzvf $BOOST_FILE \
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
	--prefix=$PWD/boost

# compile boost
./b2 \
	cxxflags="$BOOST_BUILD_CXXFLAGS" \
	linkflags="$BOOST_BUILD_CXXFLAGS" \
	--prefix=$PWD/boost \
	$BOOST_BUILD_OPTIONS \
	--build-type=minimal \
	install

# exit boost directory
cd ..
