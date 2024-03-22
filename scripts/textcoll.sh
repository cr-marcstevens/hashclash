#!/bin/bash

export BINDIR=$(dirname $0)/../bin
export HELPER=$BINDIR/md5_diffpathhelper
export TEXTCOLL=$BINDIR/md5_textcoll

export DIFFPATH=$BINDIR/../src/md5textcoll/path2.txt
MDIFF="--diffm5 11"

prefixfile=$1

if [ -z $prefixfile ]; then
	prefixfile=dummy.prefix.bin
fi
if [ ! -f $prefixfile ]; then
	touch $prefixfile
fi

ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-~=+:;|?@#^&*(){}[]<>"
#ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-~=+:;|?@#^&*"
#ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

$HELPER $MDIFF --pathfromtext --inputfile1 $DIFFPATH --outputfile1 textcoll_path.bin.gz || exit 1

$TEXTCOLL $MDIFF --prepare --pathfile textcoll_path.bin.gz --alphabet "$ALPHABET" --prefixfile ${prefixfile} || exit 1

$TEXTCOLL $MDIFF --firstblock --pathfile textcoll_path.bin.gz --alphabet "$ALPHABET" --prefixfile ${prefixfile} || exit 1

collfile=`ls textcoll1_*.txt | head -n`

cat ${prefixfile} ${collfile} > partial_solution.txt

$TEXTCOLL $MDIFF --secondblock --alphabet "$ALPHABET" --prefixfile partial_solution.txt || exit 1

