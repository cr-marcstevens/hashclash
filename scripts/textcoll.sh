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

#ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.,_-~=+:;|?@#^&*(){}[]<>"
ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.,_-~=+:;|?@#^&*"
#ALPHABET="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

# First block: byte21 has a +4 difference: you must ensure that this difference is possible in the alphabet, e.g., have both A and E(=A+4)
FIRSTBLOCKBYTES='--byte0 T --byte1 E --byte2 X --byte3 T --byte4 C --byte5 O --byte6 L --byte7 L --byte20 hH --byte21 aAeE --byte22 cC --byte23 kK'

# Second block: keep the alphabet of bytes 0-7 large: otherwise there could be no solutions
SECONDBLOCKBYTES='--byte8 B --byte9 y --byte10 M --byte11 a --byte12 r --byte13 c --byte14 S --byte15 t --byte16 e --byte17 v --byte18 e --byte19 n --byte20 s'

$HELPER $MDIFF --pathfromtext --inputfile1 $DIFFPATH --outputfile1 textcoll_path.bin.gz || exit 1

$TEXTCOLL $MDIFF $FIRSTBLOCKBYTES --prepare --pathfile textcoll_path.bin.gz --alphabet $ALPHABET --prefixfile ${prefixfile} || exit 1

$TEXTCOLL $MDIFF $FIRSTBLOCKBYTES --firstblock --pathfile textcoll_path.bin.gz --alphabet $ALPHABET --prefixfile ${prefixfile} || exit 1

collfile=`ls textcoll1_*.txt | head -n1`

cat ${prefixfile} ${collfile} > partial_solution.txt

$TEXTCOLL $MDIFF $SECONDBLOCKBYTES --secondblock --alphabet $ALPHABET --prefixfile partial_solution.txt || exit 1
