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

if [ ! -f textcoll_path.bin.gz ]; then
	$HELPER $MDIFF --pathfromtext --inputfile1 $DIFFPATH --outputfile1 textcoll_path.bin.gz || exit 1
fi

if [ ! -f Q7Q24.bin.gz ]; then
	( $TEXTCOLL $MDIFF $FIRSTBLOCKBYTES --prepare --pathfile textcoll_path.bin.gz --alphabet $ALPHABET --prefixfile ${prefixfile} | tee prepare.log ) || exit 1
	echo "Solutions stored in Q7Q24.bin.gz."
	echo "It's possible to try to control more bytes by editing this script and deleting Q7Q24.bin.gz."
fi

collfile=`ls textcoll1_block1_[0-9]*.txt | head -n1`

if [ "$collfile" = "" ]; then
	echo "Starting search for first near-collision block in 10 seconds..."
	sleep 10
	
	( $TEXTCOLL $MDIFF $FIRSTBLOCKBYTES --firstblock --pathfile textcoll_path.bin.gz --alphabet $ALPHABET --prefixfile ${prefixfile} | tee firstblock.log ) || exit 1
	collfile=`ls textcoll1_*.txt | head -n1`
fi

echo "Found first near-collision block in file: $collfile"

cat ${prefixfile} ${collfile} > partial_solution.txt

echo "Starting search for second near-collision block in 10 seconds..."
sleep 10

$TEXTCOLL $MDIFF $SECONDBLOCKBYTES --secondblock --alphabet $ALPHABET --prefixfile partial_solution.txt || exit 1

# second block file
coll2file=`ls textcoll1_block2_[0-9]*.txt | head -n1`

# other first block file
collfile2=`echo $collfile | sed "s/textcoll1/textcoll2/"`

cat ${prefixfile} ${collfile} ${coll2file} > final_collision1.txt
cat ${prefixfile} ${collfile2} ${coll2file} > final_collision2.txt

md5sum final_collision*
sha1sum final_collision*

echo "========= final_collision1.txt =========="
cat final_collision1.txt
echo "========= final_collision2.txt =========="
cat final_collision2.txt
