#!/bin/bash

export BINDIR=$(dirname $0)/../bin
export HELPER=$BINDIR/md5_diffpathhelper
export FORWARD=$BINDIR/md5_diffpathforward
export BACKWARD=$BINDIR/md5_diffpathbackward
export CONNECT=$BINDIR/md5_diffpathconnect

N=1
tconnect=12
data=200000
prefixfile=prefix.bin

if [ ! -z $1 ]; then
	prefixfile="$1"
fi
if [ ! -z $prefixfile ]; then
	if [ ! -f $prefixfile ]; then
		echo "Not using prefixfile: doesn't exist"
		prefixfile=
	fi
fi

# Note that positive message bit differences are given with bit indices 1,...,32
# and with -1,...,-32 for a negative bit difference
case $N in 

1)
diff="--diffm2 9"
tmid=32
;;

2)
diff="--diffm13 28 --diffm0 32 --diffm6 32"
tmid=32
;;

3)
diff="--diffm6 9 --diffm9 32 --diffm15 32"
tmid=32
tconnect=10
data=500000
;;

esac

opts="-a $data -e 4 --fillfraction 1 -q 8"

# automatically kill connect 10 seconds after first found differential path
# but don't if abortfile exists
# (on some machines the first NC block process can be so fast, 
#  the first auto_kill_connect will actually inadvertently kill the second NC process)
function auto_kill_connect
{
	abortfile="$1"
	while [ ! -f data/bestpath.bin.gz ]; do
		sleep 1
	done
	sleep 10
	if [ "$abortfile" != "" ]; then
		if [ ! -f "$abortfile" ]; then
			killall md5_diffpathconnect
		fi
	fi
}

rm -f *.log
rm -rf data
mkdir -p data




###################### FIRST NEAR-COLLISION BLOCK PAIR #####################

#create lower paths
if [ -z $prefixfile ]; then
cat <<EOF > path_cv.txt
Q-3:    |........ ........ ........ ........| ok p=1
Q-2:    |........ ........ ........ ........| ok p=1
Q-1:    |........ ........ ........ ........| ok p=1
Q0:     |........ ........ ........ ........| ok p=1
EOF
$HELPER $diff --pathfromtext --inputfile1 path_cv.txt --outputfile1 path_cv.bin
echo "Continuing in 3 seconds..."
sleep 3
$FORWARD $diff $opts -f path_cv.bin --normalt01 -t 0 --trange $(($tconnect-1)) 2>&1 | tee -a lower.log
else
$HELPER $diff --startpartialpathfromfile --inputfile1 $prefixfile --outputfile1 $prefixfile.base --outputfile2 path_prefix.bin 2>&1 | tee lower.log
echo "Continuing in 3 seconds..."
sleep 3
startt=`cat lower.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
$FORWARD $diff $opts -f path_prefix.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a lower.log
fi

#create upper paths
mkdir -p upper$N
if [ ! -f upper$N/done ]; then
	$BACKWARD $diff -w upper$N -n -t $tmid -b $(($tmid+1))
	$FORWARD $diff $opts -w upper$N -t $((tmid+1)) -b $tmid --trange $((61-$tmid)) 2>&1 | tee upper.log
	$FORWARD $diff -w upper$N -a 16 -t 63 -b $tmid 2>&1 | tee -a upper.log

	$BACKWARD $diff $opts -w upper$N -f upper$N/paths63_0of1.bin.gz -t $((tmid-1)) --trange $(($tmid-$tconnect-5)) 2>&1 | tee -a upper.log
	touch upper$N/done
fi

#connect
rm -f data/bestpath.bin.gz
auto_kill_connect data/connect1done &
$CONNECT $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh upper$N/paths$((tconnect+4))_0of1.bin.gz 2>&1 | tee connect.log
touch data/bestpath.bin.gz data/connect1done
sleep 2

#find example collision
$HELPER $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee collfind.log

################## SECOND NEAR-COLLISION BLOCK PAIR #########################

# swap msg1 and msg2 to use opposing signs as first block
cat $prefixfile.base ./data/coll1_* > ./data/file2.bin
cat $prefixfile.base ./data/coll2_* > ./data/file1.bin
mv ./data/coll1_* ./data/nc1_coll1.bin
mv ./data/coll2_* ./data/nc1_coll2.bin
mv data/bestpath.bin.gz data/nc1_bestpath.bin.gz

$HELPER $diff --startpartialpathfromfile --inputfile1 ./data/file1.bin --inputfile2 ./data/file2.bin --outputfile1 /dev/null --outputfile2 path_2nc.bin 2>&1 | tee lower2.log
startt=`cat lower2.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
$FORWARD $diff $opts -f path_2nc.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a lower2.log

#connect
rm -f data/bestpath.bin.gz
auto_kill_connect data/connect2done &
$CONNECT $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh upper$N/paths$((tconnect+4))_0of1.bin.gz 2>&1 | tee connect2.log
touch data/bestpath.bin.gz data/connect2done

#find example collision
$HELPER $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee collfind2.log

cat data/file1.bin data/coll1_* > collision2.bin
cat data/file2.bin data/coll2_* > collision1.bin
md5sum collision1.bin collision2.bin
sha1sum collision1.bin collision2.bin
ls -als collision1.bin collision2.bin
