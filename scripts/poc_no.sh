#!/bin/bash

export BINDIR=$(dirname $0)/../bin
export HELPER=$BINDIR/md5_diffpathhelper
export FORWARD=$BINDIR/md5_diffpathforward
export BACKWARD=$BINDIR/md5_diffpathbackward
export CONNECT=$BINDIR/md5_diffpathconnect

N=1
tconnect=12
datalower=160000
dataupper=640000
prefixfile=$1

if [ -z $prefixfile ]; then
	prefixfile=dummy.prefix.bin
fi
if [ ! -f $prefixfile ]; then
	touch $prefixfile
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
	datalower=500000
	dataupper=500000
	;;

esac

# automatically kill connect 10 seconds after first found differential path
# only if pidfile exists
function auto_kill_connect
{
	pidfile="$1"
	while [ ! -s data/bestpath.bin.gz ]; do
		sleep 1
	done
	sleep 10
	if [ "$pidfile" != "" ]; then
		if [ -s "$pidfile" ]; then
			echo "Found differential paths, killing md5diffpathconnect"
			kill $(<$pidfile)
		fi
	fi
}

#global loop that ensures that the script executes successfully
while true; do

	#clean the state
	rm -f *.log
	rm -rf data
	mkdir -p data
	mkdir -p logs
	rm md5diffpath*.cfg md5diffpath*.template

	# update configuration (datalower may be increased in the global loop)
	optsupper="-a ${dataupper} -e 4 --fillfraction 1 -q 8"
	optslower="-a ${datalower} -e 4 --fillfraction 1 -q 8"
	updir=upper_${N}_${dataupper}

###################### FIRST NEAR-COLLISION BLOCK PAIR #####################

	$HELPER $diff --startpartialpathfromfile --inputfile1 $prefixfile --outputfile1 data/$prefixfile.base --outputfile2 data/path_prefix.bin 2>&1 | tee logs/lower.log
	echo "Continuing in 3 seconds..."
	sleep 3
	startt=`cat logs/lower.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
	$FORWARD $diff $optslower -f data/path_prefix.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a logs/lower.log

	#create upper paths
	mkdir -p  $updir
	if [ ! -f $updir/done ]; then
		$BACKWARD $diff -w $updir -n -t $tmid -b $(($tmid+1))
		$FORWARD $diff $optsupper -w $updir -t $((tmid+1)) -b $tmid --trange $((61-$tmid)) 2>&1 | tee logs/upper.log
		$FORWARD $diff -w $updir -a 16 -t 63 -b $tmid 2>&1 | tee -a logs/upper.log

		$BACKWARD $diff $optsupper -w $updir -f $updir/paths63_0of1.bin.gz -t $((tmid-1)) --trange $(($tmid-$tconnect-5)) 2>&1 | tee -a logs/upper.log
		touch $updir/done
	fi

	#connect
	rm -f data/bestpath.bin.gz
	auto_kill_connect data/connect1pid &
	autokiller1=$!
	trap "kill -TERM $autokiller1; exit" TERM INT
	($CONNECT $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh $updir/paths$((tconnect+4))_0of1.bin.gz 2>&1 & echo $! >&3) 3>data/connect1pid | tee logs/connect.log
	rm -f data/connect1pid
	kill $autokiller1
	sleep 2

	#check that the first connection happened
	if [ ! -s data/bestpath.bin.gz ]; then let datalower=$((2*$datalower)); echo "First step failed. Restarting"; sleep 2; continue; fi


	#find example collision
	$HELPER $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee logs/collfind.log

################## SECOND NEAR-COLLISION BLOCK PAIR #########################

	# swap msg1 and msg2 to use opposing signs as first block
	cat data/$prefixfile.base ./data/coll1_* > ./data/file2.bin
	cat data/$prefixfile.base ./data/coll2_* > ./data/file1.bin
	mv ./data/coll1_* ./data/nc1_coll1.bin
	mv ./data/coll2_* ./data/nc1_coll2.bin
	mv data/bestpath.bin.gz data/nc1_bestpath.bin.gz

	$HELPER $diff --startpartialpathfromfile --inputfile1 ./data/file1.bin --inputfile2 ./data/file2.bin --outputfile1 /dev/null --outputfile2 data/path_2nc.bin 2>&1 | tee logs/lower2.log
	startt=`cat logs/lower2.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
	$FORWARD $diff $optslower -f data/path_2nc.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a logs/lower2.log

	#connect
	rm -f data/bestpath.bin.gz
	auto_kill_connect data/connect2pid &
	autokiller2=$!
	trap "kill -TERM $autokiller2; exit" TERM INT
	($CONNECT $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh $updir/paths$((tconnect+4))_0of1.bin.gz 2>&1 & echo $! >&3) 3>data/connect2pid | tee logs/connect2.log
	rm -f data/connect2pid
	kill $autokiller2
	sleep 2

	#check that the second connection happened
	if [ ! -s data/bestpath.bin.gz ]; then let datalower=$((2*$datalower)); echo "Second step failed. Restarting"; sleep 2; continue; fi

	#find example collision
	$HELPER $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee logs/collfind2.log

	cat data/file1.bin data/coll1_* > collision2.bin
	cat data/file2.bin data/coll2_* > collision1.bin
	md5sum collision1.bin collision2.bin
	sha1sum collision1.bin collision2.bin
	ls -als collision1.bin collision2.bin

	#check that the collision is a collision
	if [ `md5sum collision1.bin | cut -d' ' -f1` != `md5sum collision2.bin | cut -d' ' -f1` ]; then echo "Not a real collision. Restarting"; sleep 2; continue; fi

	#closing the global loop
	break

done
