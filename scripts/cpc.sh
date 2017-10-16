#!/bin/bash

export BINDIR=$(dirname $0)/../bin
export BIRTHDAYSEARCH=$BINDIR/md5_birthdaysearch
export HELPER=$BINDIR/md5_diffpathhelper
export FORWARD=$BINDIR/md5_diffpathforward
export BACKWARD=$BINDIR/md5_diffpathbackward
export CONNECT=$BINDIR/md5_diffpathconnect
export CPUS=`cat /proc/cpuinfo | grep "^processor" | wc -l`
if [ ! -z $MAXCPUS ]; then
	if [ $CPUS -gt $MAXCPUS ]; then
		export CPUS=$MAXCPUS
	fi
fi
export TTT=12

rm -r data 2>/dev/null
mkdir data || exit 1

file1=$1
file2=$2

if [ -f "$file1" -a -f "$file2" ]; then
	echo "Chosen-prefix file 1: $file1"
	echo "Chosen-prefix file 2: $file2"
else
	echo "Usage $0 <file1> <file2> [<redonearcollstep>] [<nobirthday>]"
	exit 1
fi

if [ "$3" = "" ] ; then
	$BIRTHDAYSEARCH --inputfile1 "$file1" --inputfile2 "$file2" --hybridbits 0 --pathtyperange 2 --maxblocks 9 --maxmemory 100 --threads $CPUS
else
	if [ "$4" != "" ]; then
		cp $file1 file1.bin
		cp $file2 file2.bin
	fi
fi

function doforward {
	$FORWARD -w $1 -f $1/lowerpath.bin.gz --normalt01 -t 1 --trange $(($TTT-3)) --threads $CPUS || return 1
	$FORWARD -w $1 -a 500000 -t $(($TTT-1)) --trange 0 --threads $CPUS || return 1
}

function dobackward {
	$BACKWARD -w $1 -f $1/upperpath.bin.gz -t 34 --trange 4 -a 65536 -q 128 --threads $CPUS || return 1
	$BACKWARD -w $1 -t 29 -a 100000 --trange 8 --threads $CPUS || return 1
	$BACKWARD -w $1 -t 20 -a 16384 --threads $CPUS || return 1
	$BACKWARD -w $1 -t 19 -a 500000 --trange $((18-$TTT-3)) --threads $CPUS || return 1
}

function testcoll {
	for f in $1/coll1_*; do
		if [ -e `echo $f | sed s/coll1/coll2/` ]; then return 0; fi
	done
	return 1
}

function doconnect {
	let c=0

	mkdir $1/connect
	$CONNECT -w $1/connect -t $TTT --inputfilelow $1/paths$(($TTT-1))_0of1.bin.gz --inputfilehigh $1/paths$(($TTT+4))_0of1.bin.gz --threads $CPUS &
	CPID="$!"

	let contime=0
	while true; do
		sleep 5
		ps -p $CPID &>/dev/null || break
		let contime=contime+$((5*$CPUS))
		if [ $contime -gt 10000 ]; then
			kill $CPID
			break;
		fi
	done
}

function docollfind {
	$HELPER -w $1 --findcoll $2/bestpaths.bin.gz --threads $CPUS |& tee $2/collfind.log &
	CPID=$!
	while true; do
		if testcoll $1 ; then break; fi
		sleep 1
	done
	sleep 3
	kill $CPID
}

cp file1.bin file1_0.bin
cp file2.bin file2_0.bin
let k=$3
if [ "$k" = "" ]; then
	let k=0
fi

cat <<EOF >md5diffpathbackward.cfg.template
autobalance = 1000000
estimate = 4
fillfraction = 1
maxsdrs = 1
condtend = 35
EOF

cat <<EOF >md5diffpathforward.cfg.template
autobalance = 2000000
estimate = 4
fillfraction = 1
maxsdrs = 1
minQ456tunnel = 18
minQ91011tunnel = 18
EOF

cat <<EOF >md5diffpathconnect.cfg.template
Qcondstart = 21
EOF

while true; do
	if [ -f ${file1}.coll -a -f ${file2}.coll ]; then
		if [ `cat ${file1}.coll | md5sum | cut -d' ' -f1` = `cat ${file2}.coll | md5sum | cut -d' ' -f1` ]; then
			echo "Collision generated: ${file1}.coll ${file2}.coll"
			md5sum ${file1}.coll ${file2}.coll
			exit
		fi
	fi

	rm -r workdir$k 2>/dev/null
	mkdir workdir$k
	( $HELPER -w workdir$k --startnearcollision file1_$k.bin file2_$k.bin --pathtyperange 2 |& tee workdir$k/start.log ) || exit
	cp *.cfg workdir$k

	doforward workdir$k |& tee workdir$k/forward.log
	dobackward workdir$k |& tee workdir$k/backward.log

	doconnect workdir$k |& tee workdir$k/connect.log

	docollfind workdir$k workdir$k/connect

	for f in workdir$k/coll1_*; do
		if [ -e `echo $f | sed s/coll1/coll2/` ]; then
			cat file1_$k.bin $f > file1_$(($k+1)).bin
			cat file2_$k.bin `echo $f | sed s/coll1/coll2/` > file2_$(($k+1)).bin
			cp file1_$(($k+1)).bin ${file1}.coll
			cp file2_$(($k+1)).bin ${file2}.coll
			break;
		fi
	done
	let k=k+1
done
