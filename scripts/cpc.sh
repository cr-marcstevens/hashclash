#!/bin/bash

export BINDIR=../bin-1.57.0
export BIRTHDAYSEARCH=$BINDIR/md5_birthdaysearch
export HELPER=$BINDIR/md5_diffpathhelper
export FORWARD=$BINDIR/md5_diffpathforward
export BACKWARD=$BINDIR/md5_diffpathbackward
export CONNECT=$BINDIR/md5_diffpathconnect
export CPUS=`cat /proc/cpuinfo | grep "^processor" | wc -l`
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
	$BIRTHDAYSEARCH --inputfile1 "$file1" --inputfile2 "$file2" --hybridbits 0 --pathtyperange 2 --maxblocks 19 --maxmemory 100 --threads $CPUS
else
	if [ "$4" != "" ]; then
		cp $file1 file1.bin
		cp $file2 file2.bin
	fi
fi

function doforward {
	$FORWARD -w $1 -f $1/lowerpath.bin.gz --normalt01 -t 1 --trange $(($TTT-2)) || exit 1
}

function dobackward {
	$BACKWARD -w $1 -f $1/upperpath.bin.gz -t 34 --trange 4 -a 65536 -q 128 || exit 1
	$BACKWARD -w $1 -t 29 --trange 8 || exit 1
	$BACKWARD -w $1 -t 20 -a 16384 || exit 1
	$BACKWARD -w $1 -t 19 --trange $((18-$TTT-4)) || exit 1
	$BACKWARD -w $1 -t $(($TTT+4)) --maxweight 13 --minweight 10 || exit 1
}

function testcoll {
	for f in $1/coll1_*; do
		if [ -e `echo $f | sed s/coll1/coll2/` ]; then return 0; fi
	done
	return 1
}

function doconnect {
	let c=0
	CPIDS=""
	while [ $c -lt $CPUS ]; do
		mkdir $1/connect$c
		$CONNECT -w $1/connect$c -t $TTT --inputfilelow $1/paths$(($TTT-1))_0of1.bin.gz --inputfilehigh $1/paths$(($TTT+4))_0of1.bin.gz -m $CPUS -i $c &
		CPIDS="$CPIDS $!"
		let c=c+1
	done
	let contime=0
	while true; do
		sleep 5
		let contime=contime+$((5*$CPUS))
		cstop=true;
		for cp in $CPIDS; do
			if ps -p $cp 2>/dev/null >/dev/null; then cstop=false; fi
		done
		if $cstop; then break; fi
		if testcoll $1 ; then
			sleep 5
			kill $CPIDS
			break;
		fi
		if [ $contime -gt 7200 ]; then
			kill $CPIDS
			break;
		fi
	done
}

function docollfind {
	cf=""
	for f in `ls -t $2/bestpath*.bin.gz`; do	
		$HELPER -w $1 --findcoll $f 2>&1 > $2/helper.log &
		pid=$!
		sleep 10
		kill $pid
		if [ `grep "^[.].*$" $2/helper.log | wc -l` -gt 0 ]; then cf=$f; break; fi
	done
	$HELPER -w $1 --findcoll $cf &
	CPID=$!
	while true; do
		if testcoll $1 ; then break; fi
		sleep 1
	done
	sleep 5
	kill $CPID
}

cp file1.bin file1_0.bin
cp file2.bin file2_0.bin
let k=0$3

while true; do
	rm -r workdir$k 2>/dev/null
	mkdir workdir$k
	$HELPER -w workdir$k --startnearcollision file1_$k.bin file2_$k.bin --pathtyperange 2 || exit
	cp *.cfg workdir$k

	if [ $CPUS -gt 1 ]; then
		doforward workdir$k &
		fpid=$!
		dobackward workdir$k &
		bpid=$!
		while true; do
			if ps -p $fpid 2>/dev/null >/dev/null; then continue; fi
			if ps -p $bpid 2>/dev/null >/dev/null; then continue; fi
			break;
		done
	else
		doforward workdir$k
		dobackward workdir$k
	fi

	doconnect workdir$k

	for d in workdir$k/connect*; do
		docollfind workdir$k $d &
	done
	while true; do
		if testcoll workdir$k ; then break; fi
		sleep 1
	done
	sleep 5

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
