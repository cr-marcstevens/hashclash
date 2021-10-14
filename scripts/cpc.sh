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
export EXPECTED_BLOCKTIME=$((3600 * 8 / $CPUS))
export AUTO_KILL_TIME=$((EXPECTED_BLOCKTIME * 2))

rm -r data 2>/dev/null
mkdir data || exit 1

file1=$1
file2=$2
starttime=$(date +%s)

function notify {
	msg=$1
	echo "[*] $msg"
	if [[ -x "$(command -v notify-send)" ]]; then
		notify-send "HashClash" "$msg"
	fi
}

if [ -f "$file1" -a -f "$file2" ]; then
	echo "Chosen-prefix file 1: $file1"
	echo "Chosen-prefix file 2: $file2"
else
	echo "Usage $0 <file1> <file2> [<redonearcollstep>] [<nobirthday>]"
	exit 1
fi

if [ "$3" = "" ] ; then
	# check for CUDA / AVX256 threads
	echo -n "Detecting worklevel..."
	worklevel=`$BIRTHDAYSEARCH --inputfile1 "$file1" --inputfile2 "$file2" --hybridbits 0 --pathtyperange 2 --maxblocks 9 --maxmemory 100 --threads $CPUS --cuda_enable |& grep "^Work" -m1 | head -n1 | cut -d'(' -f2 | cut -d'.' -f1`
	echo ": $worklevel"
	if [ $worklevel -ge 31 ]; then
		$BIRTHDAYSEARCH --inputfile1 "$file1" --inputfile2 "$file2" --hybridbits 0 --pathtyperange 2 --maxblocks 7 --maxmemory 4000 --threads $CPUS --cuda_enable
	else
		$BIRTHDAYSEARCH --inputfile1 "$file1" --inputfile2 "$file2" --hybridbits 0 --pathtyperange 2 --maxblocks 9 --maxmemory 100 --threads $CPUS --cuda_enable
	fi
	notify "Birthday search completed."
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
			kill $CPID &>/dev/null
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
	kill $CPID &>/dev/null
}

function dostepk {
	k=$1
	pid=$$
	echo "[*] Starting step $k"
	sleep 1

	rm -r workdir$k 2>/dev/null
	mkdir workdir$k

	echo $$ > workdir$k/pid
	set -o pipefail
	$HELPER -w workdir$k --startnearcollision file1_$k.bin file2_$k.bin --pathtyperange 2 |& tee workdir$k/start.log
	if [[ $? -ne 0 ]]; then
	    touch workdir$k/killed
	    exit
	fi

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
}

function auto_kill
{
	workdir="$1"
	pidfile="$workdir/pid"
	killfile="$workdir/killed"
	remaining=$2

	while [[ $remaining -gt 0 ]]; do
		echo "[*] Time before backtrack: $remaining s"
		sleep 10
		let remaining=$remaining-10
	done

	pid=$(<$pidfile)
	echo "[*] Timeout reached. Killing process with pid $pid"
	touch $killfile
	pkill -KILL -P $pid &>/dev/null
	kill -KILL $pid &>/dev/null
}

cp file1.bin file1_0.bin
cp file2.bin file2_0.bin
let k=0
if [ "$3" != "" ]; then
	let k=$3
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

backtracks=0
while true; do

	echo "[*] Number of backtracks until now: $backtracks"
	if [[ $backtracks -gt 20 ]]; then
		notify "More than 20 backtracks is not normal. Please restart from scratch."
		break
	fi


	# Check if the collision has been generated
	if [ -f ${file1}.coll -a -f ${file2}.coll ]; then
		if [ `cat ${file1}.coll | md5sum | cut -d' ' -f1` = `cat ${file2}.coll | md5sum | cut -d' ' -f1` ]; then
			notify "Collision generated: ${file1}.coll ${file2}.coll"
			md5sum ${file1}.coll ${file2}.coll
			break
		fi
	fi

	# Start the autokiller
	auto_kill workdir$k $AUTO_KILL_TIME &
	autokillerpid=$!
	mainpid=$$
	trap "kill $autokillerpid &>/dev/null; pkill -KILL -P $mainpid &>/dev/null; killall -r md5_ &>/dev/null; exit" TERM INT KILL

	# Start the computation
	rm -f step$k.log
	(dostepk $k 2>&1 &) | tee step$k.log
	kill $autokillerpid &>/dev/null
	killall -r md5_ &>/dev/null

	# Check if the termination was completed or killed
	if [[ -f workdir$k/killed ]]; then
		failedk=$k
		let k=$((k > 1 ? k-1 : 0))
		notify "Step $failedk failed. Backtracking to step $k"
		let backtracks=backtracks+1
	else
		notify "Step $k completed"
		let k=k+1
	fi
	sleep 2
done

runtime=$((($(date +%s)-$starttime)/60))
notify "Process completed in $runtime minutes ($backtracks backtracks)."

# kill any pending thing
rm md5diffpath*.cfg md5diffpath*.template
pkill -P $$ &>/dev/null
killall -r md5_ &>/dev/null
exit
