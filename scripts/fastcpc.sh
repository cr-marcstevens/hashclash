#!/bin/bash

PATHTYPERANGE=2
HYBRIDBITS=0
MAXBLOCKS=7
MAXMEMORY=8000
UPPERPATHCOUNT=1000000
TTT=12

# Block timeout: Once a near-collision block lasts more than this number of second
#                the whole attack is aborted
# A lower value means less time lost once a near-collision block gets 'stuck',
# but also implies a lower probability of success.
#
# Default value is for CPU i7-11800H @ 2.30GHz (8 cores / 16 threads): 2000s
# To fine-tune BLOCK_TIMEOUT for a particular machine use the following procedure:
#   1) set BLOCK_TIMEOUT to 3600 or larger
#   2) run 20 attempts and save logs:
#      $> for ((i=0; i<20; ++i)); do ./fastcpc.sh ../README.md ../LICENSE.txt > cpc$i.log 2>&1; done
#   3) Parse logs for near-collision block durations:
#      $> grep "NC-BLOCK.*FINISHED" cpc*/*.log | cut -d: -f3 | cut -ds -f1 | sort -g
#   4) Choose appropriate cut-off timeout at the end
BLOCK_TIMEOUT=2000

VERBOSE_LOG=2

BIRTHDAY_CONTROLLER_NODES=1
BIRTHDAY_GENERATOR_NODES=0

BIRTHDAY_LOCAL_CONFIG="--cuda_enable --saveloadwait 1"
BIRTHDAY_GLOBAL_CONFIG="--pathtyperange $PATHTYPERANGE --hybridbits $HYBRIDBITS --maxblocks $MAXBLOCKS --maxmemory $MAXMEMORY"

# this is the workdir in which the script was called
WORKDIR=`pwd`

# this is the directory of the script itself within 'hashclash/scripts/'
SCRIPTDIR=$(dirname $0)

# from SCRIPTDIR we derive the 'hashclash/bin' directory where all hashclash binaries reside
HASHCLASHBIN=$SCRIPTDIR/../bin

function logthis
{
	loglvl=$1
	logstep=$2
	shift 2
	if [ "$CPCTIMESTAMP" = "" ]; then
		CPCTIMESTAMP=`date +%s`
		DATESTR=`date -Is`
	else
		DATESTR=+$(( `date +%s` - $CPCTIMESTAMP ))s
	fi
	if [ $loglvl -le $VERBOSE_LOG ]; then
		echo "[LVL$loglvl][$DATESTR][$logstep] $*"
	fi
}

function precompute
{
	PRECOMPDIR=precompute/diff$1
	UPPERPATHFILE=$PRECOMPDIR/paths$(($TTT+4))_0of1_$UPPERPATHCOUNT.bin.gz
	if [ ! -f $UPPERPATHFILE ]; then
		logthis 1 "PRECOMPUTE $1" "START"
		TIMESTART=`date +%s`
		mkdir -p $PRECOMPDIR 2>/dev/null
		for f in md5diffpathbackward.cfg md5diffpathforward.cfg md5diffpathhelper.cfg; do
			rm $f 2>/dev/null
			cat ${f}.template > $f 2>/dev/null
			echo "diffm11 = $1" >> $f
			echo "workdir = $PRECOMPDIR" >> $f
		done
		[ ! -f $PRECOMPDIR/paths30_0of1.bin.gz ] && ( $HASHCLASHBIN/md5_diffpathbackward -n -t 34 --trange 4 -a 65536 -q 128 > $PRECOMPDIR/backward.log 2>&1 )
		if [ ! -f $PRECOMPDIR/paths21_0of1.bin.gz ]; then
			for ((x=40; ; ++x)); do
				$HASHCLASHBIN/md5_diffpathbackward -t 29 --trange 8 --maxQ26upcond $x >> $PRECOMPDIR/backward.log 2>&1
				if [ -f $PRECOMPDIR/paths21_0of1.bin.gz ]; then
					logthis 2 "PRECOMPUTE $1" "Using maxQ26upcond = $x"
					break
				fi
			done
		fi
		[ ! -f $PRECOMPDIR/paths20_0of1.bin.gz ] && $HASHCLASHBIN/md5_diffpathbackward -t 20 -q 1 -a 16384 >> $PRECOMPDIR/backward.log 2>&1
		[ ! -f $PRECOMPDIR/paths$(($TTT+5))_0of1.bin.gz ] && $HASHCLASHBIN/md5_diffpathbackward -t 19 --trange $((18-$TTT-4)) >> $PRECOMPDIR/backward.log 2>&1
		if $HASHCLASHBIN/md5_diffpathbackward -t $(($TTT+4)) -a $UPPERPATHCOUNT >> $PRECOMPDIR/backward.log 2>&1 ; then
			mv $PRECOMPDIR/paths$(($TTT+4))_0of1.bin.gz $UPPERPATHFILE
			TIMELAPSED=$((`date +%s` - $TIMESTART))
			logthis 1 "PRECOMPUTE $1" "DONE: $TIMELAPSED s"
		else
			TIMELAPSED=$((`date +%s` - $TIMESTART))
			logthis 1 "PRECOMPUTE $1" "FAILED: $TIMELAPSED s"
			rm $PRECOMPDIR/paths*.bin.gz
		fi
	fi
	if $2 ; then
		PRECOMPDIR2=precompute/diff-$1
		if [ ! -f $PRECOMPDIR2/paths$((TTT+4))_0of1_$UPPERPATHCOUNT.bin.gz ]; then
			logthis 1 "PRECOMPUTE -$1" "NEGATE FROM $1"
			rm -r $PRECOMPDIR2 2>/dev/null
			mkdir $PRECOMPDIR2
			$HASHCLASHBIN/md5_diffpathhelper --negatepaths --inputfile1 $PRECOMPDIR/paths$((TTT+4))_0of1_$UPPERPATHCOUNT.bin.gz --outputfile1 $PRECOMPDIR2/paths$((TTT+4))_0of1_$UPPERPATHCOUNT.bin.gz > $PRECOMPDIR2/negate.log 2>&1
		fi
	fi
}

function allprecompute
{
	logthis 0 "PRECOMPUTE" "Precomputing missing precomputation files"
	for ((diff=1; diff < 32; ++diff)); do
		precompute $diff true
	done
	precompute $diff false
}

function asyncstartprocessonhost
{
	nodeidx=$1
	logfile=$2
	shift 2

	logthis 2 "START PROCESS ON NODE$nodeidx" "$* > $logfile"

	if [ $nodeidx -gt 0 ]; then
		echo "To run this script on multiple nodes you need to manually implement job distribution and kill in this script."
		abort 1
	fi

	# Simply run the command on the local machine
	$* > $logfile 2>&1 &
}
function asynckillprocessonhosts
{
	# Run this command on all hosts: killall $1 (NOT IMPLEMENTED)
	
	# Simply run the kill command on the local machine
	logthis 2 "KILL PROCESS ON NODE0" "$1"
	killall $1
}

function asyncstartbirthdaysearch
{
	mkdir birthday_workdir 2>/dev/null
	rm birthday_workdir/* 2>/dev/null

	logthis 1 "DISTRIBUTING BIRTHDAY" "$HASHCLASHBIN/md5_birthdaysearch -w birthday_workdir -i <nodeidx> -m $BIRTHDAY_CONTROLLER_NODES $BIRTHDAY_LOCAL_CONFIG $BIRTHDAY_GLOBAL_CONFIG $*"
	
	# example single process on a single machine (uses all CUDA devices and all threads)
	for (( i=0 ; i < $BIRTHDAY_CONTROLLER_NODES ; ++i )); do
		asyncstartprocessonhost $i birthday_workdir/birthday_con_node$i.log $HASHCLASHBIN/md5_birthdaysearch -w birthday_workdir -i $i -m $BIRTHDAY_CONTROLLER_NODES $BIRTHDAY_LOCAL_CONFIG $BIRTHDAY_GLOBAL_CONFIG $*
	done
	for (( j=0 ; j < $BIRTHDAY_GENERATOR_NODES ; ++j )); do
		let i=$j+$BIRTHDAY_CONTROLLER_NODES
		asyncstartprocessonhost $i birthday_workdir/birthday_gen_node$i.log $HASHCLASHBIN/md5_birthdaysearch -w birthday_workdir --generatormode -m $BIRTHDAY_CONTROLLER_NODES $BIRTHDAY_LOCAL_CONFIG $BIRTHDAY_GLOBAL_CONFIG $*
	done
}
function asyncstopbirthdaysearch
{
	asynckillprocessonhosts md5_birthdaysearch
}

if [ "$1" = "" ]; then
	echo "Usage: $0 <prefixfile1> <prefixfile2>"
	exit 0
fi
if [ "$1" = "__birthday" ]; then
	shift 1
	startbirthdaysearch $*
	exit 0
fi
if [ "$1" = "__forwardconnect" ]; then
	shift 1
	startforwardconnect $*
	exit 0
fi

PREFIXFILE1=$1
PREFIXFILE2=$2

if [ ! -f $PREFIXFILE1 ]; then
	echo "Cannot find $PREFIXFILE1"
	exit 1
fi
if [ ! -f $PREFIXFILE2 ]; then
	echo "Cannot find $PREFIXFILE2"
	exit 1
fi

logthis 0 "START CHOSEN-PREFIX COLLISION ATTACK"
logthis 0 "WRITING CONFIGURATION FILES"

cat <<EOF > md5diffpathbackward.cfg.template
autobalance = 4000000
estimate = 4
fillfraction = 1
maxsdrs = 1
condtend = 35
EOF

cat <<EOF > md5diffpathforward.cfg.template
autobalance = 3000000
estimate = 2
fillfraction = 1
maxsdrs = 128
minQ456tunnel = 18
minQ91011tunnel = 18
EOF

cat <<EOF > md5diffpathconnect.cfg.template
Qcondstart = 21
waitinputfile = true
# Ignore low quality full paths
mintunnel=35
EOF

allprecompute

function kill_all_processes
{
asynckillprocessonhosts md5_birthdaysearch
asynckillprocessonhosts md5_diffpathhelper
asynckillprocessonhosts md5_diffpathconnect
asynckillprocessonhosts md5_diffpathforward
}
kill_all_processes
trap kill_all_processes EXIT

logthis 0 "START BIRTHDAY" "..."
TIMESTART=`date +%s`

asyncstartbirthdaysearch --inputfile1 $PREFIXFILE1 --inputfile2 $PREFIXFILE2

while true; do
	sleep 1
	FILE1=`ls birthday_workdir/birthdayblock1_*.bin 2>/dev/null| head -n1`
	if [ "$FILE1" = "" ]; then continue; fi
	FILE2=`echo $FILE1 | sed "s/birthdayblock1/birthdayblock2/"`
	if [ -f "$FILE1" -a -f "$FILE2" ]; then
		TIMELAPSED=$((`date +%s` - $TIMESTART))
		logthis 0 "BIRTHDAY DONE" "FINISHED: $TIMELAPSED s"
		break
	fi
done

asyncstopbirthdaysearch

cat file1.bin $FILE1 > file1_0.bin
cat file2.bin $FILE2 > file2_0.bin

for (( i=0; ; ++i)); do
	if [ $i -eq $MAXBLOCKS ]; then
		logthis 0 "NC-BLOCK $i" "Error: exceeded maxblocks ?!?!?"
		exit 1
	fi

	logthis 0 "NC-BLOCK $i" "Prepare near-collision files"
	TIMESTART=`date +%s`
	
	mkdir workdir$i 2>/dev/null
	rm workdir$i/* 2>/dev/null
	
	$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --startnearcollision --pathtyperange $PATHTYPERANGE --inputfile1 file1_$i.bin --inputfile2 file2_$i.bin > workdir$i/snc.log 2>&1
	cp *.cfg workdir$i/

	DIFF=`grep diffm11 md5diffpathhelper.cfg | cut -d= -f2 | tr -d ' '`
	UPPERPATH=precompute/diff$DIFF/paths$(($TTT+4))_0of1_$UPPERPATHCOUNT.bin.gz
	logthis 2 "NC-BLOCK $i" "DIFF=$DIFF UPPERPATH=$UPPERPATH"
	if [ ! -f $UPPERPATH ]; then
		echo "Error: Precomputed upperpath file does not exist:\n $UPPERPATH"
		exit 1
	fi

	logthis 0 "NC-BLOCK $i" "Starting connect process in background ..."
	HALFTHREADS=$((`nproc`/2))
	logthis 2 "NC-BLOCK $i" "Running command: $HASHCLASHBIN/md5_diffpathconnect --threads $HALFTHREADS -w workdir$i -t $TTT --inputfilelow workdir$i/paths$(($TTT-1))_0of1.bin.gz.done --inputfilehigh $UPPERPATH > workdir$i/connect.log &"
	INPUTFILELOW=workdir$i/paths$(($TTT-1))_0of1.bin.gz
	$HASHCLASHBIN/md5_diffpathconnect --threads $HALFTHREADS -w workdir$i -t $TTT --inputfilelow $INPUTFILELOW.done --inputfilehigh $UPPERPATH > workdir$i/connect.log 2>&1 &

	logthis 0 "NC-BLOCK $i" "Forward phase ..."
	logthis 2 "NC-BLOCK $i" "Running command: $HASHCLASHBIN/md5_diffpathforward -w workdir$i -f workdir$i/lowerpath.bin.gz --normalt01 -t 1 --trange $(($TTT-2)) > workdir$i/forward.log"
	$HASHCLASHBIN/md5_diffpathforward -w workdir$i -f workdir$i/lowerpath.bin.gz --normalt01 -t 1 --trange $(($TTT-2)) > workdir$i/forward.log 2>&1
	if [ ! -f $INPUTFILELOW ]; then
		logthis 0 "NC-BLOCK $i" "Error: $INPUTFILELOW not found, forward failed?"
		exit 1
	fi
	mv $INPUTFILELOW $INPUTFILELOW.done

	logthis 0 "NC-BLOCK $i" "Wait for full differential paths from connect"
	let w=0
	while [ ! -f workdir$i/bestpaths.bin.gz ]; do
		sleep 1
		let w=w+1
		if [ $w -gt ${BLOCK_TIMEOUT} ]; then
			logthis 0 "NC-BLOCK $i" "Time-out! Aborting ..."
			killall md5_diffpathconnect
			exit 1
		fi
	done
	
	PATHQUALITY=`cat workdir$i/connect.log | grep "tottunnel" | tail -n1 | grep -o "tottunnel=[0-9]*, totcond=[0-9]*"`
	cp workdir$i/bestpaths.bin.gz workdir$i/latest_bestpaths.bin.gz
	logthis 0 "NC-BLOCK $i" "Starting collision finding: $PATHQUALITY"
	$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --combinepaths --outputfile1 workdir$i/collfindpaths.bin.gz --inputfile1 workdir$i/latest_bestpaths.bin.gz --inputfile2 workdir$i/upperpath.bin.gz >> workdir$i/collfindpath.log 2>&1
	$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --inputfile1 workdir$i/collfindpaths.bin.gz --findcollision --threads $HALFTHREADS >> workdir$i/collfind.log 2>&1 &
	let lastw=0
	for (( ; ; ++w )); do
		sleep 1
		if [ $w -gt ${BLOCK_TIMEOUT} ]; then
			logthis 0 "NC-BLOCK $i" "Time-out! Aborting ..."
			killall md5_diffpathconnect
			killall md5_diffpathhelper
			exit 1
		fi
		if [ `ls workdir$i/coll1_* 2>/dev/null | wc -l` -ne 0 ]; then
			logthis 0 "NC-BLOCK $i" "Collision found"
			killall md5_diffpathconnect
			killall md5_diffpathhelper
			break
		fi
		if [ $(($w - $lastw)) -ge 10 ]; then
		if [ workdir$i/bestpaths.bin.gz -nt workdir$i/latest_bestpaths.bin.gz ]; then
			let lastw=$w
			PATHQUALITY=`cat workdir$i/connect.log | grep "tottunnel" | tail -n1 | grep -o "tottunnel=[0-9]*, totcond=[0-9]*"`
			cp workdir$i/bestpaths.bin.gz workdir$i/latest_bestpaths.bin.gz
			JOIN="-j workdir$i/latest_bestpaths.bin.gz"
			for f in workdir$i/bestpaths_t3[789]* workdir$i/bestpaths_t[456]* ; do
				if [ -f $f ]; then
					JOIN="$JOIN -j $f"
				fi
			done
			$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --outputfile1 workdir$i/joint_bestpaths.bin.gz $JOIN >> workdir$i/join.log 2>&1
			$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --combinepaths --outputfile1 workdir$i/collfindpaths.bin.gz --inputfile1 workdir$i/joint_bestpaths.bin.gz --inputfile2 workdir$i/upperpath.bin.gz >> workdir$i/collfindpath.log 2>&1
			killall md5_diffpathhelper
			logthis 0 "NC-BLOCK $i" "Starting collision finding: $PATHQUALITY"
			$HASHCLASHBIN/md5_diffpathhelper -w workdir$i --inputfile1 workdir$i/collfindpaths.bin.gz --findcollision --threads $HALFTHREADS >> workdir$i/collfind.log 2>&1 &
		fi fi
	done
	
	TIMELAPSED=$((`date +%s` - $TIMESTART))
	logthis 0 "NC-BLOCK $i" "FINISHED: $TIMELAPSED s"
	
	COLLFILE1=`ls workdir$i/coll1_* | head -n1`
	COLLFILE2=`echo $COLLFILE1 | sed "s/coll1_/coll2_/"`
	OUTFILE1=file1_$(($i+1)).bin
	OUTFILE2=file2_$(($i+1)).bin

	cat file1_$i.bin $COLLFILE1 > $OUTFILE1
	cat file2_$i.bin $COLLFILE2 > $OUTFILE2

	if [ `cat $OUTFILE1 | md5sum | cut -d' ' -f1` = `cat $OUTFILE2 | md5sum | cut -d' ' -f1` ]; then
		logthis 0 "FINISHED" "Found collision: $OUTFILE1 $OUTFILE2"
		exit 0
	fi
done
