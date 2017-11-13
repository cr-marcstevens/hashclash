#!/bin/bash

N=2
tconnect=12
tmid=32
data=200000
prefixfile=prefix.bin

if [ ! -z $1 ]; then
	prefixfile="$1"
fi
if [ ! -z $prefixfile ]; then
	if [ ! -f $prefixfile ]; then
		echo "Not using prefixfile: doesn't exist"
		pause
		prefixfile=
	fi
fi

case $N in 

1)
diff="--diffm2 9"
;;

2)
diff="--diffm13 28 --diffm0 32 --diffm6 32"
;;

3)
diff="--diffm6 9 --diffm9 32 --diffm15 32"
tconnect=10
data=500000
;;

esac

opts="-a $data -e 4 --fillfraction 1 -q 8"

# automatically kill connect 10 seconds after first found differential path
function auto_kill_connect
{
	while [ ! -f data/bestpath.bin.gz ]; do
		sleep 1
	done
	sleep 10
	killall md5_diffpathconnect
}

rm *.log
rm -r data
mkdir -p data

#create upper paths
../bin/md5_diffpathbackward $diff -n -t $tmid -b $(($tmid+1))
../bin/md5_diffpathforward $diff $opts -t $((tmid+1)) -b $tmid --trange $((61-$tmid)) 2>&1 | tee upper.log
../bin/md5_diffpathforward $diff -a 16 -t 63 -b $tmid 2>&1 | tee -a upper.log

../bin/md5_diffpathbackward $diff $opts -f data/paths63_0of1.bin.gz -t $((tmid-1)) --trange $(($tmid-$tconnect-5)) 2>&1 | tee -a upper.log



###################### FIRST NEAR-COLLISION BLOCK PAIR #####################

#create lower paths
if [ -z $prefixfile ]; then
cat <<EOF > path_cv.txt
Q-3:    |........ ........ ........ ........| ok p=1
Q-2:    |........ ........ ........ ........| ok p=1
Q-1:    |........ ........ ........ ........| ok p=1
Q0:     |........ ........ ........ ........| ok p=1
EOF
../bin/md5_diffpathhelper $diff --pathfromtext --inputfile1 path_cv.txt --outputfile1 path_cv.bin
../bin/md5_diffpathforward $diff $opts -f path_cv.bin --normalt01 -t 0 --trange $(($tconnect-1)) 2>&1 | tee -a lower.log
else
../bin/md5_diffpathhelper $diff --startpartialpathfromfile --inputfile1 $prefixfile --outputfile1 $prefixfile.base --outputfile2 path_prefix.bin 2>&1 | tee lower.log
startt=`cat lower.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
../bin/md5_diffpathforward $diff $opts -f path_prefix.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a lower.log
fi

#./md5diffpathbackward $diff -t 0 -b 0 -n 2>&1 | tee lower.log
#./md5diffpathforward $diff $opts --normalt01 -t 1 --trange $(($tconnect-2)) 2>&1 | tee -a lower.log

#connect
auto_kill_connect &
../bin/md5_diffpathconnect $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh data/paths$((tconnect+4))_0of1.bin.gz 2>&1 | tee connect.log
touch data/bestpath.bin.gz
sleep 2

#find example collision
../bin/md5_diffpathhelper $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee collfind.log

################## SECOND NEAR-COLLISION BLOCK PAIR #########################

# swap msg1 and msg2 to use opposing signs as first block
cat $prefixfile.base ./data/coll1_* > ./data/file2.bin
cat $prefixfile.base ./data/coll2_* > ./data/file1.bin
mv ./data/coll1_* ./data/nc1_coll1.bin
mv ./data/coll2_* ./data/nc1_coll2.bin
mv data/bestpath.bin.gz data/nc1_bestpath.bin.gz

../bin/md5_diffpathhelper $diff --startpartialpathfromfile --inputfile1 ./data/file1.bin --inputfile2 ./data/file2.bin --outputfile1 /dev/null --outputfile2 path_2nc.bin 2>&1 | tee lower2.log
startt=`cat lower2.log | grep "^Q" | tail -n1 | cut -d: -f1 | cut -dQ -f2`
../bin/md5_diffpathforward $diff $opts -f path_2nc.bin --normalt01 -t $startt --trange $(($tconnect-$startt-1)) 2>&1 | tee -a lower2.log

#connect
auto_kill_connect &
../bin/md5_diffpathconnect $diff -t $tconnect --inputfilelow data/paths$((tconnect-1))_0of1.bin.gz --inputfilehigh data/paths$((tconnect+4))_0of1.bin.gz 2>&1 | tee connect2.log
touch data/bestpath.bin.gz

#find example collision
../bin/md5_diffpathhelper $diff --findcollision --inputfile1 data/bestpath.bin.gz 2>&1 | tee collfind2.log

cat data/file1.bin data/coll1_* > collision1.bin
cat data/file2.bin data/coll2_* > collision2.bin
md5sum collision1.bin collision2.bin
sha1sum collision1.bin collision2.bin
ls -als collision1.bin collision2.bin
