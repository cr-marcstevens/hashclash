#!/bin/bash

if [ ! -f Makefile ]; then
	echo "First run ./configure [--with-cuda=/usr/local/cuda-x.y]"
	exit 1
fi

if [ "$1" = "" ]; then
	echo "Usage: $0 <release-name>"
	echo "e.g.:  $0 release-linux-v0"
	exit 1
fi
d="$1"

rm -r "$d"

cat Makefile \
	| sed "s/^\(LDFLAGS =.*$\)/\1 -static -static-libstdc++/g" \
	| sed "s/ -lcudart / -lcudart_static -ldl -lrt /" \
	| sed "s/ -lcuda[ ]*$/ -Wl,-Bdynamic,-lcuda/" \
	> Makefile.tmp
	
AM_MAKEFLAGS="-f Makefile.tmp" make -f Makefile.tmp clean
AM_MAKEFLAGS="-f Makefile.tmp" make -f Makefile.tmp -j4 || exit
rm Makefile.tmp

mkdir -p "$d/bin" "$d/scripts"
cp bin/* "$d/bin/"
cp scripts/*.sh "$d/scripts/"

tar -czvf "$d.tar.gz" "$d/"
