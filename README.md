# hashclash
Project HashClash - MD5 & SHA-1 cryptanalytic toolbox

## Requirements

- C++11 compiler (e.g. g++)
- make
- autoconf & automake & libtool

  `sudo apt-get install autoconf automake libtool`
  
- zlib & bzip2 libraries

  `sudo apt-get install zlib1g-dev libbz2-dev`
  
- local boost C++ libraries (preferable version 1.57)

  `./install_boost.sh`

## Building

	autoreconf --install
	./configure --with-boost=$PWD/local
	make
