# hashclash
Project HashClash - MD5 & SHA-1 cryptanalytic toolbox

## Requirements

- C++11 compiler (e.g. g++)
- make
- autoconf & automake & libtool

  `sudo apt-get install autoconf automake libtool`
  
- zlib & bzip2 libraries

  `sudo apt-get install zlib1g-dev libbz2-dev`
  
- local boost C++ libraries **compiled with -std=c++11** (preferable version 1.57)

  `./install_boost.sh` 

  Override default boost version 1.57 and/or installation directory as follows:
  
  `BOOST_VERSION=1.65.1 BOOST_INSTALL_PREFIX=$(pwd)/boost-1.65.1 ./install_boost.sh`
  
## Building

- Build configure script

  `autoreconf --install`
  
- Run configure (with boost installed in $(pwd)/local by install_boost.sh)

  `./configure --with-boost=$PWD/local`
  
- Build programs

  `make`

## Create your own chosen-prefix collisions

- Create temporary working directory (in top directory)

  `mkdir cpc_workdir`
  `cd cpc_workdir`
  
- Copy necessary files

  `cp ../scripts/cpc.sh ../scripts/*.template .`
  
- Run script

  `./cpc.sh <prefix.filename1> <prefix.filename2>`

- Monitor progress of script

  `ls -altr`
  
  If last change to any directory is several hours old then
  * kill script & any running md5_diffpathhelper programs
  * let K be the number of the last `workdir$(K)` directory
  * restart script:
    `./cpc.sh <prefix.filename1> <prefix.filename2> $((K-1))`
