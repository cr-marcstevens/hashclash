# Project HashClash - MD5 & SHA-1 cryptanalytic toolbox

## Feedback & pingback appreciated!

I'm always curious to hear where HashClash is used for! ([@realhashbreaker](https://twitter.com/realhashbreaker)  on Twitter or let me know via [Issues](https://github.com/cr-marcstevens/hashclash/issues))

If you like HashClash and use it more often then please consider becoming a sponsor.

[![](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/cr-marcstevens)

## Requirements

- C++11 compiler (e.g. g++)
- make
- autoconf & automake & libtool

  `sudo apt-get install autoconf automake libtool`
  
- zlib & bzip2 libraries

  `sudo apt-get install zlib1g-dev libbz2-dev`
  
- Optional: CUDA
  
## Building (automatic)

- Run build.sh script

   `./build.sh`

## Building (manual)

<details>
<summary>See details</summary>
<p>
  
- local boost C++ libraries

  `./install_boost.sh` 

  Override default boost version and/or installation directory as follows:
  
  `BOOST_VERSION=1.88.0 BOOST_INSTALL_PREFIX=$HOME/boost/boost-1.88.0 ./install_boost.sh`
  
- Build configure script

  `autoreconf --install`
  
- Run configure (with boost installed in `$(pwd)/boost-VERSION` by `install_boost.sh`)

  `./configure --with-boost=$(pwd)/boost-1.88.0 [--without-cuda|--with-cuda=/usr/local/cuda-X.X]`

- Build programs

  `make [-j 4]`

</p>
</details>

## Create your own chosen-prefix collisions

- Create temporary working directory (in top directory)

  `mkdir cpc_workdir`
  
  `cd cpc_workdir`
  
- Run script

  `../scripts/cpc.sh <prefix.filename1> <prefix.filename2>`

## Create you own identical-prefix collision

- Create temporary working directory

  `mkdir ipc_workdir`

  `cd ipc_workdir`

- Run script

  `echo -n "TEST" > prefix.txt`

  `../scripts/generic_ipc.sh prefix.txt`

  Note: the prefix file is expected to be a multiple of 64 bytes 
  and optionally plus 1, 2, or 3 bytes or a small multiple of 4 bytes.
  These last bytes will be used as forced message words
  inside the first near-collision block.
  Any remaining 1, 2 or 3 bytes of the prefix file are ignored.

  Note: the first time the script is run it will be slower
  as it needs to create the 'upper' differential path set for the first time.

- Example output:

```
$ xxd collision1.bin
0000000: 5445 5354 1789 6de2 c568 339d 8bbf 6269  TEST..m..h3...bi
0000010: 563a c4ab 2fba 7aba 6ec7 e182 8566 8883  V:../.z.n....f..
0000020: b0f2 d716 17d5 8f97 b244 b9ca dcaa af93  .........D......
0000030: 3a33 4552 a9fd 023b 7012 7e5c d644 9646  :3ER...;p.~\.D.F
0000040: 723e 737a 6c3a 66e5 5d51 8e7c 7bc2 ec4f  r>szl:f.]Q.|{..O
0000050: 95a2 349f 0f4b 2540 cbe6 5644 d113 104b  ..4..K%@..VD...K
0000060: 5c39 898d f19a f9e9 e2aa bbad e191 d2d3  \9..............
0000070: a25d 4df0 9058 a873 4ee8 9dc9 47a5 281c  .]M..X.sN...G.(.
```

- Make your own attack

  Inside generic_ipc.sh there are three examples of identical-prefix collision attacks,
  selected by N=1, 2 or 3.
  You can add your own choice of message word differences here
  and see if you can make your own collision attack!

## Create you own text identical-prefix collision

- Example:

```
md5("TEXTCOLLBYfGiJUETHQ4hAcKSMd5zYpgqf1YRDhkmxHkhPWptrkoyz28wnI9V0aHeAuaKnak")
=
md5("TEXTCOLLBYfGiJUETHQ4hEcKSMd5zYpgqf1YRDhkmxHkhPWptrkoyz28wnI9V0aHeAuaKnak")
```

- Create temporary working directory

  `mkdir textcoll_workdir`

  `cd textcoll_workdir`

- Run script

  `../scripts/textcoll.sh`

- Edit the script for more options: change alphabet, force specific bytes, etc...
