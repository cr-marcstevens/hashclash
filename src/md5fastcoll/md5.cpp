/*

MD5 collision generator
=======================
Source code files:
  block0.cpp
  block1.cpp
  main.cpp
  main.hpp
  md5.cpp
  block1wang.cpp
  block1stevens00.cpp
  block1stevens01.cpp
  block1stevens10.cpp
  block1stevens11.cpp
Win32 executable:
  fastcoll_v1.0.0.5.exe

Version
=======
version 1.0.0.5-1, April 2006.

Copyright
=========
© M. Stevens, 2006. All rights reserved.

Disclaimer
==========
This software is provided as is. Use is at the user's risk.
No guarantee whatsoever is given on how it may function or malfunction.
Support cannot be expected.
This software is meant for scientific and educational purposes only.
It is forbidden to use it for other than scientific or educational purposes.
In particular, commercial and malicious use is not allowed.
Further distribution of this software, by whatever means, is not allowed
without our consent.
This includes publication of source code or executables in printed form,
on websites, newsgroups, CD-ROM's, etc.
Changing the (source) code without our consent is not allowed.
In all versions of the source code this disclaimer, the copyright
notice and the version number should be present.

*/

#include "main.hpp"

void md5_compress(uint32 ihv[], const uint32 block[])
{
	uint32 a = ihv[0];
	uint32 b = ihv[1];
	uint32 c = ihv[2];
	uint32 d = ihv[3];

	MD5_STEP(FF, a, b, c, d, block[0], 0xd76aa478, 7);
	MD5_STEP(FF, d, a, b, c, block[1], 0xe8c7b756, 12);
	MD5_STEP(FF, c, d, a, b, block[2], 0x242070db, 17);
	MD5_STEP(FF, b, c, d, a, block[3], 0xc1bdceee, 22);
	MD5_STEP(FF, a, b, c, d, block[4], 0xf57c0faf, 7);
	MD5_STEP(FF, d, a, b, c, block[5], 0x4787c62a, 12);
	MD5_STEP(FF, c, d, a, b, block[6], 0xa8304613, 17);
	MD5_STEP(FF, b, c, d, a, block[7], 0xfd469501, 22);
	MD5_STEP(FF, a, b, c, d, block[8], 0x698098d8, 7);
	MD5_STEP(FF, d, a, b, c, block[9], 0x8b44f7af, 12);
	MD5_STEP(FF, c, d, a, b, block[10], 0xffff5bb1, 17);
	MD5_STEP(FF, b, c, d, a, block[11], 0x895cd7be, 22);
	MD5_STEP(FF, a, b, c, d, block[12], 0x6b901122, 7);
	MD5_STEP(FF, d, a, b, c, block[13], 0xfd987193, 12);
	MD5_STEP(FF, c, d, a, b, block[14], 0xa679438e, 17);
	MD5_STEP(FF, b, c, d, a, block[15], 0x49b40821, 22);
	MD5_STEP(GG, a, b, c, d, block[1], 0xf61e2562, 5);
	MD5_STEP(GG, d, a, b, c, block[6], 0xc040b340, 9);
	MD5_STEP(GG, c, d, a, b, block[11], 0x265e5a51, 14);
	MD5_STEP(GG, b, c, d, a, block[0], 0xe9b6c7aa, 20);
	MD5_STEP(GG, a, b, c, d, block[5], 0xd62f105d, 5);
	MD5_STEP(GG, d, a, b, c, block[10], 0x02441453, 9);
	MD5_STEP(GG, c, d, a, b, block[15], 0xd8a1e681, 14);
	MD5_STEP(GG, b, c, d, a, block[4], 0xe7d3fbc8, 20);
	MD5_STEP(GG, a, b, c, d, block[9], 0x21e1cde6, 5);
	MD5_STEP(GG, d, a, b, c, block[14], 0xc33707d6, 9);
	MD5_STEP(GG, c, d, a, b, block[3], 0xf4d50d87, 14);
	MD5_STEP(GG, b, c, d, a, block[8], 0x455a14ed, 20);
	MD5_STEP(GG, a, b, c, d, block[13], 0xa9e3e905, 5);
	MD5_STEP(GG, d, a, b, c, block[2], 0xfcefa3f8, 9);
	MD5_STEP(GG, c, d, a, b, block[7], 0x676f02d9, 14);
	MD5_STEP(GG, b, c, d, a, block[12], 0x8d2a4c8a, 20);
	MD5_STEP(HH, a, b, c, d, block[5], 0xfffa3942, 4);
	MD5_STEP(HH, d, a, b, c, block[8], 0x8771f681, 11);
	MD5_STEP(HH, c, d, a, b, block[11], 0x6d9d6122, 16);
	MD5_STEP(HH, b, c, d, a, block[14], 0xfde5380c, 23);
	MD5_STEP(HH, a, b, c, d, block[1], 0xa4beea44, 4);
	MD5_STEP(HH, d, a, b, c, block[4], 0x4bdecfa9, 11);
	MD5_STEP(HH, c, d, a, b, block[7], 0xf6bb4b60, 16);
	MD5_STEP(HH, b, c, d, a, block[10], 0xbebfbc70, 23);
	MD5_STEP(HH, a, b, c, d, block[13], 0x289b7ec6, 4);
	MD5_STEP(HH, d, a, b, c, block[0], 0xeaa127fa, 11);
	MD5_STEP(HH, c, d, a, b, block[3], 0xd4ef3085, 16);
	MD5_STEP(HH, b, c, d, a, block[6], 0x04881d05, 23);
	MD5_STEP(HH, a, b, c, d, block[9], 0xd9d4d039, 4);
	MD5_STEP(HH, d, a, b, c, block[12], 0xe6db99e5, 11);
	MD5_STEP(HH, c, d, a, b, block[15], 0x1fa27cf8, 16);
	MD5_STEP(HH, b, c, d, a, block[2], 0xc4ac5665, 23);
	MD5_STEP(II, a, b, c, d, block[0], 0xf4292244, 6);
	MD5_STEP(II, d, a, b, c, block[7], 0x432aff97, 10);
	MD5_STEP(II, c, d, a, b, block[14], 0xab9423a7, 15);
	MD5_STEP(II, b, c, d, a, block[5], 0xfc93a039, 21);
	MD5_STEP(II, a, b, c, d, block[12], 0x655b59c3, 6);
	MD5_STEP(II, d, a, b, c, block[3], 0x8f0ccc92, 10);
	MD5_STEP(II, c, d, a, b, block[10], 0xffeff47d, 15);
	MD5_STEP(II, b, c, d, a, block[1], 0x85845dd1, 21);
	MD5_STEP(II, a, b, c, d, block[8], 0x6fa87e4f, 6);
	MD5_STEP(II, d, a, b, c, block[15], 0xfe2ce6e0, 10);
	MD5_STEP(II, c, d, a, b, block[6], 0xa3014314, 15);
	MD5_STEP(II, b, c, d, a, block[13], 0x4e0811a1, 21);
	MD5_STEP(II, a, b, c, d, block[4], 0xf7537e82, 6);
	MD5_STEP(II, d, a, b, c, block[11], 0xbd3af235, 10);
	MD5_STEP(II, c, d, a, b, block[2], 0x2ad7d2bb, 15);
	MD5_STEP(II, b, c, d, a, block[9], 0xeb86d391, 21);

	ihv[0] += a;
	ihv[1] += b;
	ihv[2] += c;
	ihv[3] += d;
}
