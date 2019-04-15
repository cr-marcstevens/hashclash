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

/* Use one of the following typedef's */
//#include <stdint.h>
//typedef ::uint32_t uint32;
typedef unsigned int uint32;
//typedef unsigned __int32 uint32;

void md5_compress(uint32 ihv[], const uint32 block[]);

void find_block0(uint32 block[], const uint32 IV[]);
void find_block1(uint32 block[], const uint32 IV[]);

// very fast inlined xorshift random number generator with period 2^64 - 1
// by G. Marsaglia: http://www.jstatsoft.org/v08/i14/xorshift.pdf 
extern uint32 seed32_1, seed32_2;
inline uint32 xrng64()
{
	uint32 t = seed32_1 ^ (seed32_1 << 10);
	seed32_1 = seed32_2;
	seed32_2 = (seed32_2^(seed32_2>>10))^(t^(t>>13));
	return seed32_1;
}

inline uint32 FF(uint32 b, uint32 c, uint32 d) 
{	return d ^ (b & (c ^ d)); }
inline uint32 GG(uint32 b, uint32 c, uint32 d) 
{	return c ^ (d & (b ^ c)); }
inline uint32 HH(uint32 b, uint32 c, uint32 d) 
{	return b ^ c ^ d; }
inline uint32 II(uint32 b, uint32 c, uint32 d) 
{	return c ^ (b | ~d); }
inline uint32 RL(uint32 x, unsigned int n)
{	return (x << n) | (x >> (32-n)); }
inline uint32 RR(uint32 x, unsigned int n)
{	return (x >> n) | (x << (32-n)); }

#define MD5_STEP(f, a, b, c, d, m, ac, rc) ( \
	a += f(b, c, d) + m + ac, \
	a = (a<<rc | a>>(32-rc)) + b )
#define MD5_REVERSE_STEP(t,AC,RC) (	\
	block[t] = Q[Qoff + t + 1] - Q[Qoff + t], \
	block[t] = RR(block[t], RC) - FF(Q[Qoff + t], Q[Qoff + t - 1], Q[Qoff + t - 2]) - Q[Qoff + t - 3] - AC )

#define Qoff 3
