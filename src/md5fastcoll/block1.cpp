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

#include <iostream>
#include <vector>
#include "main.hpp"


void find_block1_wang(uint32 block[], const uint32 IV[]);
void find_block1_stevens_11(uint32 block[], const uint32 IV[]);
void find_block1_stevens_10(uint32 block[], const uint32 IV[]);
void find_block1_stevens_01(uint32 block[], const uint32 IV[]);
void find_block1_stevens_00(uint32 block[], const uint32 IV[]);

void find_block1(uint32 block[], const uint32 IV[])
{
	if ( ((IV[1]^IV[2])&(1<<31))==0 && ((IV[1]^IV[3])&(1<<31))==0
		&& (IV[3]&(1<<25))==0 && (IV[2]&(1<<25))==0 && (IV[1]&(1<<25))==0 && ((IV[2]^IV[1])&1)==0 
	   )
	{
		uint32 IV2[4] = { IV[0]+(1<<31), IV[1]+(1<<31)+(1<<25), IV[2]+(1<<31)+(1<<25), IV[3]+(1<<31)+(1<<25) };
		if ((IV[1]&(1<<6))!=0 && (IV[1]&1)!=0) {
			std::cout << "S11" << std::flush;
			find_block1_stevens_11(block, IV2);
		} else if ((IV[1]&(1<<6))!=0 && (IV[1]&1)==0) {
			std::cout << "S10" << std::flush;
			find_block1_stevens_10(block, IV2);
		} else if ((IV[1]&(1<<6))==0 && (IV[1]&1)!=0) {
			std::cout << "S01" << std::flush;
			find_block1_stevens_01(block, IV2);
		} else {
			std::cout << "S00" << std::flush;
			find_block1_stevens_00(block, IV2);
		}
		block[4] += 1<<31;
		block[11] += 1<<15;
		block[14] += 1<<31;
	} else {
		std::cout << "W" << std::flush;
		find_block1_wang(block, IV);
	}
}
