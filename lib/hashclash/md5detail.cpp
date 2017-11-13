/**************************************************************************\
|
|    Copyright (C) 2009 Marc Stevens
|
|    This program is free software: you can redistribute it and/or modify
|    it under the terms of the GNU General Public License as published by
|    the Free Software Foundation, either version 3 of the License, or
|    (at your option) any later version.
|
|    This program is distributed in the hope that it will be useful,
|    but WITHOUT ANY WARRANTY; without even the implied warranty of
|    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
|    GNU General Public License for more details.
|
|    You should have received a copy of the GNU General Public License
|    along with this program.  If not, see <http://www.gnu.org/licenses/>.
|
\**************************************************************************/

#include <stdlib.h>

#include "md5detail.hpp"

namespace hashclash {

#ifndef MD5DETAIL_INLINE_IMPL

	const uint32 md5_iv[] = { 0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476 };

	const uint32 md5_ac[] =
		{ 0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501
		, 0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821
		, 0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x2441453 , 0xd8a1e681, 0xe7d3fbc8
		, 0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a
		, 0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70
		, 0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x4881d05 , 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665
		, 0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1
		, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391 };

	const unsigned md5_rc[] = 
		{ 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22
		, 5, 9 , 14, 20, 5, 9 , 14, 20, 5, 9 , 14, 20, 5, 9 , 14, 20
		, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23
		, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21 };

	const unsigned md5_wt[] =
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
		, 1, 6, 11, 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12
		, 5, 8, 11, 14, 1, 4, 7, 10, 13, 0, 3, 6, 9, 12, 15, 2
		, 0, 7, 14, 5, 12, 3, 10, 1, 8, 15, 6, 13, 4, 11, 2, 9 };

#endif


	#define HASHCLASH_MD5COMPRESS_STEP(f, a, b, c, d, m, ac, rc) \
		a += f(b, c, d) + m + ac; a = rotate_left(a,rc); a += b;
	
	void md5compress(uint32 ihv[4], const uint32 block[16])
	{
		uint32 a = ihv[0]; uint32 b = ihv[1]; uint32 c = ihv[2]; uint32 d = ihv[3];

		HASHCLASH_MD5COMPRESS_STEP(md5_ff, a, b, c, d, block[ 0], 0xd76aa478,  7);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, d, a, b, c, block[ 1], 0xe8c7b756, 12); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, c, d, a, b, block[ 2], 0x242070db, 17); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, b, c, d, a, block[ 3], 0xc1bdceee, 22); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, a, b, c, d, block[ 4], 0xf57c0faf,  7);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, d, a, b, c, block[ 5], 0x4787c62a, 12); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, c, d, a, b, block[ 6], 0xa8304613, 17); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, b, c, d, a, block[ 7], 0xfd469501, 22); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, a, b, c, d, block[ 8], 0x698098d8,  7);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, d, a, b, c, block[ 9], 0x8b44f7af, 12); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, c, d, a, b, block[10], 0xffff5bb1, 17);
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, b, c, d, a, block[11], 0x895cd7be, 22);
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, a, b, c, d, block[12], 0x6b901122,  7); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, d, a, b, c, block[13], 0xfd987193, 12);
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, c, d, a, b, block[14], 0xa679438e, 17);
		HASHCLASH_MD5COMPRESS_STEP(md5_ff, b, c, d, a, block[15], 0x49b40821, 22);
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, a, b, c, d, block[ 1], 0xf61e2562,  5);  
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, d, a, b, c, block[ 6], 0xc040b340,  9);  
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, c, d, a, b, block[11], 0x265e5a51, 14);
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, b, c, d, a, block[ 0], 0xe9b6c7aa, 20); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, a, b, c, d, block[ 5], 0xd62f105d,  5);  
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, d, a, b, c, block[10], 0x02441453,  9); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, c, d, a, b, block[15], 0xd8a1e681, 14);
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, b, c, d, a, block[ 4], 0xe7d3fbc8, 20); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, a, b, c, d, block[ 9], 0x21e1cde6,  5);  
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, d, a, b, c, block[14], 0xc33707d6,  9); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, c, d, a, b, block[ 3], 0xf4d50d87, 14); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, b, c, d, a, block[ 8], 0x455a14ed, 20); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, a, b, c, d, block[13], 0xa9e3e905,  5); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, d, a, b, c, block[ 2], 0xfcefa3f8,  9);  
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, c, d, a, b, block[ 7], 0x676f02d9, 14); 
		HASHCLASH_MD5COMPRESS_STEP(md5_gg, b, c, d, a, block[12], 0x8d2a4c8a, 20);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, a, b, c, d, block[ 5], 0xfffa3942,  4); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, d, a, b, c, block[ 8], 0x8771f681, 11); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, c, d, a, b, block[11], 0x6d9d6122, 16);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, b, c, d, a, block[14], 0xfde5380c, 23);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, a, b, c, d, block[ 1], 0xa4beea44,  4);  
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, d, a, b, c, block[ 4], 0x4bdecfa9, 11); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, c, d, a, b, block[ 7], 0xf6bb4b60, 16); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, b, c, d, a, block[10], 0xbebfbc70, 23);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, a, b, c, d, block[13], 0x289b7ec6,  4); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, d, a, b, c, block[ 0], 0xeaa127fa, 11); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, c, d, a, b, block[ 3], 0xd4ef3085, 16); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, b, c, d, a, block[ 6], 0x04881d05, 23); 
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, a, b, c, d, block[ 9], 0xd9d4d039,  4);  
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, d, a, b, c, block[12], 0xe6db99e5, 11);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, c, d, a, b, block[15], 0x1fa27cf8, 16);
		HASHCLASH_MD5COMPRESS_STEP(md5_hh, b, c, d, a, block[ 2], 0xc4ac5665, 23); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, block[ 0], 0xf4292244,  6);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, block[ 7], 0x432aff97, 10); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, block[14], 0xab9423a7, 15);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, block[ 5], 0xfc93a039, 21); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, block[12], 0x655b59c3,  6); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, block[ 3], 0x8f0ccc92, 10); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, block[10], 0xffeff47d, 15);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, block[ 1], 0x85845dd1, 21); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, block[ 8], 0x6fa87e4f,  6);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, block[15], 0xfe2ce6e0, 10);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, block[ 6], 0xa3014314, 15); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, block[13], 0x4e0811a1, 21);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, block[ 4], 0xf7537e82,  6);  
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, block[11], 0xbd3af235, 10);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, block[ 2], 0x2ad7d2bb, 15); 
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, block[ 9], 0xeb86d391, 21); 

		ihv[0] += a; ihv[1] += b; ihv[2] += c; ihv[3] += d;
	}

} // namespace hashclash
