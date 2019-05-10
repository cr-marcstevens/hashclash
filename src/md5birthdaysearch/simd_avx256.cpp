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

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdexcept>
#include <boost/cstdint.hpp>

using namespace std;

typedef boost::uint32_t uint32;
typedef boost::uint64_t uint64;

#include <hashclash/config.h>
#include "birthday_types.hpp"

#ifndef HASHCLASH_HAVE_AVX2
bool simd_device_avx256::init(const uint32 ihv1b[4], const uint32 ihv2b[4], const uint32 ihv2modb[4], const uint32 precomp1b[4], const uint32 precomp2b[4], const uint32 msg1b[16], const uint32 msg2b[16], uint32 hmask, uint32 dpmask, uint32 maxlen)
{
	return false;
}
void simd_device_avx256::fill_trail_buffer(uint64 seed, vector<trail_type>& buf, bool mod)
{
}
#else

#define SHA1DC_HAVE_AVX256
#include <hashclash/simd/simd_avx256.h>

union simd_word_t {
	SIMD_WORD v;
	uint32 w[SIMD_VECSIZE];
};

struct simd_avx256_detail
{
	SIMD_WORD msg1[16];
	SIMD_WORD msg2[16];
	SIMD_WORD ihv1[4];
	SIMD_WORD ihv2[4];
	SIMD_WORD ihv2mod[4];
	SIMD_WORD precomp1[4];
	SIMD_WORD precomp2[4];
	SIMD_WORD hybridmask, distinguishedpointmask, maximumpathlength;
	uint32 hmask, dpmask, maxlen;
	
	simd_word_t s0, s1, s2, len, e0, e1, e2;
};

bool simd_device_avx256::init(const uint32 ihv1b[4], const uint32 ihv2b[4], const uint32 ihv2modb[4], const uint32 precomp1b[4], const uint32 precomp2b[4], const uint32 msg1b[16], const uint32 msg2b[16], uint32 hmask, uint32 dpmask, uint32 maxlen)
{
	detail = new simd_avx256_detail;
	for (unsigned i = 0; i < 16; ++i)
	{
		detail->msg1[i] = SIMD_WTOV(msg1b[i]);
		detail->msg2[i] = SIMD_WTOV(msg2b[i]);
	}
	for (unsigned i = 0; i < 4; ++i)
	{
		detail->ihv1[i] = SIMD_WTOV(ihv1b[i]);
		detail->ihv2[i] = SIMD_WTOV(ihv2b[i]);
		detail->ihv2mod[i] = SIMD_WTOV(ihv2modb[i]);
		detail->precomp1[i] = SIMD_WTOV(precomp1b[i]);
		detail->precomp2[i] = SIMD_WTOV(precomp2b[i]);
	}
	detail->hybridmask = SIMD_WTOV(hmask);
	detail->distinguishedpointmask = SIMD_WTOV(dpmask);
	detail->maximumpathlength = SIMD_WTOV(maxlen);
	detail->hmask = hmask;
	detail->dpmask = dpmask;
	detail->maxlen = maxlen;
	detail->len.v = SIMD_WTOV(0);
	return true;
}

void birthday_step_avx256(simd_avx256_detail& detail, SIMD_WORD& x, SIMD_WORD& y, SIMD_WORD& z, bool mod);
void simd_device_avx256::fill_trail_buffer(uint64 seed, vector<trail_type>& buf, bool mod)
{
	// generate starting points where necessary
	for (unsigned i = 0; i < SIMD_VECSIZE; ++i)
	{
		if (detail->len.w[i] == 0)
		{
			detail->s0.w[i] = detail->e0.w[i] = uint32(++seed);
			detail->s1.w[i] = detail->e1.w[i] = uint32((seed+=(uint64(1)<<32))>>32);
			detail->s2.w[i] = detail->e2.w[i] = 0;
		}
	}
	for (unsigned k = 0; k < 0x400; ++k)
	{
		birthday_step_avx256(*detail, detail->e0.v, detail->e1.v, detail->e2.v, mod);
		simd_word_t l, t;
		l.v = detail->len.v = SIMD_ADD_VW(detail->len.v,1);
		t.v = SIMD_AND_VV(detail->e0.v, detail->distinguishedpointmask);
		for (unsigned i = 0; i < SIMD_VECSIZE; ++i)
		{
			if (t.w[i] == 0 || l.w[i] > detail->maxlen)
			{
				// if valid point storte then append to buf
				if (t.w[i] == 0)
				{
					buf.emplace_back();
					buf.back().start[0] = detail->s0.w[i];
					buf.back().start[1] = detail->s1.w[i];
					buf.back().start[2] = detail->s2.w[i];
					buf.back().end[0] = detail->e0.w[i];
					buf.back().end[1] = detail->e1.w[i];
					buf.back().end[2] = detail->e2.w[i];
					buf.back().len = detail->len.w[i];
				}
				// generate new starting point
				detail->s0.w[i] = detail->e0.w[i] = uint32(++seed);
				detail->s1.w[i] = detail->e1.w[i] = uint32((seed+=(uint64(1)<<32))>>32);
				detail->s2.w[i] = detail->e2.w[i] = 0;
				detail->len.w[i] = 0;
			}
		}
	}
}

void birthday_step_avx256(simd_avx256_detail& detail, SIMD_WORD& x, SIMD_WORD& y, SIMD_WORD& z, bool mod)
{
	SIMD_WORD block[16];
	SIMD_WORD precomp[4];
	SIMD_WORD ihv[4];
	SIMD_WORD mask = SIMD_EQ_VV(x, SIMD_MIN_VV(x,y));
	for (unsigned i = 0; i < 16; ++i)
		block[i] = SIMD_SEL_VVV(mask, detail.msg1[i], detail.msg2[i]);
	for (unsigned i = 0; i < 4; ++i)
	{
		precomp[i] = SIMD_SEL_VVV(mask, detail.precomp1[i], detail.precomp2[i]);
		ihv[i] = SIMD_SEL_VVV(mask, detail.ihv1[i], detail.ihv2mod[i]);
	}


#define SIMD_MD5_FF(b,c,d) \
		( SIMD_XOR_VV(d, SIMD_AND_VV(b, SIMD_XOR_VV(c, d))) )
//		( SIMD_OR_VV( SIMD_AND_VV(b, c), SIMD_ANDNOT_VV(b, d) ) )
#define SIMD_MD5_GG(b,c,d) \
		( SIMD_XOR_VV(c, SIMD_AND_VV(d, SIMD_XOR_VV(b, c))) )
//		( SIMD_OR_VV( SIMD_AND_VV(d, b), SIMD_ANDNOT_VV(d, c) ) )
#define SIMD_MD5_HH(b,c,d) \
		( SIMD_XOR_VV(b, SIMD_XOR_VV(c, d)) )
#define SIMD_MD5_II(b,c,d) \
		( SIMD_XOR_VV(c, SIMD_OR_VV(b, SIMD_NOT_V(d))) )
#define SIMD_HASHCLASH_MD5COMPRESS_STEP_ff(a, b, c, d, w, ac, rc) \
		a = SIMD_ADD_VV(a, SIMD_ADD_VV( SIMD_ADD_VW(w,ac), SIMD_MD5_FF(b,c,d) ) ); \
		a = SIMD_ROL_V(a, rc); a = SIMD_ADD_VV(a, b);
#define SIMD_HASHCLASH_MD5COMPRESS_STEP_gg(a, b, c, d, w, ac, rc) \
		a = SIMD_ADD_VV(a, SIMD_ADD_VV( SIMD_ADD_VW(w,ac), SIMD_MD5_GG(b,c,d) ) ); \
		a = SIMD_ROL_V(a, rc); a = SIMD_ADD_VV(a, b);
#define SIMD_HASHCLASH_MD5COMPRESS_STEP_hh(a, b, c, d, w, ac, rc) \
		a = SIMD_ADD_VV(a, SIMD_ADD_VV( SIMD_ADD_VW(w,ac), SIMD_MD5_HH(b,c,d) ) ); \
		a = SIMD_ROL_V(a, rc); a = SIMD_ADD_VV(a, b);
#define SIMD_HASHCLASH_MD5COMPRESS_STEP_ii(a, b, c, d, w, ac, rc) \
		a = SIMD_ADD_VV(a, SIMD_ADD_VV( SIMD_ADD_VW(w,ac), SIMD_MD5_II(b,c,d) ) ); \
		a = SIMD_ROL_V(a, rc); a = SIMD_ADD_VV(a, b);

	SIMD_WORD a = precomp[0], b = precomp[1], c = precomp[2], d = precomp[3];
		
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ff( d, a, b, c, z, 0xfd987193, 12);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ff( c, d, a, b, x, 0xa679438e, 17);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ff( b, c, d, a, y, 0x49b40821, 22);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( a, b, c, d, block[ 1], 0xf61e2562,  5);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( d, a, b, c, block[ 6], 0xc040b340,  9);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( c, d, a, b, block[11], 0x265e5a51, 14);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( b, c, d, a, block[ 0], 0xe9b6c7aa, 20); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( a, b, c, d, block[ 5], 0xd62f105d,  5);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( d, a, b, c, block[10], 0x02441453,  9); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( c, d, a, b, y, 0xd8a1e681, 14);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( b, c, d, a, block[ 4], 0xe7d3fbc8, 20); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( a, b, c, d, block[ 9], 0x21e1cde6,  5);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( d, a, b, c, x, 0xc33707d6,  9); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( c, d, a, b, block[ 3], 0xf4d50d87, 14); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( b, c, d, a, block[ 8], 0x455a14ed, 20); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( a, b, c, d, z, 0xa9e3e905,  5); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( d, a, b, c, block[ 2], 0xfcefa3f8,  9);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( c, d, a, b, block[ 7], 0x676f02d9, 14); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_gg( b, c, d, a, block[12], 0x8d2a4c8a, 20);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( a, b, c, d, block[ 5], 0xfffa3942,  4); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( d, a, b, c, block[ 8], 0x8771f681, 11); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( c, d, a, b, block[11], 0x6d9d6122, 16);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( b, c, d, a, x, 0xfde5380c, 23);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( a, b, c, d, block[ 1], 0xa4beea44,  4);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( d, a, b, c, block[ 4], 0x4bdecfa9, 11); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( c, d, a, b, block[ 7], 0xf6bb4b60, 16); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( b, c, d, a, block[10], 0xbebfbc70, 23);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( a, b, c, d, z, 0x289b7ec6,  4); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( d, a, b, c, block[ 0], 0xeaa127fa, 11); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( c, d, a, b, block[ 3], 0xd4ef3085, 16); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( b, c, d, a, block[ 6], 0x04881d05, 23); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( a, b, c, d, block[ 9], 0xd9d4d039,  4);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( d, a, b, c, block[12], 0xe6db99e5, 11);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( c, d, a, b, y, 0x1fa27cf8, 16);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_hh( b, c, d, a, block[ 2], 0xc4ac5665, 23); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( a, b, c, d, block[ 0], 0xf4292244,  6);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( d, a, b, c, block[ 7], 0x432aff97, 10); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( c, d, a, b, x, 0xab9423a7, 15);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( b, c, d, a, block[ 5], 0xfc93a039, 21); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( a, b, c, d, block[12], 0x655b59c3,  6); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( d, a, b, c, block[ 3], 0x8f0ccc92, 10); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( c, d, a, b, block[10], 0xffeff47d, 15);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( b, c, d, a, block[ 1], 0x85845dd1, 21); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( a, b, c, d, block[ 8], 0x6fa87e4f,  6);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( d, a, b, c, y, 0xfe2ce6e0, 10);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( c, d, a, b, block[ 6], 0xa3014314, 15); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( b, c, d, a, z, 0x4e0811a1, 21);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( a, b, c, d, block[ 4], 0xf7537e82,  6);  
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( d, a, b, c, block[11], 0xbd3af235, 10);
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( c, d, a, b, block[ 2], 0x2ad7d2bb, 15); 
		SIMD_HASHCLASH_MD5COMPRESS_STEP_ii( b, c, d, a, block[ 9], 0xeb86d391, 21); 
		
	a = SIMD_ADD_VV(a, ihv[0]);
	b = SIMD_ADD_VV(b, ihv[1]);
	c = SIMD_ADD_VV(c, ihv[2]);
	d = SIMD_ADD_VV(d, ihv[3]);

		if (!mod) {
			// standard multi-block cpc birthday search
			x = a;
			y = SIMD_SUB_VV(d,c);
			z = SIMD_AND_VV(SIMD_SUB_VV(d,b), detail.hybridmask);
		} else {
			// special 1-block cpc birthday search
			x = a;
			y = d;
			z = SIMD_AND_VV(c, detail.hybridmask);
		}

}
#endif
