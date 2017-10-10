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
#include <cutil_inline.h>
#include <cutil.h>
#include <cuda.h>
#include <boost/cstdint.hpp>

using namespace std;

typedef boost::uint32_t uint32;
typedef boost::uint64_t uint64;

#define TRAIL_NOCONSTRUCTOR
#include "birthday_types.hpp"

class cuda_device_detail {
public:
	uint32 device;
	uint32 blocks;
	trail_type* buffer_host;
};

/* We assume that these are _thread specific_ (instead of global) storage managed by the cuda realtime libraries */
__device__ trail_type working_states2[122880]; 
__device__ trail_type buffer2[122880];

__constant__ uint32 msg1[16], msg2[16], ihv1[4], ihv2[4], ihv2mod[4];
__constant__ uint32 precomp1[4], precomp2[4];
__constant__ uint32 hybridmask, distinguishedpointmask, maximumpathlength;


/* F, G and H are basic MD5 functions: selection, majority, parity */
#define MD5_F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define MD5_G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define MD5_H(x, y, z) ((x) ^ (y) ^ (z))
#define MD5_I(x, y, z) ((y) ^ ((x) | (~z))) 

/* ROTATE_LEFT rotates x left n bits */
#define ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32-(n))))

/* FF, GG, HH, and II transformations for rounds 1, 2, 3, and 4 */
/* Rotation is separate from addition to prevent recomputation */
#define MD5_FF(a, b, c, d, x, s, ac) \
  {(a) += MD5_F ((b), (c), (d)) + (x) + (uint32)(ac); \
   (a) = ROTATE_LEFT ((a), (s)); \
   (a) += (b); \
  }
#define MD5_GG(a, b, c, d, x, s, ac) \
  {(a) += MD5_G ((b), (c), (d)) + (x) + (uint32)(ac); \
   (a) = ROTATE_LEFT ((a), (s)); \
   (a) += (b); \
  }
#define MD5_HH(a, b, c, d, x, s, ac) \
  {(a) += MD5_H ((b), (c), (d)) + (x) + (uint32)(ac); \
   (a) = ROTATE_LEFT ((a), (s)); \
   (a) += (b); \
  }
#define MD5_II(a, b, c, d, x, s, ac) \
  {(a) += MD5_I ((b), (c), (d)) + (x) + (uint32)(ac); \
   (a) = ROTATE_LEFT ((a), (s)); \
   (a) += (b); \
  }


__global__ void cuda_md5_init()
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	working_states2[idx].len = 0;
	buffer2[idx].len = 0;
}

bool cuda_device::init(uint32 device, const uint32 ihv1b[4], const uint32 ihv2b[4], const uint32 ihv2modb[4], const uint32 msg1b[16], const uint32 msg2b[16], uint32 hmask, uint32 dpmask, uint32 maxlen)
{
	detail = new cuda_device_detail;
	detail->device = device;
	
    int deviceCount;
    CUDA_SAFE_CALL( cudaGetDeviceCount(&deviceCount) );
    if (deviceCount == 0) {
        cout << "There is no device supporting CUDA!" << endl;
        return false;
	}
    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL( cudaGetDeviceProperties(&deviceProp, device) );
    if (deviceProp.major == 9999) {
		cout << "Emulation device found." << endl;
		return false;
	}
	cout << "CUDA device " << device << ": " << deviceProp.name << " (" << 8 * deviceProp.multiProcessorCount << " cores)" << endl;
	detail->blocks = 16 * deviceProp.multiProcessorCount;
	
	CUDA_SAFE_CALL( cudaSetDevice(device) );
	CUDA_SAFE_CALL( cudaSetDeviceFlags( cudaDeviceBlockingSync ) );

	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->buffer_host)), 122880 * sizeof(trail_type) ) );
	
	uint32 pc1[4], pc2[4];
	uint32 a = ihv1b[0], b = ihv1b[1], c = ihv1b[2], d = ihv1b[3];
	MD5_FF ( a, b, c, d, msg1b[ 0],  7, 3614090360); /* 1 */
	MD5_FF ( d, a, b, c, msg1b[ 1], 12, 3905402710); /* 2 */
	MD5_FF ( c, d, a, b, msg1b[ 2], 17,  606105819); /* 3 */
	MD5_FF ( b, c, d, a, msg1b[ 3], 22, 3250441966); /* 4 */
	MD5_FF ( a, b, c, d, msg1b[ 4],  7, 4118548399); /* 5 */
	MD5_FF ( d, a, b, c, msg1b[ 5], 12, 1200080426); /* 6 */
	MD5_FF ( c, d, a, b, msg1b[ 6], 17, 2821735955); /* 7 */
	MD5_FF ( b, c, d, a, msg1b[ 7], 22, 4249261313); /* 8 */
	MD5_FF ( a, b, c, d, msg1b[ 8],  7, 1770035416); /* 9 */
	MD5_FF ( d, a, b, c, msg1b[ 9], 12, 2336552879); /* 10 */
	MD5_FF ( c, d, a, b, msg1b[10], 17, 4294925233); /* 11 */
	MD5_FF ( b, c, d, a, msg1b[11], 22, 2304563134); /* 12 */
	MD5_FF ( a, b, c, d, msg1b[12],  7, 1804603682); /* 13 */
	pc1[0] = a; pc1[1] = b; pc1[2] = c; pc1[3] = d;
	a = ihv2b[0]; b = ihv2b[1]; c = ihv2b[2]; d = ihv2b[3];
	MD5_FF ( a, b, c, d, msg2b[ 0],  7, 3614090360); /* 1 */
	MD5_FF ( d, a, b, c, msg2b[ 1], 12, 3905402710); /* 2 */
	MD5_FF ( c, d, a, b, msg2b[ 2], 17,  606105819); /* 3 */
	MD5_FF ( b, c, d, a, msg2b[ 3], 22, 3250441966); /* 4 */
	MD5_FF ( a, b, c, d, msg2b[ 4],  7, 4118548399); /* 5 */
	MD5_FF ( d, a, b, c, msg2b[ 5], 12, 1200080426); /* 6 */
	MD5_FF ( c, d, a, b, msg2b[ 6], 17, 2821735955); /* 7 */
	MD5_FF ( b, c, d, a, msg2b[ 7], 22, 4249261313); /* 8 */
	MD5_FF ( a, b, c, d, msg2b[ 8],  7, 1770035416); /* 9 */
	MD5_FF ( d, a, b, c, msg2b[ 9], 12, 2336552879); /* 10 */
	MD5_FF ( c, d, a, b, msg2b[10], 17, 4294925233); /* 11 */
	MD5_FF ( b, c, d, a, msg2b[11], 22, 2304563134); /* 12 */
	MD5_FF ( a, b, c, d, msg2b[12],  7, 1804603682); /* 13 */
	pc2[0] = a; pc2[1] = b; pc2[2] = c; pc2[3] = d;

	CUDA_SAFE_CALL( cudaMemcpyToSymbol(msg1, msg1b, sizeof(msg1)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(msg2, msg2b, sizeof(msg2)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(ihv1, ihv1b, sizeof(ihv1)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(ihv2, ihv2b, sizeof(ihv2)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(ihv2mod, ihv2modb, sizeof(ihv2mod)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(precomp1, pc1, sizeof(pc1)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(precomp2, pc2, sizeof(pc2)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(hybridmask, &hmask, sizeof(hmask)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(distinguishedpointmask, &dpmask, sizeof(dpmask)) );
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(maximumpathlength, &maxlen, sizeof(maxlen)) );
	
	cuda_md5_init<<<detail->blocks, 256>>>();
	
	return true;
}

__global__ void cuda_md5_work(uint64 seed)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	buffer2[idx].len = 0;
	uint32 len = working_states2[idx].len;	
	uint32 x = working_states2[idx].end[0];
	uint32 y = working_states2[idx].end[1];
	uint32 z = working_states2[idx].end[2];
	if (len >= maximumpathlength || len == 0) {
		x = uint32(seed>>32) ^ threadIdx.x;
		y = uint32(seed) ^ blockIdx.x;
		z = 0;
		working_states2[idx].start[0] = x;
		working_states2[idx].start[1] = y;
		working_states2[idx].start[2] = z;
		len = 0;
	}	
	__syncthreads();	

	for (unsigned j = 0; j < 0x100; ++j)
	{
		{
			uint32* in = msg1;
			uint32 a = precomp1[0], b = precomp1[1], c = precomp1[2], d = precomp1[3];
			if (x > y) {
				in = msg2;
				a = precomp2[0]; b = precomp2[1]; c = precomp2[2]; d = precomp2[3];
			}
			MD5_FF ( d, a, b, c, z, 12, 4254626195); /* 14 */
			MD5_FF ( c, d, a, b, x, 17, 2792965006); /* 15 */
			MD5_FF ( b, c, d, a, y, 22, 1236535329); /* 16 */

			MD5_GG ( a, b, c, d, in[ 1],  5, 4129170786); /* 17 */
			MD5_GG ( d, a, b, c, in[ 6],  9, 3225465664); /* 18 */
			MD5_GG ( c, d, a, b, in[11], 14,  643717713); /* 19 */
			MD5_GG ( b, c, d, a, in[ 0], 20, 3921069994); /* 20 */
			MD5_GG ( a, b, c, d, in[ 5],  5, 3593408605); /* 21 */
			MD5_GG ( d, a, b, c, in[10],  9,   38016083); /* 22 */
			MD5_GG ( c, d, a, b, y, 14, 3634488961); /* 23 */
			MD5_GG ( b, c, d, a, in[ 4], 20, 3889429448); /* 24 */
			MD5_GG ( a, b, c, d, in[ 9],  5,  568446438); /* 25 */
			MD5_GG ( d, a, b, c, x,  9, 3275163606); /* 26 */
			MD5_GG ( c, d, a, b, in[ 3], 14, 4107603335); /* 27 */
			MD5_GG ( b, c, d, a, in[ 8], 20, 1163531501); /* 28 */
			MD5_GG ( a, b, c, d, z,  5, 2850285829); /* 29 */
			MD5_GG ( d, a, b, c, in[ 2],  9, 4243563512); /* 30 */
			MD5_GG ( c, d, a, b, in[ 7], 14, 1735328473); /* 31 */
			MD5_GG ( b, c, d, a, in[12], 20, 2368359562); /* 32 */

			MD5_HH ( a, b, c, d, in[ 5],  4, 4294588738); /* 33 */
			MD5_HH ( d, a, b, c, in[ 8], 11, 2272392833); /* 34 */
			MD5_HH ( c, d, a, b, in[11], 16, 1839030562); /* 35 */
			MD5_HH ( b, c, d, a, x, 23, 4259657740); /* 36 */
			MD5_HH ( a, b, c, d, in[ 1],  4, 2763975236); /* 37 */
			MD5_HH ( d, a, b, c, in[ 4], 11, 1272893353); /* 38 */
			MD5_HH ( c, d, a, b, in[ 7], 16, 4139469664); /* 39 */
			MD5_HH ( b, c, d, a, in[10], 23, 3200236656); /* 40 */
			MD5_HH ( a, b, c, d, z,  4,  681279174); /* 41 */
			MD5_HH ( d, a, b, c, in[ 0], 11, 3936430074); /* 42 */
			MD5_HH ( c, d, a, b, in[ 3], 16, 3572445317); /* 43 */
			MD5_HH ( b, c, d, a, in[ 6], 23,   76029189); /* 44 */
			MD5_HH ( a, b, c, d, in[ 9],  4, 3654602809); /* 45 */
			MD5_HH ( d, a, b, c, in[12], 11, 3873151461); /* 46 */
			MD5_HH ( c, d, a, b, y, 16,  530742520); /* 47 */
			MD5_HH ( b, c, d, a, in[ 2], 23, 3299628645); /* 48 */

			MD5_II ( a, b, c, d, in[ 0],  6, 4096336452); /* 49 */
			MD5_II ( d, a, b, c, in[ 7], 10, 1126891415); /* 50 */
			MD5_II ( c, d, a, b, x, 15, 2878612391); /* 51 */
			MD5_II ( b, c, d, a, in[ 5], 21, 4237533241); /* 52 */
			MD5_II ( a, b, c, d, in[12],  6, 1700485571); /* 53 */
			MD5_II ( d, a, b, c, in[ 3], 10, 2399980690); /* 54 */
			MD5_II ( c, d, a, b, in[10], 15, 4293915773); /* 55 */
			MD5_II ( b, c, d, a, in[ 1], 21, 2240044497); /* 56 */
			MD5_II ( a, b, c, d, in[ 8],  6, 1873313359); /* 57 */
			MD5_II ( d, a, b, c, y, 10, 4264355552); /* 58 */
			MD5_II ( c, d, a, b, in[ 6], 15, 2734768916); /* 59 */
			MD5_II ( b, c, d, a, z, 21, 1309151649); /* 60 */
			MD5_II ( a, b, c, d, in[ 4],  6, 4149444226); /* 61 */
			MD5_II ( d, a, b, c, in[11], 10, 3174756917); /* 62 */
			MD5_II ( c, d, a, b, in[ 2], 15,  718787259); /* 63 */
			MD5_II ( b, c, d, a, in[ 9], 21, 3951481745); /* 64 */

			if (x <= y) {
				a += ihv1[0];
				b += ihv1[1];
				c += ihv1[2];
				d += ihv1[3];
			} else {
				a += ihv2mod[0];
				b += ihv2mod[1];
				c += ihv2mod[2];
				d += ihv2mod[3];
			}
			x = a;
			y = d - c;
			z = (d - b) & hybridmask;
			++len;
		}
		
		{
			if (0 == (x & distinguishedpointmask)) {
				buffer2[idx].end[0] = x;
				buffer2[idx].end[1] = y;
				buffer2[idx].end[2] = z;
				buffer2[idx].len = len;
				buffer2[idx].start[0] = working_states2[idx].start[0];
				buffer2[idx].start[1] = working_states2[idx].start[1];
				buffer2[idx].start[2] = working_states2[idx].start[2];
				x = uint32(seed>>32) ^ (threadIdx.x<<16) + len;
				y = uint32(seed) ^ blockIdx.x;
				z = 0;
				len = 0;
				working_states2[idx].start[0] = x;
				working_states2[idx].start[1] = y;
				working_states2[idx].start[2] = z;
			}
		}
		__syncthreads();
	}

	working_states2[idx].end[0] = x;
	working_states2[idx].end[1] = y;
	working_states2[idx].end[2] = z;
	working_states2[idx].len = len;
}


__global__ void cuda_md5_workmod(uint64 seed)
{
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	buffer2[idx].len = 0;
	uint32 len = working_states2[idx].len;	
	uint32 x = working_states2[idx].end[0];
	uint32 y = working_states2[idx].end[1];
	uint32 z = working_states2[idx].end[2];
	if (len >= maximumpathlength || len == 0) {
		x = uint32(seed>>32) ^ threadIdx.x;
		y = uint32(seed) ^ blockIdx.x;
		z = 0;
		working_states2[idx].start[0] = x;
		working_states2[idx].start[1] = y;
		working_states2[idx].start[2] = z;
		len = 0;
	}	
	__syncthreads();	

	for (unsigned j = 0; j < 0x100; ++j)
	{
		{
			uint32* in = msg1;
			uint32 a = precomp1[0], b = precomp1[1], c = precomp1[2], d = precomp1[3];
			if (x > y) {
				in = msg2;
				a = precomp2[0]; b = precomp2[1]; c = precomp2[2]; d = precomp2[3];
			}
			MD5_FF ( d, a, b, c, z, 12, 4254626195); /* 14 */
			MD5_FF ( c, d, a, b, x, 17, 2792965006); /* 15 */
			MD5_FF ( b, c, d, a, y, 22, 1236535329); /* 16 */

			MD5_GG ( a, b, c, d, in[ 1],  5, 4129170786); /* 17 */
			MD5_GG ( d, a, b, c, in[ 6],  9, 3225465664); /* 18 */
			MD5_GG ( c, d, a, b, in[11], 14,  643717713); /* 19 */
			MD5_GG ( b, c, d, a, in[ 0], 20, 3921069994); /* 20 */
			MD5_GG ( a, b, c, d, in[ 5],  5, 3593408605); /* 21 */
			MD5_GG ( d, a, b, c, in[10],  9,   38016083); /* 22 */
			MD5_GG ( c, d, a, b, y, 14, 3634488961); /* 23 */
			MD5_GG ( b, c, d, a, in[ 4], 20, 3889429448); /* 24 */
			MD5_GG ( a, b, c, d, in[ 9],  5,  568446438); /* 25 */
			MD5_GG ( d, a, b, c, x,  9, 3275163606); /* 26 */
			MD5_GG ( c, d, a, b, in[ 3], 14, 4107603335); /* 27 */
			MD5_GG ( b, c, d, a, in[ 8], 20, 1163531501); /* 28 */
			MD5_GG ( a, b, c, d, z,  5, 2850285829); /* 29 */
			MD5_GG ( d, a, b, c, in[ 2],  9, 4243563512); /* 30 */
			MD5_GG ( c, d, a, b, in[ 7], 14, 1735328473); /* 31 */
			MD5_GG ( b, c, d, a, in[12], 20, 2368359562); /* 32 */

			MD5_HH ( a, b, c, d, in[ 5],  4, 4294588738); /* 33 */
			MD5_HH ( d, a, b, c, in[ 8], 11, 2272392833); /* 34 */
			MD5_HH ( c, d, a, b, in[11], 16, 1839030562); /* 35 */
			MD5_HH ( b, c, d, a, x, 23, 4259657740); /* 36 */
			MD5_HH ( a, b, c, d, in[ 1],  4, 2763975236); /* 37 */
			MD5_HH ( d, a, b, c, in[ 4], 11, 1272893353); /* 38 */
			MD5_HH ( c, d, a, b, in[ 7], 16, 4139469664); /* 39 */
			MD5_HH ( b, c, d, a, in[10], 23, 3200236656); /* 40 */
			MD5_HH ( a, b, c, d, z,  4,  681279174); /* 41 */
			MD5_HH ( d, a, b, c, in[ 0], 11, 3936430074); /* 42 */
			MD5_HH ( c, d, a, b, in[ 3], 16, 3572445317); /* 43 */
			MD5_HH ( b, c, d, a, in[ 6], 23,   76029189); /* 44 */
			MD5_HH ( a, b, c, d, in[ 9],  4, 3654602809); /* 45 */
			MD5_HH ( d, a, b, c, in[12], 11, 3873151461); /* 46 */
			MD5_HH ( c, d, a, b, y, 16,  530742520); /* 47 */
			MD5_HH ( b, c, d, a, in[ 2], 23, 3299628645); /* 48 */

			MD5_II ( a, b, c, d, in[ 0],  6, 4096336452); /* 49 */
			MD5_II ( d, a, b, c, in[ 7], 10, 1126891415); /* 50 */
			MD5_II ( c, d, a, b, x, 15, 2878612391); /* 51 */
			MD5_II ( b, c, d, a, in[ 5], 21, 4237533241); /* 52 */
			MD5_II ( a, b, c, d, in[12],  6, 1700485571); /* 53 */
			MD5_II ( d, a, b, c, in[ 3], 10, 2399980690); /* 54 */
			MD5_II ( c, d, a, b, in[10], 15, 4293915773); /* 55 */
			MD5_II ( b, c, d, a, in[ 1], 21, 2240044497); /* 56 */
			MD5_II ( a, b, c, d, in[ 8],  6, 1873313359); /* 57 */
			MD5_II ( d, a, b, c, y, 10, 4264355552); /* 58 */
			MD5_II ( c, d, a, b, in[ 6], 15, 2734768916); /* 59 */
			MD5_II ( b, c, d, a, z, 21, 1309151649); /* 60 */
			MD5_II ( a, b, c, d, in[ 4],  6, 4149444226); /* 61 */
			MD5_II ( d, a, b, c, in[11], 10, 3174756917); /* 62 */
			MD5_II ( c, d, a, b, in[ 2], 15,  718787259); /* 63 */			

			if (x <= y) {
				x = a + ihv1[0];
				y = d + ihv1[3];
				z = (c + ihv1[2]) & hybridmask;
			} else {
				x = a + ihv2mod[0];
				y = d + ihv2mod[3];
				z = (c + ihv2mod[2]) & hybridmask;
			}
			++len;
		}
		
		{
			if (0 == (x & distinguishedpointmask)) {
				buffer2[idx].end[0] = x;
				buffer2[idx].end[1] = y;
				buffer2[idx].end[2] = z;
				buffer2[idx].len = len;
				buffer2[idx].start[0] = working_states2[idx].start[0];
				buffer2[idx].start[1] = working_states2[idx].start[1];
				buffer2[idx].start[2] = working_states2[idx].start[2];
				x = uint32(seed>>32) ^ (threadIdx.x<<16) + len;
				y = uint32(seed) ^ blockIdx.x;
				z = 0;
				len = 0;
				working_states2[idx].start[0] = x;
				working_states2[idx].start[1] = y;
				working_states2[idx].start[2] = z;
			}
		}
		__syncthreads();
	}

	working_states2[idx].end[0] = x;
	working_states2[idx].end[1] = y;
	working_states2[idx].end[2] = z;
	working_states2[idx].len = len;
}

void cuda_device::cuda_fill_trail_buffer(uint32 id, uint64 seed,
							vector<trail_type>& buf,
							vector< pair<trail_type,trail_type> >& collisions, bool mod)
{	
	// transfer results
	cudaMemcpyFromSymbol(detail->buffer_host, buffer2, sizeof(trail_type)*detail->blocks*256);
	
	// start new cuda computation
	if (!mod)
		cuda_md5_work<<<detail->blocks, 256>>>(seed);
	else
		cuda_md5_workmod<<<detail->blocks, 256>>>(seed);
	
	// process and return results
	buf.clear();
	for (unsigned i = 0; i < detail->blocks*256; ++i)
		if (detail->buffer_host[i].len)
			buf.push_back(detail->buffer_host[i]);
}							










#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

	class timer_detail;
	class timer {
	public:
		timer(bool direct_start = false);
		~timer();
		void start();
		void stop();
		double time() const;// get time between start and stop (or now if still running) in seconds
		bool isrunning() const { return running; } // check if timer is running

	private:
		timer_detail* detail;
		bool running;
	};
	class timer_detail {
	public:
#ifdef _WIN32
		LARGE_INTEGER tstart, tend;
		double freq;
#else
		struct timeval tstart, tend;
		struct timezone tz;
#endif
	};

	timer::~timer()
	{
		delete detail;
	}

	timer::timer(bool direct_start): running(false) 
	{ 
		detail = new timer_detail;
#ifdef _WIN32
		LARGE_INTEGER tmp_freq;
		QueryPerformanceFrequency(&tmp_freq);
		detail->freq = double(tmp_freq.QuadPart);
#endif
		if (direct_start)
			start();
	}

#ifdef _WIN32

	void timer::start()
	{
		running = true;
		QueryPerformanceCounter(&detail->tstart);
	}

	void timer::stop()
	{
		QueryPerformanceCounter(&detail->tend);
		running = false;
	}

	double timer::time() const
	{
		if (running)
		{
			LARGE_INTEGER tmp_end;
			QueryPerformanceCounter(&tmp_end);
			return (double(tmp_end.QuadPart) - double(detail->tstart.QuadPart))/detail->freq;
		} else 
			return (double(detail->tend.QuadPart) - double(detail->tstart.QuadPart))/detail->freq;
	}

#else

	void timer::start()
	{
		running = true;
		gettimeofday(&detail->tstart, &detail->tz);
	}

	void timer::stop()
	{
		gettimeofday(&detail->tend, &detail->tz);
		running = false;
	}

	double timer::time() const
	{
		double t1 = double(detail->tstart.tv_sec) + (double(detail->tstart.tv_usec)/1e6);
		if (running)
		{
			struct timeval tmp_end;
			gettimeofday(&tmp_end, &detail->tz);
			return double(tmp_end.tv_sec) + (double(tmp_end.tv_usec)/1e6) - t1;
		} else
			return double(detail->tend.tv_sec) + (double(detail->tend.tv_usec)/1e6) - t1;
	}

#endif

void cuda_device::benchmark()
{
	timer sw;
	for (int blocksize = 4; blocksize <= 256; ++blocksize)
	for (int threadsize = 250; threadsize <= 257; ++threadsize)
	{
		sw.start();
		uint64 work = 0;
		while (sw.time() < 10) {
			cuda_md5_work<<<blocksize, threadsize>>>(0);
			cudaMemcpyFromSymbol(detail->buffer_host, buffer2, sizeof(trail_type)*blocksize*threadsize);
			++work;
		}
		uint64 ow = work;
		work *= 0x400 * blocksize * threadsize;
		cout << blocksize << "x" << threadsize << ":\t" << work << " (" << ow << ")" << endl;
	}
}














int get_num_cuda_devices()
{
	int deviceCount = 0;
	cudaGetDeviceCount(&deviceCount);
	return deviceCount;
}

void cuda_device_query() 
{
    int deviceCount;
    cutilSafeCall(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    int dev;
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cutilSafeCall(cudaGetDeviceProperties(&deviceProp, dev));
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
        }
        printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
        printf("  Major revision number:                         %d\n",
               deviceProp.major);
        printf("  Minor revision number:                         %d\n",
               deviceProp.minor);
        printf("  Total amount of global memory:                 %u bytes\n",
               deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
        printf("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
        printf("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
    #endif
        printf("  Total amount of constant memory:               %u bytes\n",
               deviceProp.totalConstMem); 
        printf("  Total amount of shared memory per block:       %u bytes\n",
               deviceProp.sharedMemPerBlock);
        printf("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
        printf("  Warp size:                                     %d\n",
               deviceProp.warpSize);
        printf("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf("  Maximum memory pitch:                          %u bytes\n",
               deviceProp.memPitch);
        printf("  Texture alignment:                             %u bytes\n",
               deviceProp.textureAlignment);
        printf("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        printf("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    }

}
