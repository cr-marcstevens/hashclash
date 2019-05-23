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
#include <cuda.h>
#include <stdexcept>
#include <boost/cstdint.hpp>
#include "cuda_cyclicbuffer.hpp"

using namespace std;

typedef boost::uint32_t uint32;
typedef boost::uint64_t uint64;

#define MAX_CUDA_THREADS (1<<20)
#define MAX_CUDA_BLOCKS 256
#define MAX_CUDA_THREADS_PER_BLOCK 2048
#define REGISTERS_PER_CUDA_THREAD 64

#define TRAIL_NOCONSTRUCTOR
#include "birthday_types.hpp"

#ifndef CUDA_SAFE_CALL
#define CUDA_SAFE_CALL(s) { auto ce = s; if (ce != cudaSuccess) { throw std::runtime_error("CUDA API Error:\n" + std::string(cudaGetErrorName(ce)) + ":\n" + std::string(cudaGetErrorString(ce))); } }
#endif

#ifndef cutilSafeCall
#define cutilSafeCall(s) (s)
#endif


/****
  NOTE WARNING: 
  We assume that all global __device__ variables below are *thread* *specific*
 (instead of global) storage managed by the cuda realtime libraries 
*****/
// last template parameter is fence type: 0=none, 1=block, 2=gpu
typedef cyclic_buffer_cas_t<MAX_CUDA_THREADS,uint32,7,cyclic_buffer_control_cas_t<MAX_CUDA_THREADS>,2> state_buffer_t;
typedef cyclic_buffer_mask_t<MAX_CUDA_THREADS_PER_BLOCK,uint32,7,cyclic_buffer_control_mask_t<MAX_CUDA_THREADS_PER_BLOCK>,1> work_buffer_t;
typedef work_buffer_t::control_t work_control_t;
typedef cyclic_buffer_cas_t<MAX_CUDA_THREADS,uint32,15,cyclic_buffer_control_cas_t<MAX_CUDA_THREADS>,2> collisions_buffer_t;
typedef collisions_buffer_t::control_t collisions_control_t;

// static gpu buffer that always stays on GPU
__device__ state_buffer_t gworking_states;
__device__ collisions_buffer_t gcollision_states;

// per-block buffer for trails
__device__ work_buffer_t gtrailsout_buf[MAX_CUDA_BLOCKS];
__device__ work_control_t gtrailsout_ctl[MAX_CUDA_BLOCKS];
__shared__ work_control_t gtrailsout_ctlblock;

// gpu-wide in- and out-put buffers collisions
__device__ collisions_buffer_t gcollisionsin_buf;
__device__ collisions_buffer_t gcollisionsout_buf;
__device__ collisions_control_t gcollisionsin_ctl;
__device__ collisions_control_t gcollisionsout_ctl;
__device__ volatile uint32 halt_flag;

__constant__ uint32 msg1[16], msg2[16], ihv1[4], ihv2[4], ihv2mod[4];
__constant__ uint32 precomp1[4], precomp2[4];
__constant__ uint32 hybridmask, distinguishedpointmask, maximumpathlength;


class cuda_device_detail {
public:
	uint32 device;
	uint32 blocks;
	uint32 threadsperblock;
	work_buffer_t* trailsout_buf;
	work_control_t* trailsout_ctl;

	// host-side buffer
	size_t nrcollisions_on_gpu;
	vector< pair<trail_type,trail_type> > collisions;
	collisions_buffer_t* collisionsin_buf;
	collisions_control_t* collisionsin_ctl;
	collisions_buffer_t* collisionsout_buf;
	collisions_control_t* collisionsout_ctl;
};


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


__device__ void backup_controls()
{
	__syncthreads();
	if (threadIdx.x == 0)
	{
		gtrailsout_ctlblock = gtrailsout_ctl[blockIdx.x];
	}
	__syncthreads();
}

__device__ void restore_controls()
{
	__syncthreads();
	if (threadIdx.x == 0)
	{
		 gtrailsout_ctl[blockIdx.x] = gtrailsout_ctlblock;
	}
	__syncthreads();
}


__global__ void cuda_md5_init()
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	gworking_states.get_ref<6>(idx) = 0; // len = 0
	gcollision_states.get_ref<14>(idx) = 1; // bad = 1
	if (threadIdx.x == 0)
	{
		gtrailsout_buf[blockIdx.x].reset(gtrailsout_ctl[blockIdx.x]);
		gcollisionsin_buf.reset(gcollisionsin_ctl);
		gcollisionsout_buf.reset(gcollisionsout_ctl);
	}

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
	cout << "CUDA device " << device << ": " << deviceProp.name << " (" << deviceProp.multiProcessorCount << " MPs)" << endl;
	unsigned maxthreadspermp = deviceProp.maxThreadsPerMultiProcessor;
	if (maxthreadspermp > MAX_CUDA_THREADS)
		maxthreadspermp = (MAX_CUDA_THREADS/32)*32;
	while (maxthreadspermp > deviceProp.regsPerMultiprocessor * REGISTERS_PER_CUDA_THREAD)
		maxthreadspermp -= 32;
	unsigned minblockspermp = 1;
	while (maxthreadspermp > minblockspermp * deviceProp.maxThreadsPerBlock)
		minblockspermp += 1;
	while (maxthreadspermp * REGISTERS_PER_CUDA_THREAD > minblockspermp * deviceProp.regsPerBlock)
		minblockspermp += 1;

	detail->threadsperblock = ((maxthreadspermp / minblockspermp) / 32) * 32;
	detail->blocks = minblockspermp * deviceProp.multiProcessorCount;
	cout << "Using " << detail->blocks << " blocks with " << detail->threadsperblock << " threads each: total " << detail->blocks * detail->threadsperblock << " threads." << endl;

	CUDA_SAFE_CALL( cudaSetDevice(device) );
//	CUDA_SAFE_CALL( cudaSetDeviceFlags( cudaDeviceBlockingSync ) );

//	work_buffer_t* trailsout_buf;//[MAX_CUDA_BLOCKS];
//	work_control_t* trailsout_ctl;//[MAX_CUDA_BLOCKS];
//	collisions_buffer_t* collisionsin_buf;//[MAX_CUDA_BLOCKS];
//	collisions_control_t* collisionsin_ctl;//[MAX_CUDA_BLOCKS];
//	collisions_buffer_t* collisionsout_buf;//[MAX_CUDA_BLOCKS];
//	collisions_control_t* collisionsout_ctl;//[MAX_CUDA_BLOCKS];

	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->trailsout_buf)), detail->blocks * sizeof(work_buffer_t) ) );
	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->trailsout_ctl)), detail->blocks * sizeof(work_control_t) ) );
	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->collisionsin_buf)),  sizeof(collisions_buffer_t) ) );
	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->collisionsin_ctl)),  sizeof(collisions_control_t) ) );
	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->collisionsout_buf)), sizeof(collisions_buffer_t) ) );
	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->collisionsout_ctl)), sizeof(collisions_control_t) ) );

	for (unsigned b = 0; b < detail->blocks; ++b)
		detail->trailsout_buf[b].reset(detail->trailsout_ctl[b]);
	detail->collisionsin_buf->reset(*(detail->collisionsin_ctl));
	detail->collisionsout_buf->reset(*(detail->collisionsout_ctl));

	detail->nrcollisions_on_gpu = 0;

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


	cuda_md5_init<<<detail->blocks, detail->threadsperblock>>>();

	return true;
}




template<bool mod = false>
__device__ void cuda_md5_work2(uint64 seed)
{
	halt_flag = 0;
	/********************* GENERATE TRAILS ***********************/
	restore_controls();
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	uint32 len = gworking_states.get<6>(idx);
	uint32 x = gworking_states.get<3>(idx); //end[0]
	uint32 y = gworking_states.get<4>(idx); //end[1]
	uint32 z = gworking_states.get<5>(idx); //end[2]
	if (len >= maximumpathlength || len == 0) {
		x = uint32(seed>>32) ^ threadIdx.x;
		y = uint32(seed) ^ blockIdx.x;
		z = 0;
		gworking_states.get_ref<0>(idx) = x;
		gworking_states.get_ref<1>(idx) = y;
		gworking_states.get_ref<2>(idx) = z;
		len = 0;
	}

	for (unsigned j = 0; j < (1<<12); ++j)
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

			if (mod)
			{
				if (x <= y) {
					x = a + ihv1[0];
					y = d + ihv1[3];
					z = (c + ihv1[2]) & hybridmask;
				} else {
					x = a + ihv2mod[0];
					y = d + ihv2mod[3];
					z = (c + ihv2mod[2]) & hybridmask;
				}
			}
			else
			{
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
			}
			++len;
		}

		{
			// conditionally write
			bool done = (0 == (x & distinguishedpointmask));
			gtrailsout_buf[blockIdx.x].write(gtrailsout_ctlblock, done,
				gworking_states.get_ref<0>(idx),
				gworking_states.get_ref<1>(idx),
				gworking_states.get_ref<2>(idx),
				x, y, z, len
				);
			if (done)
			{
				x = uint32(seed>>32) ^ (threadIdx.x<<16) + len;
				y = uint32(seed) ^ blockIdx.x;
				z = 0;
				len = 0;
				gworking_states.get_ref<0>(idx) = x;
				gworking_states.get_ref<1>(idx) = y;
				gworking_states.get_ref<2>(idx) = z;
			}
		}
//		__syncthreads();
	}

	gworking_states.get_ref<3>(idx) = x;
	gworking_states.get_ref<4>(idx) = y;
	gworking_states.get_ref<5>(idx) = z;
	gworking_states.get_ref<6>(idx) = len;
	backup_controls();
	halt_flag = 1;
}

template<bool mod = false>
__global__ void cuda_md5_work(uint64 seed)
{
	cuda_md5_work2<mod>(seed);
}


template<bool mod = false>
__global__ void cuda_md5_collisions(uint64 seed)
{
	halt_flag = 0;
	/********** PROCESS COLLIDING TRAILS INTO COLLISIONS ***************/
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	uint32 bad = gcollision_states.get<14>(idx);
	uint32 len = gcollision_states.get<6>(idx);
	uint32 len2 = gcollision_states.get<7+6>(idx);
	// if collision state is empty then go read a collision
	if (len == 0 || len2 == 0)
		bad = 1;
	uint32 readidx = gcollisionsin_buf.getreadidx(gcollisionsin_ctl,bad);
	if (bad && readidx < 0xEEEEEEEE)
	{
		len  = gcollisionsin_buf.get<6>(readidx);
		len2 = gcollisionsin_buf.get<7+6>(readidx);
		gcollision_states.get_ref<0>(idx) = gcollisionsin_buf.get<0>(readidx); //start[0]
		gcollision_states.get_ref<1>(idx) = gcollisionsin_buf.get<1>(readidx); //start[1]
		gcollision_states.get_ref<2>(idx) = gcollisionsin_buf.get<2>(readidx); //start[2]
		gcollision_states.get_ref<3>(idx) = gcollisionsin_buf.get<0>(readidx); //start[0]
		gcollision_states.get_ref<4>(idx) = gcollisionsin_buf.get<1>(readidx); //start[1]
		gcollision_states.get_ref<5>(idx) = gcollisionsin_buf.get<2>(readidx); //start[2]
		gcollision_states.get_ref<6>(idx) = gcollisionsin_buf.get<6>(readidx); //len
		gcollision_states.get_ref<7+0>(idx) = gcollisionsin_buf.get<7+0>(readidx);
		gcollision_states.get_ref<7+1>(idx) = gcollisionsin_buf.get<7+1>(readidx);
		gcollision_states.get_ref<7+2>(idx) = gcollisionsin_buf.get<7+2>(readidx);
		gcollision_states.get_ref<7+3>(idx) = gcollisionsin_buf.get<7+0>(readidx);
		gcollision_states.get_ref<7+4>(idx) = gcollisionsin_buf.get<7+1>(readidx);
		gcollision_states.get_ref<7+5>(idx) = gcollisionsin_buf.get<7+2>(readidx);
		gcollision_states.get_ref<7+6>(idx) = gcollisionsin_buf.get<7+6>(readidx);
		gcollision_states.get_ref<14>(idx) = bad = 0;
	}
	if (__all_sync(WARP_FULL_MASK,bad))
	{
//		cuda_md5_work2<mod>(seed);
		return;
	}
	for (unsigned j = 0; j < (1<<12); ++j)
	{
		// always process the longest
		uint32 x, y, z;
		if (len >= len2)
		{
			// process trail1
			// load start+1, write to start
			x = gcollision_states.get<3>(idx);
			y = gcollision_states.get<4>(idx);
			z = gcollision_states.get<5>(idx);
			gcollision_states.get_ref<0>(idx) = x;
			gcollision_states.get_ref<1>(idx) = y;
			gcollision_states.get_ref<2>(idx) = z;
		}
		else
		{
			// process trail2
			// load start+1, write to start
			x = gcollision_states.get<7+3>(idx);
			y = gcollision_states.get<7+4>(idx);
			z = gcollision_states.get<7+5>(idx);
			gcollision_states.get_ref<7+0>(idx) = x;
			gcollision_states.get_ref<7+1>(idx) = y;
			gcollision_states.get_ref<7+2>(idx) = z;
		}

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

			if (mod)
			{
				if (x <= y) {
					x = a + ihv1[0];
					y = d + ihv1[3];
					z = (c + ihv1[2]) & hybridmask;
				} else {
					x = a + ihv2mod[0];
					y = d + ihv2mod[3];
					z = (c + ihv2mod[2]) & hybridmask;
				}
			}
			else
			{
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
			}
		}

		if (len >= len2)
		{
			// processed trail1
			// write to end
			gcollision_states.get_ref<3>(idx) = x;
			gcollision_states.get_ref<4>(idx) = y;
			gcollision_states.get_ref<5>(idx) = z;
			if (len > 0)
				--len;
		}
		else
		{
			// processed trail2
			// write to end
			gcollision_states.get_ref<7+3>(idx) = x;
			gcollision_states.get_ref<7+4>(idx) = y;
			gcollision_states.get_ref<7+5>(idx) = z;
			if (len2 > 0)
				--len2;
		}

		bool done = (bad == 0)
			&& (len == 0 || len2 == 0 || 
				(
				  (gcollision_states.get<3>(idx) == gcollision_states.get<7+3>(idx))
				  && (gcollision_states.get<4>(idx) == gcollision_states.get<7+4>(idx))
				  && (gcollision_states.get<5>(idx) == gcollision_states.get<7+5>(idx))
				));

		{
			if (done)
			{
				if (len > 0) len = 1;
				if (len2 > 0) len2 = 1;
				bad = 1;
			}
			// conditionally write result and load a new one
			gcollisionsout_buf.write(gcollisionsout_ctl, done,
				gcollision_states.get<0>(idx),
				gcollision_states.get<1>(idx),
				gcollision_states.get<2>(idx),
				gcollision_states.get<3>(idx),
				gcollision_states.get<4>(idx),
				gcollision_states.get<5>(idx),
				len,
				gcollision_states.get<7+0>(idx),
				gcollision_states.get<7+1>(idx),
				gcollision_states.get<7+2>(idx),
				gcollision_states.get<7+3>(idx),
				gcollision_states.get<7+4>(idx),
				gcollision_states.get<7+5>(idx),
				len2);
		}

		if (4 <= __popc(__ballot_sync(WARP_FULL_MASK,bad)))
		{
			uint32 readidx = gcollisionsin_buf.getreadidx(gcollisionsin_ctl, bad);
			if (bad && readidx < 0xEEEEEEEE)
			{
				len  = gcollisionsin_buf.get<6>(readidx);
				len2 = gcollisionsin_buf.get<7+6>(readidx);
				gcollision_states.get_ref<0>(idx) = gcollisionsin_buf.get<0>(readidx);
				gcollision_states.get_ref<1>(idx) = gcollisionsin_buf.get<1>(readidx);
				gcollision_states.get_ref<2>(idx) = gcollisionsin_buf.get<2>(readidx);
				gcollision_states.get_ref<3>(idx) = gcollisionsin_buf.get<0>(readidx);
				gcollision_states.get_ref<4>(idx) = gcollisionsin_buf.get<1>(readidx);
				gcollision_states.get_ref<5>(idx) = gcollisionsin_buf.get<2>(readidx);
				gcollision_states.get_ref<6>(idx) = gcollisionsin_buf.get<6>(readidx);
				gcollision_states.get_ref<7+0>(idx) = gcollisionsin_buf.get<7+0>(readidx);
				gcollision_states.get_ref<7+1>(idx) = gcollisionsin_buf.get<7+1>(readidx);
				gcollision_states.get_ref<7+2>(idx) = gcollisionsin_buf.get<7+2>(readidx);
				gcollision_states.get_ref<7+3>(idx) = gcollisionsin_buf.get<7+0>(readidx);
				gcollision_states.get_ref<7+4>(idx) = gcollisionsin_buf.get<7+1>(readidx);
				gcollision_states.get_ref<7+5>(idx) = gcollisionsin_buf.get<7+2>(readidx);
				gcollision_states.get_ref<7+6>(idx) = gcollisionsin_buf.get<7+6>(readidx);
				gcollision_states.get_ref<14>(idx) = bad = 0;
			}
		}
		if (__all_sync(WARP_FULL_MASK,bad))
			break;
		if (__shfl_sync(WARP_FULL_MASK,halt_flag,0)) // read global halt flag together and halt if set
			break;
//		__syncthreads();
	}
	gcollision_states.get_ref<6>(idx) = len;
	gcollision_states.get_ref<7+6>(idx) = len2;
	gcollision_states.get_ref<14>(idx) = bad;
}



void cuda_device::cuda_fill_trail_buffer(uint32 id, uint64 seed,
							vector<trail_type>& buf,
							vector< pair<trail_type,trail_type> >& collisions, bool mod)
{
//	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->trailsout_buf)), detail->blocks * sizeof(work_buffer_t) ) );
//	CUDA_SAFE_CALL( cudaMallocHost( (void**)(&(detail->trailsout_ctl)), detail->blocks * sizeof(work_control_t) ) );

	// move all collisions into buffer
	for (auto& c : collisions)
		detail->collisions.emplace_back(c);
	collisions.clear();

	// if collisions buffer is big enough then actually launch it
	uint32 collisionblocks = 0;
	if (detail->collisions.size())
	{
		size_t oldsize = detail->collisions.size();
		// store input collisions to GPU by writing to host buffer
		// and sending it to GPU, we only move the control back and forth
		uint32 count = detail->collisions.size();
		// don't overwrite collision data still in the buffer
		if (count >= detail->collisionsin_ctl->free_count())
			count = detail->collisionsin_ctl->free_count();
		if (count > 0)
			count -= 1;
		for (std::size_t i = 0; i < count; ++i)
		{
				detail->collisionsin_buf->write(*(detail->collisionsin_ctl), true,
					detail->collisions[i].first.start[0], detail->collisions[i].first.start[1], detail->collisions[i].first.start[2],
					detail->collisions[i].first.end[0], detail->collisions[i].first.end[1], detail->collisions[i].first.end[2],
					detail->collisions[i].first.len,
					detail->collisions[i].second.start[0], detail->collisions[i].second.start[1], detail->collisions[i].second.start[2],
					detail->collisions[i].second.end[0], detail->collisions[i].second.end[1], detail->collisions[i].second.end[2],
					detail->collisions[i].second.len);
		}
		detail->collisions.erase(detail->collisions.begin(), detail->collisions.begin() + count);
		detail->nrcollisions_on_gpu += count;

		// determine how many cuda blocks to start for collision
		collisionblocks = (detail->nrcollisions_on_gpu / detail->threadsperblock)/2;
		if (collisionblocks > detail->blocks)
			collisionblocks = detail->blocks;
		// only copy data to GPU when we're actually going to run GPU code
		if (collisionblocks > 0)
		{
			// send control and buffer structures to GPU
			cudaMemcpyToSymbol(gcollisionsin_ctl, detail->collisionsin_ctl, sizeof(collisions_control_t));
			cudaMemcpyToSymbol(gcollisionsin_buf, detail->collisionsin_buf, sizeof(collisions_buffer_t));
		}
		if (0) std::cout << "C: " << oldsize << " " << detail->collisions.size() 
			<< " " << detail->collisionsin_ctl->used_count()
			<< " " << detail->collisionsin_ctl->free_count()
			<< " " << collisionblocks << " " << detail->blocks
			<< std::endl;
		
	}

	// send control structures to GPU
	cudaMemcpyToSymbol(gtrailsout_ctl, detail->trailsout_ctl, detail->blocks * sizeof(work_control_t));
	// retrieve store buffers from GPU
	cudaMemcpyToSymbol(gtrailsout_buf, detail->trailsout_buf, detail->blocks * sizeof(work_buffer_t));

	// run GPU code
	if (mod)
	{
		cuda_md5_work<true><<<detail->blocks - collisionblocks, detail->threadsperblock>>>(seed);
		cuda_md5_collisions<true><<< collisionblocks, detail->threadsperblock>>>(seed);
	}
	else
	{
		cuda_md5_work<false><<<detail->blocks - collisionblocks, detail->threadsperblock>>>(seed);
		cuda_md5_collisions<false><<< collisionblocks, detail->threadsperblock>>>(seed);
	}

	// retrieve store buffers from GPU
	cudaMemcpyFromSymbol(detail->trailsout_buf, gtrailsout_buf, detail->blocks * sizeof(work_buffer_t));
	// retrieve control structures from GPU
	cudaMemcpyFromSymbol(detail->trailsout_ctl, gtrailsout_ctl, detail->blocks * sizeof(work_control_t));

/*	std::cout << detail->trailsout_ctl[0].write_idx << " " <<
		detail->trailsout_ctl[0].read_idx << std::endl;
*/
	// if we started a collision processing cuda job then process its output
	if (collisionblocks > 0)
	{
		cudaMemcpyFromSymbol(detail->collisionsout_buf, gcollisionsout_buf, sizeof(collisions_buffer_t));
		cudaMemcpyFromSymbol(detail->collisionsin_ctl, gcollisionsin_ctl, sizeof(collisions_control_t));
		cudaMemcpyFromSymbol(detail->collisionsout_ctl, gcollisionsout_ctl, sizeof(collisions_control_t));
		uint32 readidx;
		while ((readidx=detail->collisionsout_buf->getreadidx(*(detail->collisionsout_ctl)))
			< 0xEEEEEEEE)
		{
			--detail->nrcollisions_on_gpu;
			collisions.emplace_back();
			trail_type& first = collisions.back().first;
			trail_type& second = collisions.back().second;
			first.start[0]  = detail->collisionsout_buf->get<0>(readidx);
			first.start[1]  = detail->collisionsout_buf->get<1>(readidx);
			first.start[2]  = detail->collisionsout_buf->get<2>(readidx);
			first.end[0]    = detail->collisionsout_buf->get<3>(readidx);
			first.end[1]    = detail->collisionsout_buf->get<4>(readidx);
			first.end[2]    = detail->collisionsout_buf->get<5>(readidx);
			first.len       = detail->collisionsout_buf->get<6>(readidx);
			second.start[0] = detail->collisionsout_buf->get<7+0>(readidx);
			second.start[1] = detail->collisionsout_buf->get<7+1>(readidx);
			second.start[2] = detail->collisionsout_buf->get<7+2>(readidx);
			second.end[0]   = detail->collisionsout_buf->get<7+3>(readidx);
			second.end[1]   = detail->collisionsout_buf->get<7+4>(readidx);
			second.end[2]   = detail->collisionsout_buf->get<7+5>(readidx);
			second.len      = detail->collisionsout_buf->get<7+6>(readidx);
		}
		cudaMemcpyToSymbol(gcollisionsout_ctl, detail->collisionsout_ctl, sizeof(collisions_control_t));
	}

	// process and return results
	buf.clear();
	for (unsigned b = 0; b < detail->blocks; ++b)
	{
		uint32 readidx;
		trail_type trail;
		while ((readidx=detail->trailsout_buf[b].getreadidx(detail->trailsout_ctl[b])) != 0xFFFFFFFF)
		{
			trail.start[0] = detail->trailsout_buf[b].get<0>(readidx);
			trail.start[1] = detail->trailsout_buf[b].get<1>(readidx);
			trail.start[2] = detail->trailsout_buf[b].get<2>(readidx);
			trail.end[0]   = detail->trailsout_buf[b].get<3>(readidx);
			trail.end[1]   = detail->trailsout_buf[b].get<4>(readidx);
			trail.end[2]   = detail->trailsout_buf[b].get<5>(readidx);
			trail.len      = detail->trailsout_buf[b].get<6>(readidx);
			buf.push_back(trail);
		}
	}
//	std::cout << "B " << buf.size() << std::endl;
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
/*
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
*/
}














int get_num_cuda_devices()
{
	int deviceCount = 0;
	cutilSafeCall(cudaGetDeviceCount(&deviceCount));
	return deviceCount;
}

void cuda_device_query()
{
    int deviceCount = 0;
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
