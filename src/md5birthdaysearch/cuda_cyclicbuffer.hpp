/*****
  Copyright (C) 2015 Pierre Karpman, INRIA France/Nanyang Technological University Singapore (-2016), CWI (2016/2017), L'Universite Grenoble Alpes (2017-)
            (C) 2015 Marc Stevens, Centrum Wiskunde & Informatica (CWI), Amsterdam.

  This file is part of sha1_gpu_nearcollisionattacks source-code and released under the MIT License
*****/

#include <stdint.h>

#define WARP_FULL_MASK 0xFFFFFFFF
// fencetype
// - 0 none
// - 1 block-wide
// - 2 gpu-wide
template<int fencetype> __device__ void memoryfence() { __threadfence(); }
template<>              __device__ void memoryfence<0>() {}
template<>              __device__ void memoryfence<1>() { __threadfence_block(); }

//// cyclic buffer mask version control struct: to be placed in fast CUDA 'shared' memory
// template<size_t N> class cyclic_buffer_control_mask_t;

//// cyclic buffer mask version storage struct: to be placed in main CUDA 'global' memory
// template<size_t N, typename val_t = uint32_t, size_t val_cnt = 1, typename control_type = cyclic_buffer_control_mask_t<N>, int fencetype = 2 >
// class cyclic_buffer_mask_t;

//// cyclic buffer cas version control struct: to be placed in fast CUDA 'shared' memory
// template<size_t N> class cyclic_buffer_control_cas_t;

//// cyclic buffer cas version storage struct: to be placed in main CUDA 'global' memory
// template<size_t N, typename val_t = uint32_t, size_t val_cnt = 1, typename control_type = cyclic_buffer_control_cas_t<N>, int fencetype = 2 >
// class cyclic_buffer_cas_t;

//// temporary buffer of 64*2 words in fast CUDA 'shared' memory
//// to ensure writes to cyclic buffer (with N=1,2) happen in blocks of 32 words
// class warp_tmp_buf_t;


/* UTILITY TYPES */
template<bool b>
struct assert_compile_time_bool
{
};
template<>
struct assert_compile_time_bool<true>
{
	typedef bool compile_time_bool_is_true;
};
#define ASSERT_COMPILE_TIME(b) { typename assert_compile_time_bool< (b) >::compile_time_bool_is_true ___assert = true; if (!___assert) { printf("boooh!"); } }

template<size_t i, bool pred = true>
struct pred_t {};



/******************************************** MASK VERSION *****************************************/
/******************************************** MASK VERSION *****************************************/
/******************************************** MASK VERSION *****************************************/

// control logic using mask version
template<size_t N>
class cyclic_buffer_control_mask_t
{
	public:

	static const uint32_t size = N;
	volatile uint32_t write_idx;
	volatile uint32_t read_idx;

	__host__ __device__ inline void reset()
	{
		write_idx = 0;
		read_idx = 0;
	}
	__host__ __device__ inline uint32_t used_count()
	{
		return write_idx - read_idx;
	}
	__host__ __device__ inline uint32_t free_count()
	{
		uint32_t c = (N + read_idx - write_idx);
		return c==0 ? N : c;
	}

	// assumes entire warp calls this function with identical warp_to_write_mask
	// returns per-thread idx to write to
	__device__ inline uint32_t warp_write_idx(uint32_t warp_to_write_mask)
	{
		// warp: determine count and offset
		uint32_t count  = __popc(warp_to_write_mask);
		// thread ZERO has offset 0
		uint32_t offset = count - __popc(warp_to_write_mask >> (threadIdx.x&31));
		uint32_t baseidx = 0;
		// thread ZERO increases write_idx with count and obtains old value in baseidx
		if ((threadIdx.x&31)==0)
		{
			baseidx = atomicAdd((uint32_t*)&write_idx,count);
		}
		// thread ZERO passes baseidx to entire warp
		baseidx = __shfl_sync(WARP_FULL_MASK,baseidx,0);
		return (baseidx + offset) % N;
	}

	// assumes entire warp calls this function with identical inputs
	// if read is safe return per-thread idx to read from, returns 0xFFFFFFFF otherwise
	template<typename safereadval_t>
	__device__ inline uint32_t warp_read_idx(safereadval_t safereadvals[N])
	{
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx,0);
			uint32_t thrdidx = (baseidx+(threadIdx.x&31)) % N;
			// if not safe to read then return false
			if (!__all_sync(WARP_FULL_MASK,0 != safereadvals[ thrdidx ] ))
			{
				return 0xFFFFFFFF;
			}
			// thread ZERO tries to increase read_idx with 32, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x&31)==0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidx+32);
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// if read_idx was successfully increased return (true,baseidx+offset)
			if (baseidx2 == baseidx)
			{
				safereadvals[thrdidx] = 0;
				return thrdidx;
			}
			// increase failed due to race, try again
		}
	}

	// assumes entire warp calls this function with identical inputs
	// if read is safe return per-thread idx to read from, returns 0xFFFFFFFF otherwise
	template<typename safereadval_t>
	__device__ inline uint32_t warp_read_idx(uint32_t warp_to_read_mask, safereadval_t safereadvals[N])
	{
		// warp: determine count and offset
		uint32_t count  = __popc(warp_to_read_mask);
		// thread ZERO has offset 0
		uint32_t offset = count - __popc(warp_to_read_mask >> (threadIdx.x&31));
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx,0);
			uint32_t thrdidx = (baseidx+offset) % N;
			// if not safe to read then return false
			if (!__all_sync(WARP_FULL_MASK,0 != safereadvals[ thrdidx ] ))
			{
				return 0xFFFFFFFF;
			}
			// thread ZERO tries to increase read_idx with count, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x&31)==0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidx+count);
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// if read_idx was successfully increased return (true,baseidx+offset)
			if (baseidx2 == baseidx)
			{
				safereadvals[thrdidx] = 0;
				return thrdidx;
			}
			// increase failed due to race, try again
		}
	}

	
	// variant partial warp read
	// assumes entire warp calls this function with identical inputs
	// if there is at least 1 element and read is safe return per-thread idx to read from or 0xEEEEEEEE if there was no input for that thread, returns 0xFFFFFFFF otherwise (no work at all)
	template<typename safereadval_t>
	__device__ inline uint32_t warp_read_any_idx(safereadval_t safereadvals[N])
	{
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx, 0);
			uint32_t thrdidx = (baseidx + (threadIdx.x & 31)) % N;
			bool issafetoread = (0 != safereadvals[thrdidx]);
			uint32_t safereadmask = __ballot_sync(WARP_FULL_MASK,issafetoread);
			if (safereadmask == 0)
				return 0xFFFFFFFF;
			// safereadmask needs to be all 1-bits for some lower half, the remaining upper half needs to be all 0-bits
			// otherwise there is some work but not yet safe to read, simply check again
			if ((safereadmask | (safereadmask >> 1)) != safereadmask)
				continue;

			// thread ZERO tries to increase read_idx with 32, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x & 31) == 0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidx + __popc(safereadmask));
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// increase failed due to race, try again
			if (baseidx2 != baseidx)
				continue;
			// if read_idx was successfully increased return baseidx+offset or 0xEEEEEEEE
			if (issafetoread)
			{
				safereadvals[thrdidx] = 0;
				return thrdidx;
			}
			else
				return 0xEEEEEEEE;
		}
	}

	/**** HOST FUNCTIONS ASSUME EXCLUSIVE ACCESS FROM 1 HOST THREAD ****/
	template<typename safereadval_t>
	__host__ inline uint32_t host_read_idx(safereadval_t safereadvals[N])
	{
		uint32_t i = read_idx % N;
		if (safereadvals[i])
		{
			safereadvals[i] = 0;
			++read_idx;
			return i;
		}
		else
		{
			return 0xFFFFFFFF;
		}
	}
	__host__ inline uint32_t host_write_idx()
	{
		uint32_t i = write_idx % N;
		++write_idx;
		return i;
	}
};



template<size_t N, typename val_t = uint32_t, size_t val_cnt = 1, typename control_type = cyclic_buffer_control_mask_t<N>, int fencetype = 2 >
class cyclic_buffer_mask_t
{
	public:

	typedef control_type control_t;

	static const uint32_t size = N;
	static const uint32_t val_size = val_cnt;

	val_t val[val_cnt][N];
	volatile char safe_to_read[N];
	
	template<size_t i, bool pred>
	__host__ __device__ inline void cond_write2(uint32_t j, const val_t& v, pred_t<i, pred>)
	{
		val[i][j % N] = v;
	}
	template<size_t i>
	__host__ __device__ inline void cond_write2(uint32_t j, const val_t& v, pred_t<i, false>)
	{
	}
	template<size_t i, bool pred>
	__host__ __device__ inline void cond_write(uint32_t j, const val_t& v, pred_t<i, pred>)
	{
		cond_write2(j, v, pred_t<i,(i<val_cnt)>());
	}
	
	// called by entire warp
	__host__ __device__ inline void reset(control_t& control)
	{
#ifdef __CUDA_ARCH__
		for (uint32_t i = (threadIdx.x&31); i < N; i+=32)
			safe_to_read[i] = 0;
#else
		for (uint32_t i = 0; i < N; ++i)
			safe_to_read[i] = 0;
#endif
		control.reset();
	}

	// device: called by entire warp
	// host  : assumes exclusive access !!
	__host__ __device__ inline void write(control_t& control, bool dowrite
		, val_t _val0   , val_t _val1 =0, val_t _val2 =0, val_t _val3 =0, val_t _val4 =0, val_t _val5 =0, val_t _val6 =0, val_t _val7 =0, val_t _val8 =0, val_t _val9 =0
		, val_t _val10=0, val_t _val11=0, val_t _val12=0, val_t _val13=0, val_t _val14=0, val_t _val15=0, val_t _val16=0, val_t _val17=0, val_t _val18=0, val_t _val19=0
		, val_t _val20=0, val_t _val21=0, val_t _val22=0, val_t _val23=0, val_t _val24=0, val_t _val25=0, val_t _val26=0, val_t _val27=0, val_t _val28=0, val_t _val29=0
		, val_t _val30=0, val_t _val31=0
		)
	{
#ifdef __CUDA_ARCH__
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,dowrite);
		if (mask == 0) 
		{
			return;
		}
		uint32_t wi = control.warp_write_idx(mask);
		if (dowrite)
		{
#else
		if (dowrite)
		{
			uint32_t wi = control.host_write_idx();
#endif
			cond_write(wi, _val0, pred_t<0>());
			cond_write(wi, _val1, pred_t<1>());
			cond_write(wi, _val2, pred_t<2>());
			cond_write(wi, _val3, pred_t<3>());
			cond_write(wi, _val4, pred_t<4>());
			cond_write(wi, _val5, pred_t<5>());
			cond_write(wi, _val6, pred_t<6>());
			cond_write(wi, _val7, pred_t<7>());
			cond_write(wi, _val8, pred_t<8>());
			cond_write(wi, _val9, pred_t<9>());
			cond_write(wi, _val10, pred_t<10>());
			cond_write(wi, _val11, pred_t<11>());
			cond_write(wi, _val12, pred_t<12>());
			cond_write(wi, _val13, pred_t<13>());
			cond_write(wi, _val14, pred_t<14>());
			cond_write(wi, _val15, pred_t<15>());
			cond_write(wi, _val16, pred_t<16>());
			cond_write(wi, _val17, pred_t<17>());
			cond_write(wi, _val18, pred_t<18>());
			cond_write(wi, _val19, pred_t<19>());
			cond_write(wi, _val20, pred_t<20>());
			cond_write(wi, _val21, pred_t<21>());
			cond_write(wi, _val22, pred_t<22>());
			cond_write(wi, _val23, pred_t<23>());
			cond_write(wi, _val24, pred_t<24>());
			cond_write(wi, _val25, pred_t<25>());
			cond_write(wi, _val26, pred_t<26>());
			cond_write(wi, _val27, pred_t<27>());
			cond_write(wi, _val28, pred_t<28>());
			cond_write(wi, _val29, pred_t<29>());
			cond_write(wi, _val30, pred_t<30>());
			cond_write(wi, _val31, pred_t<31>());
#ifdef __CUDA_ARCH__
			memoryfence<fencetype>();
#endif
			safe_to_read[wi % N] = 1;
		}
	}

	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe
	__host__ __device__ inline uint32_t getreadidx(control_t& control)
	{
#ifdef __CUDA_ARCH__
		uint32_t ret = control.warp_read_idx(safe_to_read);
		if (ret != 0xFFFFFFFF)
			memoryfence<fencetype>();
		return ret;
#else
		return control.host_read_idx(safe_to_read);
#endif
	}

	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe
	__host__ __device__ inline uint32_t getreadidx(control_t& control, bool doread)
	{
#ifdef __CUDA_ARCH__
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,doread);
		if (mask == 0) 
			return 0xFFFFFFFF;
		uint32_t ret = control.warp_read_idx(mask,safe_to_read);
		if (ret != 0xFFFFFFFF)
			memoryfence<fencetype>();
		return doread ? ret : 0xFFFFFFFF;
#else
		if (!doread)
			return 0xFFFFFFFF;
		return control.host_read_idx(safe_to_read);
#endif
	}

	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe
	__host__ __device__ inline uint32_t getreadanyidx(control_t& control)
	{
#ifdef __CUDA_ARCH__
		return control.warp_read_any_idx(safe_to_read);
#else
		return control.host_read_any_idx(safe_to_read);
#endif
	}


	template<size_t i>
	__host__ __device__ inline val_t get(uint32_t idx)
	{
		ASSERT_COMPILE_TIME(i<val_cnt);
		return val[i][idx % N];
	}
	template<size_t i>
	__host__ __device__ inline val_t& get_ref(uint32_t idx)
	{
		ASSERT_COMPILE_TIME(i<val_cnt);
		return val[i][idx % N];
	}
};

















/******************************************** CAS VERSION *****************************************/
/******************************************** CAS VERSION *****************************************/
/******************************************** CAS VERSION *****************************************/

// control logic using CAS version
template<size_t N>
class cyclic_buffer_control_cas_t
{
	public:

	static const uint32_t size = N;
	volatile uint32_t write_idx;
	volatile uint32_t written_idx;
	volatile uint32_t read_idx;

	__host__ __device__ inline void reset()
	{
		write_idx = 0;
		written_idx = 0;
		read_idx = 0;
	}
	__host__ __device__ inline uint32_t used_count()
	{
		return write_idx - read_idx;
	}
	__host__ __device__ inline uint32_t free_count()
	{
		uint32_t c = (N + read_idx - write_idx);
		return c;
//		return c==0 ? N : c;
	}

	// assumes entire warp calls this function with identical warp_to_write_mask
	// returns per-thread idx to write to
	__device__ inline uint32_t warp_write_idx(uint32_t warp_to_write_mask)
	{
		// warp: determine count and offset
		uint32_t count  = __popc(warp_to_write_mask);
		// thread ZERO has offset 0
		uint32_t offset = count - __popc(warp_to_write_mask >> (threadIdx.x&31));
		uint32_t baseidx = 0;
		// thread ZERO increases write_idx with count and obtains old value in baseidx
		if ((threadIdx.x&31)==0)
		{
			baseidx = atomicAdd((uint32_t*)&write_idx,count);
		}
		// thread ZERO passes baseidx to entire warp
		baseidx = __shfl_sync(WARP_FULL_MASK,baseidx,0);
		return (baseidx + offset) % N;
	}

	// assumes entire warp calls this function with their per-thread idx to write to
	__device__ inline void warp_write_finish(uint32_t idx, uint32_t mask)
	{
		uint32_t count = __popc(mask);
		// thread ZERO updates written_idx
		// idx was old value of write_idx
		// waits till written_idx has old value of write_idx
		// and then increases it by count
		if ((threadIdx.x&31)==0)
		{
			while (idx != atomicCAS((uint32_t*)&written_idx, idx, idx+count) )
				;
		}
	}

	// assumes entire warp calls this function with identical inputs
	// if read is safe return per-thread idx to read from, returns 0xFFFFFFFF otherwise
	__device__ inline uint32_t warp_read_idx()
	{
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx,0);
			uint32_t baseidxdist = __shfl_sync(WARP_FULL_MASK,written_idx,0);
			// if not safe to read then return false
			if ( baseidxdist-baseidx < 32)
			{
				return 0xFFFFFFFF;
			}
			uint32_t thrdidx = (baseidx+(threadIdx.x&31));
			// thread ZERO tries to increase read_idx with 32, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x&31)==0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidx+32);
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// if read_idx was successfully increased return (true,baseidx+offset)
			if (baseidx2 == baseidx)
			{
				return thrdidx % N;
			}
			// increase failed due to race, try again
		}
	}

	// assumes entire warp calls this function with identical inputs
	// if read is safe return per-thread idx to read from, returns 0xFFFFFFFF otherwise
	__device__ inline uint32_t warp_read_idx(uint32_t warp_to_read_mask)
	{
		// warp: determine count and offset
		uint32_t count  = __popc(warp_to_read_mask);
		// thread ZERO has offset 0
		uint32_t offset = count - __popc(warp_to_read_mask >> (threadIdx.x&31));
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx,0);
			uint32_t baseidxdist = __shfl_sync(WARP_FULL_MASK,written_idx,0);
			// if not safe to read then return false
			if ( baseidxdist-baseidx < count)
			{
				return 0xFFFFFFFF;
			}
			uint32_t thrdidx = (baseidx+offset);
			// thread ZERO tries to increase read_idx with 32, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x&31)==0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidx+count);
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// if read_idx was successfully increased return (true,baseidx+offset)
			if (baseidx2 == baseidx)
			{
				return thrdidx % N;
			}
			// increase failed due to race, try again
		}
	}

	// variant partial warp read
	// assumes entire warp calls this function with identical inputs
	// if read is safe return per-thread idx to read from, 0xEEEEEEEE for no thread-specific input, or returns 0xFFFFFFFF otherwise (no work at all)
	__device__ inline uint32_t warp_read_any_idx()
	{
		while (true)
		{
			// obtain warp-wide unique race-free value of read_idx
			uint32_t baseidx = __shfl_sync(WARP_FULL_MASK,read_idx, 0);
			uint32_t baseidxdist = __shfl_sync(WARP_FULL_MASK,written_idx, 0);
			if (baseidxdist == baseidx)
				return 0xFFFFFFFF;

			uint32_t thrdidx = (baseidx + (threadIdx.x & 31));
			// thread ZERO tries to increase read_idx with 32, baseidx2==baseidx in case of success, != otherwise
			uint32_t baseidx2;
			if ((threadIdx.x & 31) == 0)
			{
				baseidx2 = atomicCAS((uint32_t*)&read_idx, baseidx, baseidxdist);
			}
			// thread ZERO passes baseidx2 to entire warp
			baseidx2 = __shfl_sync(WARP_FULL_MASK,baseidx2, 0);
			// increase failed due to race, try again
			if (baseidx2 != baseidx)
				continue;
			// if read_idx was successfully increased return (true,baseidx+offset)
			if (thrdidx < baseidxdist)
				return thrdidx % N;
			return 0xEEEEEEEE;
		}
	}


	/**** HOST FUNCTIONS ASSUME EXCLUSIVE ACCESS FROM 1 HOST THREAD ****/
	__host__ inline uint32_t host_write_idx()
	{
		uint32_t i = write_idx % N;
		++write_idx;
		return i;
	}
	__host__ inline void host_write_finished()
	{
		++written_idx;
	}
	__host__ inline uint32_t host_read_idx()
	{
		if (read_idx == write_idx)
		{
			return 0xFFFFFFFF;
		}
		uint32_t i = read_idx % N;
		++read_idx;
		return i;
	}
};



template<size_t N, typename val_t = uint32_t, size_t val_cnt = 1, typename control_type = cyclic_buffer_control_cas_t<N>, int fencetype = 2 >
class cyclic_buffer_cas_t
{
	public:

	typedef control_type control_t;

	static const uint32_t size = N;
	static const uint32_t val_size = val_cnt;

	val_t val[val_cnt][N];

	template<size_t i, bool pred>
	__host__ __device__ inline void cond_write2(uint32_t j, const val_t& v, pred_t<i, pred>)
	{
		val[i][j % N] = v;
	}
	template<size_t i>
	__host__ __device__ inline void cond_write2(uint32_t j, const val_t& v, pred_t<i, false>)
	{
	}
	template<size_t i, bool pred>
	__host__ __device__ inline void cond_write(uint32_t j, const val_t& v, pred_t<i, pred>)
	{
		cond_write2(j, v, pred_t<i,(i<val_cnt)>());
	}
	
	// called by entire warp
	__host__ __device__ inline void reset(control_t& control)
	{
		control.reset();
	}

	// called by entire warp
	__host__ __device__ inline void write(control_t& control, bool dowrite
		, val_t _val0   , val_t _val1 =0, val_t _val2 =0, val_t _val3 =0, val_t _val4 =0, val_t _val5 =0, val_t _val6 =0, val_t _val7 =0, val_t _val8 =0, val_t _val9 =0
		, val_t _val10=0, val_t _val11=0, val_t _val12=0, val_t _val13=0, val_t _val14=0, val_t _val15=0, val_t _val16=0, val_t _val17=0, val_t _val18=0, val_t _val19=0
		, val_t _val20=0, val_t _val21=0, val_t _val22=0, val_t _val23=0, val_t _val24=0, val_t _val25=0, val_t _val26=0, val_t _val27=0, val_t _val28=0, val_t _val29=0
		, val_t _val30=0, val_t _val31=0
		)
	{
#ifdef __CUDA_ARCH__
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,dowrite);
		if (mask == 0) 
		{
			return;
		}
		uint32_t wi = control.warp_write_idx(mask);
		if (dowrite)
		{
#else
		if (dowrite)
		{
			uint32_t wi = control.host_write_idx();
#endif
			cond_write(wi, _val0, pred_t<0>());
			cond_write(wi, _val1, pred_t<1>());
			cond_write(wi, _val2, pred_t<2>());
			cond_write(wi, _val3, pred_t<3>());
			cond_write(wi, _val4, pred_t<4>());
			cond_write(wi, _val5, pred_t<5>());
			cond_write(wi, _val6, pred_t<6>());
			cond_write(wi, _val7, pred_t<7>());
			cond_write(wi, _val8, pred_t<8>());
			cond_write(wi, _val9, pred_t<9>());
			cond_write(wi, _val10, pred_t<10>());
			cond_write(wi, _val11, pred_t<11>());
			cond_write(wi, _val12, pred_t<12>());
			cond_write(wi, _val13, pred_t<13>());
			cond_write(wi, _val14, pred_t<14>());
			cond_write(wi, _val15, pred_t<15>());
			cond_write(wi, _val16, pred_t<16>());
			cond_write(wi, _val17, pred_t<17>());
			cond_write(wi, _val18, pred_t<18>());
			cond_write(wi, _val19, pred_t<19>());
			cond_write(wi, _val20, pred_t<20>());
			cond_write(wi, _val21, pred_t<21>());
			cond_write(wi, _val22, pred_t<22>());
			cond_write(wi, _val23, pred_t<23>());
			cond_write(wi, _val24, pred_t<24>());
			cond_write(wi, _val25, pred_t<25>());
			cond_write(wi, _val26, pred_t<26>());
			cond_write(wi, _val27, pred_t<27>());
			cond_write(wi, _val28, pred_t<28>());
			cond_write(wi, _val29, pred_t<29>());
			cond_write(wi, _val30, pred_t<30>());
			cond_write(wi, _val31, pred_t<31>());
#ifdef __CUDA_ARCH__
			memoryfence<fencetype>();
		}
		control.warp_write_finish(wi,mask);
#else
			control.host_write_finished();
		}
#endif
	}

	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe
	__host__ __device__ inline uint32_t getreadidx(control_t& control)
	{
#ifdef __CUDA_ARCH__
		uint32_t ret = control.warp_read_idx();
		if (ret != 0xFFFFFFFF)
			memoryfence<fencetype>();
		return ret;
#else
		return control.host_read_idx();
#endif
	}


	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe
	__host__ __device__ inline uint32_t getreadidx(control_t& control, bool doread)
	{
#ifdef __CUDA_ARCH__
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,doread);
		if (mask == 0) 
			return 0xFFFFFFFF;
		uint32_t ret = control.warp_read_idx(mask);
		if (ret != 0xFFFFFFFF)
			memoryfence<fencetype>();
		return doread ? ret : 0xFFFFFFFF;
#else
		if (!doread)
			return 0xFFFFFFFF;
		return control.host_read_idx();
#endif
	}

	// called by entire warp
	// returns 0xFFFFFFFF if read is not possible
	// returns per-thread read index if read is safe or 0xEEEEEEEE if no read for that thread
	__host__ __device__ inline uint32_t getreadanyidx(control_t& control)
	{
#ifdef __CUDA_ARCH__
		return control.warp_read_any_idx();
#else
		return control.host_read_idx();
#endif
	}

	template<size_t i>
	__host__ __device__ inline val_t get(uint32_t idx)
	{
		ASSERT_COMPILE_TIME(i<val_cnt);
		return val[i][idx % N];
	}
	template<size_t i>
	__host__ __device__ inline val_t& get_ref(uint32_t idx)
	{
		ASSERT_COMPILE_TIME(i<val_cnt);
		return val[i][idx % N];
	}
};






/******************************************** TEMPORARY BUFFER *****************************************/
/******************************************** TEMPORARY BUFFER *****************************************/
/******************************************** TEMPORARY BUFFER *****************************************/

class warp_tmp_buf_t
{
	public:
	// temporary buffer for 64 2-word elems: two halves of 32 elems
	// whenever one half is full it is flushed to the main buffer
	uint32_t val1[64];
	uint32_t val2[64];
	volatile char idx;

	__device__ void reset()
	{
		idx = 0;
	}

	// store 1-word elements in temporary shared buffer, flush to global buffer for each 32 values
	template<typename buffer_t, typename control_t>
	__device__ inline void write1(bool dowrite, uint32_t _val1, buffer_t& buf, control_t& ctrl)
	{
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,dowrite);
		if (mask == 0)
		{
			return;
		}

		uint32_t count = __popc(mask);
		uint32_t offset = count - __popc(mask >> (threadIdx.x&31));
		uint32_t baseidx = idx;
		if ((threadIdx.x&31)==0)
		{
			idx = (idx + (char)(count)) % 64;
		}

		if (dowrite)
		{
			val1[(baseidx + offset) % 64] = _val1;
		}

		// flush 'full' half if we cross over halves
		if ((idx^baseidx)&32)
		{
			baseidx &= 32; // point to start of 'full' half
			buf.write(ctrl, true, val1[baseidx + (threadIdx.x&31)] );
		}
	}
	
	template<class buffer_t, class control_t>
	__device__ inline void write2(bool dowrite, uint32_t _val1, uint32_t _val2, buffer_t& buf, control_t& ctrl)
	{
		uint32_t mask = __ballot_sync(WARP_FULL_MASK,dowrite);
		if (mask == 0)
		{
			return;
		}

		uint32_t count = __popc(mask);
		uint32_t offset = count - __popc(mask >> (threadIdx.x&31));
		uint32_t baseidx = idx;
		if ((threadIdx.x&31)==0)
		{
			idx = (idx + (char)(count)) % 64;
		}

		if (dowrite)
		{
			val1[(baseidx + offset) % 64] = _val1;
			val2[(baseidx + offset) % 64] = _val2;
		}

		// flush 'full' half if we cross over halves
		if ((idx^baseidx)&32)
		{
			baseidx &= 32; // point to start of 'full' half
			buf.write(ctrl, true, val1[baseidx + (threadIdx.x&31)], val2[baseidx + (threadIdx.x&31)] );
		}
	}

	template<typename buffer_t, typename control_t>
	__device__ inline void flush1(buffer_t& buf, control_t& ctrl)
	{
		uint32_t baseidx = (uint32_t)(idx) & 32;
		if (idx != baseidx)
		{
			buf.write(ctrl, (threadIdx.x&31)<(idx-baseidx), val1[baseidx + (threadIdx.x&31)] );
		}
		reset();
	}		

	template<typename buffer_t, typename control_t>
	__device__ inline void flush2(buffer_t& buf, control_t& ctrl)
	{
		uint32_t baseidx = (uint32_t)(idx) & 32;
		if (idx != baseidx)
		{
			buf.write(ctrl, (threadIdx.x&31)<(idx-baseidx), val1[baseidx + (threadIdx.x&31)], val2[baseidx + (threadIdx.x&31)] );
		}
		reset();
	}		


};
