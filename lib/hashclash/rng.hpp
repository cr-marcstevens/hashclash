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

/* 

very fast inlined xorshift random number generators with periods 2^32 - 1, 2^64 - 1, 2^96 - 1 and 2^128 - 1
by G. Marsaglia: http://www.jstatsoft.org/v08/i14/xorshift.pdf 

*/

#ifndef HASHCLASH_XORSHIFT_RNG_HPP
#define HASHCLASH_XORSHIFT_RNG_HPP

#include "types.hpp"

namespace hashclash {

	
	// seed all generators using 32-bit values
	// seed state is deterministically dependent on the given values
	void seed(uint32 s);
	void seed(const uint32* sbuf, unsigned len);
	// add seed to the generators
	// seed state is changed by given values
	void addseed(uint32 s);
	void addseed(const uint32* sbuf, unsigned len);

	// seeds used, these are initialized to random values based on the time
	extern uint32 seedd;
	extern uint32 seed32_1;
	extern uint32 seed32_2;
	extern uint32 seed32_3;
	extern uint32 seed32_4;

	/******** Random generator with perdiod (2^32 - 1)*2^32 **********/
	inline uint32 xrng32()
	{
		seed32_1 ^= seed32_1 << 13;
		seed32_1 ^= seed32_1 >> 17;
		return (seed32_1 ^= seed32_1 << 5) + (seedd += 789456123);
	}

	/******** Random generator with perdiod (2^64 - 1)*2^32 **********/
	inline uint32 xrng64()
	{
		uint32 t = seed32_1 ^ (seed32_1 << 10);
		seed32_1 = seed32_2;
		seed32_2 = (seed32_2^(seed32_2>>10))^(t^(t>>13));
		return seed32_1 + (seedd += 789456123);
	}

	/******** Random generator with perdiod (2^96 - 1)*2^32 **********/
	inline uint32 xrng96()
	{
		uint32 t = seed32_1 ^ (seed32_1 << 10);
		seed32_1 = seed32_2;
		seed32_2 = seed32_3;
		seed32_3 = (seed32_3^(seed32_3>>26))^(t^(t>>5));
		return seed32_1 + (seedd += 789456123);
	}

	/******** Random generator with perdiod (2^128 - 1)*2^32 **********/
	inline uint32 xrng128()
	{
		uint32 t = seed32_1 ^ (seed32_1 << 5);
		seed32_1 = seed32_2;
		seed32_2 = seed32_3;
		seed32_3 = seed32_4;
		seed32_4 = (seed32_4^(seed32_4>>1))^(t^(t>>14));
		return seed32_1 + (seedd += 789456123);
	}

	// call this in a other global constructor using any xrng (without seeding itself)
	// due to the inpredictable order of which global constructors are called
	// otherwise the correct seeding of xrng is not guaranteed
	void hashclash_rng_hpp_init();

} // namespace hash

#endif // HASHCLASH_XORSHIFT_RNG_HPP

