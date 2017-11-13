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

#ifndef HASHCLASH_SIGNED_DIGIT_REPRESENTATION_HPP
#define HASHCLASH_SIGNED_DIGIT_REPRESENTATION_HPP

#include <vector>
#include <ostream>
#include <functional>

#ifndef NOSERIALIZATION
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#endif // NOSERIALIZATION

#include "types.hpp"

namespace hashclash {

	// A class that holds a binary Signed Digit Representation (SDR) (modulo 2^32):
	class sdr;
	// calculate hamming weight:
	inline unsigned hw(uint32 n);
//	inline unsigned hw(unsigned n) { return hw(uint32(n)); }
	inline unsigned hw(int n) { return hw(uint32(n)); }
	inline unsigned hw(uint64 n) {
		return hw(uint32(n>>32)) + hw(uint32(n));
	}
	// calculate hamming weight of the NAF of n - is the lowest hamming weight among all SDR's:
	inline unsigned hwnaf(uint32 n);
	// calculate the NAF of n:
	sdr naf(uint32 n);

	// when rotating a difference, 4 possible differences exist after rotation. Use this function to find them:
	// returns possible differences with probability >= than minprob/2^32
	void rotate_difference(uint32 diff, int rc, std::vector<uint32>& rotateddiff, uint32 minprob = uint32(1)<<30);
	void rotate_difference(uint32 diff, int rc, std::vector<std::pair<uint32,double> >& rotateddiff);
	uint32 best_rotated_difference(uint32 diff, int rc);

	class sdr {
	public:
		/**** constructors ****/
		sdr(): mask(0), sign(0) {}
		sdr(uint32 n): mask(n), sign(n) {}
		sdr(uint32 n1, uint32 n2): mask(n1^n2), sign(n2&~n1) {}
		sdr(const sdr& r): mask(r.mask), sign(r.sign) {}

		/**** assign functions/operators ****/
		sdr& set(const sdr& r)
		{
			mask = r.mask;
			sign = r.sign;
			return *this;
		}
		sdr& set(uint32 n)
		{
			mask = sign = n;
			return *this;
		}
		sdr& set(uint32 n1, uint32 n2)
		{
			mask = n1^n2;
			sign = n2 & ~n1;
			return *this;
		}
		sdr& operator= (const sdr& r)
		{ return set(r); }
		sdr& operator= (const uint32 n)
		{ return set(n); }
		sdr& clear()
		{ return set(0); }		
		
		/**** compare operators ****/
		bool operator== (const sdr& r) const
		{ return mask == r.mask && sign == r.sign; }
		bool operator< (const sdr& r) const
		{ return mask < r.mask || (mask == r.mask && sign < r.sign); }
		bool operator!= (const sdr& r) const
		{ return !(*this == r); }
		bool operator> (const sdr& r) const
		{ return r < *this; }
		bool operator<= (const sdr& r) const
		{ return !(*this > r); }
		bool operator>= (const sdr& r) const
		{ return !(*this < r); }

		/**** arithmetic operators ****/
		sdr operator-() const
		{
			sdr tmp(*this);
			tmp.sign ^= tmp.mask;
			return tmp;
		}

		sdr& operator+= (const sdr& r)
		{
			uint32 set1 = set1conditions() + r.set1conditions();
			return set(set1, set1 + adddiff() + r.adddiff());
		}
		sdr operator+ (const sdr& r) const
		{ return sdr(*this) += r; }

		sdr& operator-= (const sdr& r)
		{
			uint32 set1 = set1conditions() + r.sign;
			return set(set1, set1 + adddiff() - r.adddiff());
		}
		sdr operator- (const sdr& r) const
		{ return sdr(*this) -= r; }

		sdr& operator^= (const sdr& r)
		{
			mask ^= r.mask;
			sign = mask & ~(sign^r.sign);
			return *this;
		}
		sdr operator^ (const sdr& r) const
		{ return sdr(*this) ^= r; }

		sdr& operator<<= (unsigned n)
		{
			mask <<= n;
			sign <<= n;
			return *this;
		}
		sdr operator<< (unsigned n) const
		{ return sdr(*this) <<= n; }

		sdr& operator>>= (unsigned n)
		{
			mask >>= n;
			sign >>= n;
			return *this;
		}
		sdr operator>> (unsigned n) const
		{ return sdr(*this) >>= n; }

		sdr rotate_left(unsigned n) const
		{
			sdr tmp(*this);
			tmp.mask = hashclash::rotate_left(mask, n);
			tmp.sign = hashclash::rotate_left(sign, n);
			return tmp;
		}
		sdr rotate_right(unsigned n) const
		{
			sdr tmp(*this);
			tmp.mask = hashclash::rotate_right(mask, n);
			tmp.sign = hashclash::rotate_right(sign, n);
			return tmp;
		}

		/**** other functions ****/
		uint32 adddiff() const
		{ return sign - (sign^mask); }
		uint32 xordiff() const
		{ return mask; }
		uint32 set0conditions() const
		{ return ~sign; }
		uint32 set1conditions() const
		{ return sign^mask; }

		int get(unsigned b) const 
		{
			return int((sign>>b)&1) - int(((sign^mask)>>b)&1);
		}
		int operator[] (unsigned b) const
		{ return get(b); }

		unsigned hw() const
		{ return hashclash::hw(mask); }
		unsigned hwnaf() const
		{ return hashclash::hwnaf(adddiff()); }
		sdr naf() const
		{ return hashclash::naf(adddiff()); }

		/**** members ****/
		uint32 mask, sign;
	};

	inline void swap(sdr& l, sdr& r)
	{
		std::swap(l.mask, r.mask);
		std::swap(l.sign, r.sign);
	}

	std::ostream& operator<<(std::ostream& o, const sdr& n);
	std::istream& operator>>(std::istream& i, sdr& n);

	extern unsigned hw_table[0x800];
	inline unsigned hw(uint32 n)
	{
		unsigned w = hw_table[n & 0x7FF];
		w += hw_table[(n >> 11) & 0x7FF];
		w += hw_table[n >> 22];
		return w;
	}

	inline unsigned hwnaf(uint32 n)
	{
		uint32 a = n>>1;
		uint32 w = a ^ (n+a);
		return hw(w);
	}

	inline sdr naf(uint32 n)
	{
		uint32 a = n>>1;
		return sdr(a,a+n);
	}

	inline unsigned hw(const sdr& n)
	{ return n.hw(); }
	inline unsigned hwnaf(const sdr& n)
	{ return n.hwnaf(); }
	inline sdr naf(const sdr& n)
	{ return n.naf(); }

	unsigned count_sdrs(uint32 n, unsigned maxweight = 32);
	void table_sdrs(std::vector<sdr>& result, uint32 n, unsigned maxweight);

	unsigned count_sdrs(uint32 n, unsigned weight, bool signpos);
	void table_sdrs(std::vector<sdr>& result, uint32 n, unsigned weight, bool signpos);

	unsigned count_sdrs(sdr n, unsigned maxweight, unsigned rot);
	void table_sdrs(std::vector<sdr>& result, sdr n, unsigned maxweight, unsigned rot);

	// call this in a other global constructor using hw or hwnaf 
	// due to the inpredictable order of which global constructors are called
	// otherwise the correct functioning of hw and hwnaf are not guaranteed
	void hashclash_sdr_hpp_init();

} // namespace hashclash


#ifndef NOSERIALIZATION
namespace boost {
	namespace serialization {

		template<class Archive>
		void serialize(Archive& ar, hashclash::sdr& d, const unsigned int file_version)
		{
			ar & make_nvp("mask", d.mask);
			ar & make_nvp("sign", d.sign);
		}

	}
}
#endif // NOSERIALIZATION

#endif // HASHCLASH_SIGNED_DIGIT_REPRESENTATION_HPP
