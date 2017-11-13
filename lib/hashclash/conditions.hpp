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

#ifndef HASHCLASH_CONDITIONS_HPP
#define HASHCLASH_CONDITIONS_HPP


#include <iostream>

#ifndef NOSERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // NOSERIALIZATION

#include "types.hpp"
#include "sdr.hpp"

namespace hashclash {

	enum bitcondition {
		bc_constant, bc_plus, bc_minus, bc_zero, bc_one,
		bc_prev, bc_prevn, bc_prev2, bc_prev2n, bc_or2,
		bc_next, bc_nextn, bc_next2, bc_next2n, bc_or2b, // counterparts of prev,...,or2
		bc_max
	};
	inline bitcondition fromdiff(int d) { 
		if (d == -1) return bc_minus;
		if (d == +1) return bc_plus;
		return bc_constant;
	}
			
	std::ostream& operator<<(std::ostream& o, const bitcondition& v);
	std::istream& operator>>(std::istream& i, bitcondition& v);

	inline bitcondition diffbitcondition(int diff)
	{
		if (diff == 0) return bc_constant;
		if (diff == -1) return bc_minus;
		return bc_plus;
	}

	inline bool isdirect(bitcondition cond) {
		if (cond <= bc_one) return true;
		return false;
	}
	inline bool isforward(bitcondition cond) {
		if (cond >= bc_prev && cond <= bc_or2)
			return false;
		return true;
	}
	inline bool isbackward(bitcondition cond) {
		if (cond >= bc_next)
			return false;
		return true;
	}
	inline bool isindirect1(bitcondition cond) {
		if (cond == bc_next || cond == bc_nextn
			|| cond == bc_prev || cond == bc_prevn)
			return true;
		return false;
	}
	inline bool isindirect2(bitcondition cond) {
		if ((cond >= bc_prev2 && cond <= bc_or2)
			|| cond >= bc_next2)
			return true;
		return false;
	}

	struct byteconditions {
		uint32 val;

		byteconditions(): val(0) {}
		byteconditions(uint32 value): val(value) {}
		byteconditions(const byteconditions& r): val(r.val) {}
		byteconditions(bitcondition b0, bitcondition b1 = bc_constant
			, bitcondition b2 = bc_constant, bitcondition b3 = bc_constant
			, bitcondition b4 = bc_constant, bitcondition b5 = bc_constant
			, bitcondition b6 = bc_constant, bitcondition b7 = bc_constant)
			: val(b0 + (b1<<4) + (b2<<8) + (b3<<12) + (b4<<16) + (b5<<20) + (b6<<24) + (b7<<28))
		{}

		void clear()
		{ val = 0; }

		byteconditions& set(uint32 value)
		{ val = value; return *this; }
		byteconditions& set(const byteconditions& r)
		{ val = r.val; return *this; }
		byteconditions& set(bitcondition b0, bitcondition b1 = bc_constant
			, bitcondition b2 = bc_constant, bitcondition b3 = bc_constant
			, bitcondition b4 = bc_constant, bitcondition b5 = bc_constant
			, bitcondition b6 = bc_constant, bitcondition b7 = bc_constant)
		{ 
			val = b0 + (b1<<4) + (b2<<8) + (b3<<12) + (b4<<16) + (b5<<20) + (b6<<24) + (b7<<28);
			return *this;
		}
		byteconditions& operator= (const byteconditions& r)
		{ return set(r); }

		byteconditions& set(unsigned b, bitcondition cond)
		{
			b <<= 2;
			val &= ~(0xF << b);
			val += cond << b;
			return *this;
		}
		bitcondition get(unsigned b) const
		{ return bitcondition( (val >> (b<<2)) & 0xF ); }
		bitcondition operator[] (unsigned b) const
		{ return get(b); }

		bool operator== (const byteconditions& r) const
		{ return val==r.val; }
		bool operator< (const byteconditions& r) const
		{ return val<r.val; }
		bool operator!= (const byteconditions& r) const
		{ return !(*this == r); }
		bool operator> (const byteconditions& r) const
		{ return r < *this; }
		bool operator<= (const byteconditions& r) const
		{ return !(*this > r); }
		bool operator>= (const byteconditions& r) const
		{ return !(*this < r); }

		uint32 mask() const
		{
			uint32 m = 0;
			for (unsigned b = 0; b < 8; ++b)
				if (get(b) != bc_constant) m |= 1<<b;
			return m;
		}
		unsigned hw() const
		{ return hashclash::hw(mask()); }
		uint32 mask(bitcondition cond) const
		{
			uint32 m = 0;
			for (unsigned b = 0; b < 8; ++b)
				if (get(b) == cond) m |= 1<<b;
			return m;
		}
		uint32 diff() const
		{ return mask(bc_plus) - mask(bc_minus); }
		uint32 set0() const
		{ return ~(mask(bc_zero) | mask(bc_plus)); }
		uint32 set1() const
		{ return mask(bc_one) | mask(bc_minus); }
		uint32 prev() const
		{ return mask(bc_prev); }
		uint32 prevn() const
		{ return mask(bc_prevn); }
		uint32 prev2() const
		{ return mask(bc_prev2); }
		uint32 prev2n() const
		{ return mask(bc_prev2n); }
		uint32 next() const
		{ return mask(bc_next); }
		uint32 nextn() const
		{ return mask(bc_nextn); }
		uint32 next2() const
		{ return mask(bc_next2); }
		uint32 next2n() const
		{ return mask(bc_next2n); }
		uint32 or2() const
		{ return mask(bc_or2); }
		uint32 or2b() const
		{ return mask(bc_or2b); }
	};
	std::ostream& operator<<(std::ostream& o, const byteconditions& b);
	std::istream& operator>>(std::istream& i, byteconditions& b);

	struct wordconditions {
		byteconditions bytes[4];
		
		wordconditions() {}
		wordconditions(const wordconditions& r)
		{ set(r); }
		wordconditions(const sdr& r)
		{ set(r); }
		wordconditions(uint32 diff, uint32 mset0, uint32 mset1)
		{ set(diff, mset0, mset1); }

		void clear()
		{ bytes[0].clear(); bytes[1].clear(); bytes[2].clear(); bytes[3].clear(); }

		wordconditions& set(const wordconditions& r)
		{ 
			bytes[0] = r.bytes[0]; 
			bytes[1] = r.bytes[1];
			bytes[2] = r.bytes[2];
			bytes[3] = r.bytes[3];
			return *this;
		}
		wordconditions& set(const sdr& r)
		{
			clear();
			for (unsigned b = 0; b < 32; ++b)
				if (r.mask & (1<<b)) {
					if (r.sign & (1<<b))
						bytes[b>>3].set(b&7, bc_plus);
					else
						bytes[b>>3].set(b&7, bc_minus);
				}
			return *this;
		}
		wordconditions& set(uint32 diff, uint32 mset0, uint32 mset1)
		{
			sdr tmp = sdr(mset1, mset1+diff);
			set(tmp);
			mset0 = ~(mset0 | tmp.mask);
			mset1 &= ~tmp.mask;
			for (unsigned b = 0; b < 32; ++b)
			{
				if ((mset0>>b)&1) bytes[b>>3].set(b&7, bc_zero);
				else if ((mset1>>b)&1) bytes[b>>3].set(b&7, bc_one);
			}
			return *this;
		}
		wordconditions& operator= (const wordconditions& r)
		{ return set(r); }
		wordconditions& operator= (const sdr& r)
		{ return set(r); }

		bool operator== (const wordconditions& r) const
		{ return bytes[0]==r.bytes[0] && bytes[1]==r.bytes[1] && bytes[2]==r.bytes[2] && bytes[3]==r.bytes[3]; }
		bool operator< (const wordconditions& r) const
		{
			for (unsigned k = 0; k < 3; ++k)
			{
				if (bytes[k] < r.bytes[k]) return true;
				if (bytes[k] > r.bytes[k]) return false;
			}
			if (bytes[3] < r.bytes[3]) return true;
			return false;
		}
		bool operator!= (const wordconditions& r) const
		{ return !(*this == r); }
		bool operator> (const wordconditions& r) const
		{ return r < *this; }
		bool operator<= (const wordconditions& r) const
		{ return !(*this > r); }
		bool operator>= (const wordconditions& r) const
		{ return !(*this < r); }

		uint32 mask() const
		{ return bytes[0].mask() | (bytes[1].mask()<<8) | (bytes[2].mask()<<16) | (bytes[3].mask()<<24); }
		unsigned hw() const
		{ return bytes[0].hw() + bytes[1].hw() + bytes[2].hw() + bytes[3].hw(); }
		uint32 diff() const
		{ return bytes[0].diff() + (bytes[1].diff()<<8) + (bytes[2].diff()<<16) + (bytes[3].diff()<<24); }
		sdr getsdr() const
		{
			uint32 mset1 = set1();
			return sdr(mset1, mset1 + diff());
		}

		uint32 set0() const
		{ return ~((~bytes[0].set0()) | ((~bytes[1].set0())<<8) | ((~bytes[2].set0())<<16) | ((~bytes[3].set0())<<24)); }
		uint32 set1() const
		{ return bytes[0].set1() | (bytes[1].set1()<<8) | (bytes[2].set1()<<16) | (bytes[3].set1()<<24); }
		uint32 prev() const
		{ return bytes[0].prev() | (bytes[1].prev()<<8) | (bytes[2].prev()<<16) | (bytes[3].prev()<<24); }
		uint32 prevn() const
		{ return bytes[0].prevn() | (bytes[1].prevn()<<8) | (bytes[2].prevn()<<16) | (bytes[3].prevn()<<24); }
		uint32 prev2() const
		{ return bytes[0].prev2() | (bytes[1].prev2()<<8) | (bytes[2].prev2()<<16) | (bytes[3].prev2()<<24); }
		uint32 prev2n() const
		{ return bytes[0].prev2n() | (bytes[1].prev2n()<<8) | (bytes[2].prev2n()<<16) | (bytes[3].prev2n()<<24); }
		uint32 next() const
		{ return bytes[0].next() | (bytes[1].next()<<8) | (bytes[2].next()<<16) | (bytes[3].next()<<24); }
		uint32 nextn() const
		{ return bytes[0].nextn() | (bytes[1].nextn()<<8) | (bytes[2].nextn()<<16) | (bytes[3].nextn()<<24); }
		uint32 next2() const
		{ return bytes[0].next2() | (bytes[1].next2()<<8) | (bytes[2].next2()<<16) | (bytes[3].next2()<<24); }
		uint32 next2n() const
		{ return bytes[0].next2n() | (bytes[1].next2n()<<8) | (bytes[2].next2n()<<16) | (bytes[3].next2n()<<24); }
		uint32 or2() const
		{ return bytes[0].or2() | (bytes[1].or2()<<8) | (bytes[2].or2()<<16) | (bytes[3].or2()<<24); }
		uint32 or2b() const
		{ return bytes[0].or2b() | (bytes[1].or2b()<<8) | (bytes[2].or2b()<<16) | (bytes[3].or2b()<<24); }

		bitcondition get(unsigned b) const
		{ return bytes[b>>3].get(b&7); }
		bitcondition operator[] (unsigned b) const
		{ return get(b); }
		wordconditions& set(unsigned b, bitcondition cond)
		{
			bytes[b>>3].set(b&7, cond);
			return *this;
		}
	};
	std::ostream& operator<<(std::ostream& o, const wordconditions& w);
	std::istream& operator>>(std::istream& i, wordconditions& w);

} // namespace hashclash

#ifndef NOSERIALIZATION
namespace boost {
	namespace serialization {

		template<class Archive>
		void serialize(Archive& ar, hashclash::byteconditions& b, const unsigned int file_version)
		{
			ar & make_nvp("val", b.val);
		}

		template<class Archive>
		void serialize(Archive& ar, hashclash::wordconditions& w, const unsigned int file_version)
		{
			ar & make_nvp("byte0", w.bytes[0]);
			ar & make_nvp("byte1", w.bytes[1]);
			ar & make_nvp("byte2", w.bytes[2]);
			ar & make_nvp("byte3", w.bytes[3]);
		}

	}
}
#endif // NOSERIALIZATION

#endif //HASHCLASH_CONDITIONS_HPP
