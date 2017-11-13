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

#ifndef HASHCLASH_TYPES_HPP
#define HASHCLASH_TYPES_HPP

#include "config.h"

#include <utility>
#include <algorithm>
#include <vector>

#ifdef __CUDACC__
#pragma message("CUDA compiler detected.")
#define NOSERIALIZATION
#define FUNC_PREFIX __device__ __host__
#include <boost/cstdint.hpp>
#else
#define FUNC_PREFIX
#include <boost/cstdint.hpp>
#endif // __CUDACC__

// Since including serialization routines can be intrusive,
// especially if it is not used at all,
// one can disable it by defining NOSERIALIZATION
#ifndef NOSERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // NOSERIALIZATION


namespace hashclash {

	typedef boost::uint32_t uint32;
	typedef boost::uint64_t uint64;
//	typedef boost::uint16_t uint16;
//	typedef boost::uint8_t uint8;
	typedef boost::int8_t int8;

	FUNC_PREFIX inline uint32 rotate_right(const uint32 x, const unsigned n)
	{ return (x>>n) | (x<<(32-n)); }
	FUNC_PREFIX inline uint32 rotate_left(const uint32 x, const unsigned n)
	{ return (x<<n) | (x>>(32-n)); }


	/**** class triple ****/
	template<class Ty1, class Ty2, class Ty3>
	struct triple {
		typedef triple<Ty1, Ty2, Ty3> MyType;
		typedef Ty1 first_type;
		typedef Ty2 second_type;
		typedef Ty3 third_type;

		triple()
			: first(Ty1()), second(Ty2()), third(Ty3()) 
		{}

		triple(const Ty1& val1, const Ty2& val2, const Ty3& val3)
			: first(val1), second(val2), third(val3)
		{}

		template<class O1, class O2, class O3>
		triple(const triple<O1, O2, O3>& r)
			: first(r.first), second(r.second), third(r.third)
		{}

		void swap(MyType& r)
		{
			std::swap(first, r.first);
			std::swap(second, r.second);
			std::swap(third, r.third);
		}

		first_type first;
		second_type second;
		third_type third;
	};

	template<class Ty1, class Ty2, class Ty3>
	inline triple<Ty1,Ty2,Ty3> make_triple(Ty1 v1, Ty2 v2, Ty3 v3)
	{ 
		return triple<Ty1,Ty2,Ty3>(v1,v2,v3); 
	}

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator==(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{
		return l.first==r.first && l.second==r.second && l.third==r.third;
	}

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator!=(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{	return !(l == r); }

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator<(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{
		return l.first<r.first || (l.first==r.first && (l.second<r.second || (l.second==r.second && l.third<r.third)));
	}

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator>(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{	return r<l;	}

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator<=(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{	return !(r<l); }

	template<class Ty1, class Ty2, class Ty3>
	inline bool operator>=(const triple<Ty1, Ty2, Ty3>& l, const triple<Ty1, Ty2, Ty3>& r)
	{	return !(l<r); }

	template<class Ty1, class Ty2, class Ty3>
	inline void swap(triple<Ty1,Ty2,Ty3>& l, triple<Ty1,Ty2,Ty3>& r)
	{	l.swap(r);	}

	template<class Ty>
	struct sortindices {
		uint32 index;
		const Ty* ptr;
		sortindices(): index(~0), ptr(0) {}
		sortindices(uint32 idx, const Ty& val): index(idx), ptr(&val) {}
		bool operator<(const sortindices<Ty>& r) const { return *ptr < *r.ptr; }
		bool operator==(const sortindices<Ty>& r) const { return *ptr == *r.ptr; }
	};

	template<class Ty1, class Ty2>
	inline void sortbyindex(std::vector<Ty2>& tosort, const std::vector< sortindices<Ty1> >& indices)
	{
		std::vector<Ty2> tosorttmp(indices.size());
		for (unsigned i = 0; i < indices.size(); ++i)
			std::swap(tosorttmp[i], tosort[ indices[i].index ]);
		std::swap(tosorttmp, tosort);
	}

	template<class Ty1, class Ty2>
	void friendsort(std::vector<Ty1>& tosort, std::vector<Ty2>& tosortfriend)
	{
		std::vector< sortindices<Ty1> > indices;
		indices.reserve(tosort.size());
		for (unsigned i = 0; i < tosort.size(); ++i)
			indices.push_back( sortindices<Ty1>(i, tosort[i]) );
		std::sort( indices.begin(), indices.end() );
		indices.erase( std::unique(indices.begin(), indices.end()), indices.end());
		sortbyindex(tosort, indices);
		sortbyindex(tosortfriend, indices);		
	}

} // namespace hashclash

#ifndef NOSERIALIZATION
namespace boost {
	namespace serialization {

		template<class Archive, class Ty1, class Ty2, class Ty3>
		void serialize(Archive& ar, hashclash::triple<Ty1,Ty2,Ty3>& t, const unsigned int file_version)
		{
			ar & make_nvp("first", t.first);
			ar & make_nvp("second", t.second);
			ar & make_nvp("third", t.third);
		}

	}
}
#endif // NOSERIALIZATION

#endif // HASHCLASH_TYPES_HPP
