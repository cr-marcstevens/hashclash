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

#ifndef HASHCLASH_SHA1DIFFERENTIALPATH_HPP
#define HASHCLASH_SHA1DIFFERENTIALPATH_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <map>
#include <string>

#ifndef NOSERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // NOSERIALIZATION

#include "types.hpp"
#include "sdr.hpp"
#include "conditions.hpp"
#include "booleanfunction.hpp"
#include "sha1detail.hpp"
#include "saveload_bz2.hpp"

namespace hashclash {

	extern booleanfunction SHA1_F1_data, SHA1_F2_data, SHA1_F3_data, SHA1_F4_data, BF_simplify;

	class sha1differentialpath;
	void show_path(const sha1differentialpath& path, std::ostream& o = std::cout);
	bool test_path(const sha1differentialpath& path);
	void cleanup_path(sha1differentialpath& path);

//	inline bitcondition operator^(const bitcondition& l, const unsigned r)
//	{ return bitcondition(unsigned(l)^r); }
//	inline bitcondition& operator^=(bitcondition& l, const unsigned r)
//	{ return l = bitcondition(unsigned(l)^r); }
	// neglects bf_outcome::fminus when fplus is set
	inline bf_outcome msb_bf_outcome(booleanfunction& bf, bitcondition input1, bitcondition input2, bitcondition input3)
	{
		bf_outcome o = bf.outcome(input1, input2, input3);
		o.c = o.c & ~((o.c & bf_outcome::fplus)<<1);
		return o;
	}
	inline bf_outcome msb_bf_outcome(booleanfunction& bf, const bf_conditions& c) {
		return msb_bf_outcome(bf, c.first, c.second, c.third);
	}
	inline bf_conditions msb_bf_forwardconditions(booleanfunction& bf, bitcondition input1, bitcondition input2, bitcondition input3, bitcondition outcome)
	{
		if (outcome == bc_constant) return bf.forwardconditions(input1, input2, input3, bc_constant);
		bf_outcome o = bf.outcome(input1, input2, input3);
		if (!o.constant()) return BF_simplify.forwardconditions(input1, input2, input3, bc_constant); //bf_conditions(input1, input2, input3);
		if (o.plus() && o.minus()) {
			bf_conditions c = bf.forwardconditions(input1, input2, input3, bc_constant);
			if (c.third == input3) {
				std::cout << "[" << input1 << input2 << input3 << "](" << outcome << ")=>[" << c.first << c.second << c.third << "]" << std::flush;
				throw std::runtime_error("msb_bf_forwardconditions: c.third == input3 ?!?");
			}
			if (c.third >= bc_next && c.third <= bc_next2n) {
				c.third = bitcondition(unsigned(c.third)^unsigned(bc_next)^unsigned(bc_nextn)); // bc_next <=> bc_nextn, bc_next2 <=> bc_next2n
				return c;
			}
			throw std::runtime_error("msb_bf_forwardconditions: unpredicted case ?!?");
		} else return bf.forwardconditions(input1, input2, input3, outcome);
	}
	inline bf_conditions msb_bf_forwardconditions(booleanfunction& bf, const bf_conditions& c, bitcondition outcome)
	{ return msb_bf_forwardconditions(bf, c.first, c.second, c.third, outcome); }
	inline bf_conditions msb_bf_backwardconditions(booleanfunction& bf, bitcondition input1, bitcondition input2, bitcondition input3, bitcondition outcome)
	{
		if (outcome == bc_constant) return bf.backwardconditions(input1, input2, input3, bc_constant);
		bf_outcome o = bf.outcome(input1, input2, input3);
		if (!o.constant()) return BF_simplify.backwardconditions(input1, input2, input3, bc_constant); // bf_conditions(input1, input2, input3);
		if (o.plus() && o.minus()) {
			bf_conditions c = bf.backwardconditions(input1, input2, input3, bc_constant);
			if (c.second != input2) {
				if (c.second >= bc_prev && c.second <= bc_prevn) {
					c.second = bitcondition(unsigned(c.second)^unsigned(bc_prev)^unsigned(bc_prevn));
					return c;
				} else {
					std::cout << "[" << input1 << input2 << input3 << "](" << outcome << ")=>[" << c.first << c.second << c.third << "]" << std::flush;
					throw std::runtime_error("msb_bf_backwardconditions: c.second != {input2,^,!}  ?!?");
				}
			} else if (c.first != input1) {
				if (c.first >= bc_prev2 && c.first <= bc_prev2n) {
					c.first = bitcondition(unsigned(c.first)^unsigned(bc_prev2)^unsigned(bc_prev2n));
					return c;
				} else {
					std::cout << "[" << input1 << input2 << input3 << "](" << outcome << ")=>[" << c.first << c.second << c.third << "]" << std::flush;
					throw std::runtime_error("msb_bf_backwardconditions: c.first != {input1,M,#}  ?!?");
				}
			} else {
				std::cout << "[" << input1 << input2 << input3 << "](" << outcome << ")=>[" << c.first << c.second << c.third << "]" << std::flush;
				throw std::runtime_error("msb_bf_backwardconditions: c.second == input2 && c.first == input1 ?!?");
			}
		} else return bf.backwardconditions(input1, input2, input3, outcome);
	}
	inline bf_conditions msb_bf_backwardconditions(booleanfunction& bf, const bf_conditions& c, bitcondition outcome)
	{ return msb_bf_backwardconditions(bf, c.first, c.second, c.third, outcome); }
	

	// checks the probability of achieving the given dQts with specified msgdiffs in sha1path over steps t= tbegin, ..., tend-1
	// path contains msgdiffs and required bitconditions on Qt and Qt'
	double deep_analysis_path(const sha1differentialpath& sha1path, const uint32 dQt[80], unsigned tbegin, unsigned tend);

	// cleans up differential path to backward conditions only
//	void cleanup(differentialpath& path);

	class sha1differentialpath {
	public:
		void clear()
		{
			path.clear();
			me.clear();
		}
		int tbegin() const
		{ return -offset; }
		int tend() const
		{ return int(path.size())-offset; }
		unsigned nrcond() const
		{
			unsigned cond = 0;
			for (int t = tbegin(); t < tend(); ++t)
				cond += path[offset+t].hw();
			return cond;
		}
		void compress()
		{
			// erase trivial wordconditions at the end
			while (path.size() > 0 && path[path.size()-1] == wordconditions()) {
				path.pop_back();
				me.pop_back();
			}
			// erase trivial wordconditions at the begin
			while (path.size() > 0 && path[0] == wordconditions())
			{
				me.erase(me.begin());
				path.erase(path.begin());
				--offset;
			}
		}

		wordconditions& get(int t)		
		{
			if (path.size() == 0) offset = -t;
			if (offset+t < 0)
			{
				unsigned oldsize = unsigned(path.size());
				int toadd = -(offset+t);				
				path.resize(path.size() + unsigned(toadd));
				me.resize(me.size() + unsigned(toadd));
				for (unsigned k = unsigned(path.size())-1; k >= unsigned(toadd); --k) {
					me[k] = me[k-toadd];
					path[k] = path[k-toadd];
				}
				for (unsigned k = 0; k < unsigned(toadd); ++k) {
					path[k].clear();
					me[k] = sdr();
				}
				offset = -t;
			} else if (offset+t >= int(path.size())) {
				path.resize(offset+t+1);
				me.resize(offset+t+1);
			}
			return path[offset+t];			
		}
		sdr& getme(int t)
		{
			if (path.size() == 0 || offset+t < 0 || offset+t >= int(path.size()))
				wordconditions tmp = get(t);
			return me[offset+t];
		}
		const sdr& getme(int t) const
		{
			if (offset+t<0 || offset+t >= int(me.size()))
				throw std::out_of_range("sha1differentialpath::getme(int t) const: t is out of bounds");
			return me[offset+t];
		}

		wordconditions& operator[] (int t)
		{ return get(t); }
		const wordconditions& operator[] (int t) const
		{ 
#if DEBUG
			if (offset+t < 0 || offset+t >= int(path.size())) 
				throw std::out_of_range("sha1differentialpath::operator[] (int t) const: t is out of bounds");
#endif //DEBUG
			return path[offset+t]; 
		}
		bitcondition operator() (int t, unsigned b) const
		{ 
			if (offset+t < 0 || offset+t >= int(path.size()))
				return bc_constant;
			return path[offset+t].get(b);
		}

		void setbitcondition(int t, unsigned b, bitcondition cond)
		{
			get(t).set(b, cond);
		}
		
		void swap(sha1differentialpath& r)
		{
			std::swap(offset, r.offset);
			std::swap(path, r.path);
			std::swap(me, r.me);
		}

	//private:
		int offset;
		std::vector<wordconditions> path;
		std::vector<sdr> me;
	};

	inline bool operator== (const sha1differentialpath& l, const sha1differentialpath& r)
	{
		int tbegin = l.tbegin() < r.tbegin() ? l.tbegin() : r.tbegin();
		for (int t = tbegin; t < l.tend() || t < r.tend(); ++t)
		{
			if (t >= l.tbegin() && t < l.tend()) {
				if (t >= r.tbegin() && t < r.tend()) {
					if (l[t] != r[t] || l.getme(t) != r.getme(t)) 
						return false;
				} else
					if (l[t] != wordconditions()) 
						return false;
			} else
			{
				if (t >= r.tbegin() && t < r.tend())
					if (wordconditions() != r[t]) 
						return false;
			}
		}
		return true;
	}

	inline bool operator< (const sha1differentialpath& l, const sha1differentialpath& r)
	{
		int tbegin = l.tbegin() < r.tbegin() ? l.tbegin() : r.tbegin();
		for (int t = tbegin; t < l.tend() || t < r.tend(); ++t)
		{
			if (t >= l.tbegin() && t < l.tend()) {
				if (t >= r.tbegin() && t < r.tend()) {
					if (l[t] < r[t]) return true;
					if (l[t] > r[t]) return false;
					if (l.getme(t) < r.getme(t)) return true;
					if (l.getme(t) > r.getme(t)) return false;
				} else {
					if (l[t] < wordconditions()) return true;
					if (l[t] > wordconditions()) return false;
				}
			} else
			{
				if (t >= r.tbegin() && t < r.tend()) {
					if (wordconditions() < r[t]) return true;
					if (wordconditions() > r[t]) return false;
				} 
			}
		}
		return false;
	}

} // namespace hash

#ifndef NOSERIALIZATION

namespace hashclash {
	/* specialization types to distinguish between differentialpath & vector<differentialpath> */
	struct saveload_sha1diffpath {
		sha1differentialpath& p;
		uint32 magic;
		saveload_sha1diffpath(const sha1differentialpath& path): p(const_cast<sha1differentialpath&>(path)), magic(0x23580003) {}
		template<class Archive>
		void serialize(Archive& ar, const unsigned int file_version)
		{			
			ar & boost::serialization::make_nvp("magic", magic);
			if (magic != 0x23580003) throw std::runtime_error("serialize(): type is not sha1differentialpath");
			ar & boost::serialization::make_nvp("path", p);
		}
	};
	template<>
	inline void save(const sha1differentialpath& val, archive_type artype, const path& filepath)
	{
		const saveload_sha1diffpath tmp(val);
		save(tmp, artype, filepath);
	}
	template<>
	inline void load(sha1differentialpath& val, archive_type artype, const path& filepath)
	{
		saveload_sha1diffpath tmp(val);
		load(tmp, artype, filepath);
	}

	struct saveload_vecsha1diffpath {
		std::vector<sha1differentialpath>& p;
		uint32 magic;
		saveload_vecsha1diffpath(const std::vector<sha1differentialpath>& path): p(const_cast<std::vector<sha1differentialpath>&>(path)), magic(0x23580004) {}
		template<class Archive>
		void serialize(Archive& ar, const unsigned int file_version)
		{			
			ar & boost::serialization::make_nvp("magic", magic);
			if (magic == 0x23580003) {
				p.resize(1);
				ar & boost::serialization::make_nvp("path", p[0]);
				return;
			}
			if (magic != 0x23580004) throw std::runtime_error("serialize(): type is not vector<sha1differentialpath>");
			ar & boost::serialization::make_nvp("path", p);
		}
	};
	template<>
	inline void save(const std::vector<sha1differentialpath>& val, archive_type artype, const path& filepath)
	{
		const saveload_vecsha1diffpath tmp(val);
		save(tmp, artype, filepath);
	}
	template<>
	inline void load(std::vector<sha1differentialpath>& val, archive_type artype, const path& filepath)
	{
		saveload_vecsha1diffpath tmp(val);
		load(tmp, artype, filepath);
	}
	template<>
	inline void save_bz2(const sha1differentialpath& val, archive_type artype, const path& filepath)
	{
		const saveload_sha1diffpath tmp(val);
		save_bz2(tmp, artype, filepath);
	}
	template<>
	inline void load_bz2(sha1differentialpath& val, archive_type artype, const path& filepath)
	{
		saveload_sha1diffpath tmp(val);
		load_bz2(tmp, artype, filepath);
	}
	template<>
	inline void save_bz2(const std::vector<sha1differentialpath>& val, archive_type artype, const path& filepath)
	{
		const saveload_vecsha1diffpath tmp(val);
		save_bz2(tmp, artype, filepath);
	}
	template<>
	inline void load_bz2(std::vector<sha1differentialpath>& val, archive_type artype, const path& filepath)
	{
		saveload_vecsha1diffpath tmp(val);
		load_bz2(tmp, artype, filepath);
	}
} // namespace hashclash

namespace boost {
	namespace serialization {

		template<class Archive>
		void serialize(Archive& ar, hashclash::sha1differentialpath& p, const unsigned int file_version)
		{			
			ar & make_nvp("offset", p.offset);
			ar & make_nvp("path", p.path);
			ar & make_nvp("me", p.me);
		}

	}
}
#endif // NOSERIALIZATION

#endif //HASHCLASH_SHA1DIFFERENTIALPATH_HPP
