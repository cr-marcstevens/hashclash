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
#include <cstdio>
#include <set>
#include <algorithm>
#include <cmath>

#include <boost/thread.hpp>

#include "saveload_bz2.hpp"
#include "saveload_gz.hpp"
#include "sdr.hpp"

using namespace std;

namespace hashclash {

	std::ostream& operator<<(std::ostream& o, const sdr& n)
	{
		o << "[!";
		bool first = true;
		for (unsigned b = 0; b < 32; ++b)
		{
			if (n[b] != 0)
			{
				if (first)
					first = false;
				else
					o << ",";
				if (n[b] == -1)
					o << "-";
				o << b;
			}
		}
		o << "!]";
		return o;
	}

	std::istream& operator>>(std::istream& i, sdr& n)
	{
		n.clear();
		char c;
		if (!(i >> c)) return i;
		if (c != '[') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		if (!(i >> c)) return i;
		if (c != '!') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		// check for empty sdr
		if (!(i >> c)) return i;
		if (c == '!') {
			if (!(i >> c)) return i;
			if (c != ']') {
				i.putback(c);
				i.setstate(std::ios::failbit);
			}
			return i;
		}
		i.putback(c);
			
		char s;
		bool neg = true;
		if (!(i >> s)) return i;
		if (s != '-') {
			neg = false;
			i.putback(s);
		}
		unsigned bit;
		while (i >> bit >> c) {
			if (!(bit < 32 && (c == ',' || c == '!'))) {
				i.putback(c);
				i.setstate(std::ios::failbit);
				return i;
			}
			n.mask |= 1<<bit;
			if (!neg)
				n.sign |= 1<<bit;

			if (c == '!') break;
			neg = true;
			if (!(i >> s)) return i;
			if (s != '-') {
				neg = false;
				i.putback(s);
			}
		}
		if (!(i >> c)) return i;
		if (c != ']') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		return i;
	}

	uint32 best_rotated_difference(uint32 diff, int rc)
	{
		if (diff == 0 || (rc&31)==0)
			return diff;

		int rc2 = 32 - rc;
		uint32 bound = 1 << rc2;
		uint32 bound2 = 1 << rc;
		uint32 y = diff >> rc2;
		uint32 x = diff - (y<<rc2);
		uint32 p1 = (bound-x)*y;
		uint32 d1 = ((x<<rc)|y);
		if ((y<<1) > bound2)
			d1 -= bound2;
		else
			p1 = ((bound-x)<<rc) - p1;

		if (x == 0)
			return d1;

		y += 1;
		y &= ~bound2; // erase possible carry at rc-th bit
		uint32 p2 = x*y;
		uint32 d2 = ((x<<rc)|y);
		if ((y<<1) > bound2)
			d2 -= bound2;
		else
			p2 = (x<<rc) - p2;
		if (p1 > p2)
			return d1;
		return d2;
	}

	void rotate_difference(uint32 diff, int rc, std::vector<uint32>& rotateddiff, uint32 minprob)
	{
		rotateddiff.clear();
		if (diff == 0 || (rc&31)==0) {
			rotateddiff.push_back(diff);
			return;
		}
		// now p1,p2,p3,p4 < 1 * 2^32, no overflows

		int rc2 = 32 - rc;
		uint32 bound = 1 << rc2;
		uint32 bound2 = 1 << rc;
		uint32 y = diff >> rc2;
		uint32 x = diff - (y<<rc2);

		uint32 p1 = (bound-x) * y;
		if (p1 >= minprob)
			rotateddiff.push_back( ((x<<rc)|y) - bound2 );
		uint32 p2 = ((bound-x)<<rc) - p1;
		if (p2 >= minprob)
			rotateddiff.push_back( (x<<rc)|y );

		if (x != 0) {
			y += 1; 
			y &= ~bound2; // erase possible carry at rc-th bit

			uint32 p3 = x * y;
			if (p3 >= minprob)
				rotateddiff.push_back( ((x<<rc)|y) - bound2 );
			uint32 p4 = (x<<rc) - p1;
			if (p4 >= minprob)
				rotateddiff.push_back( (x<<rc)|y );
		}
	}

	void rotate_difference(uint32 diff, int rc, std::vector<std::pair<uint32,double> >& rotateddiff)
	{
		double pinv = pow(double(2),-32);

		rotateddiff.clear();
		if (diff == 0 || (rc&31)==0) {
			rotateddiff.push_back(std::pair<uint32,double>(diff,1));
			return;
		}

		int rc2 = 32 - rc;
		uint32 bound = 1 << rc2;
		uint32 bound2 = 1 << rc;
		uint32 y = diff >> rc2;
		uint32 x = diff - (y<<rc2);

		uint32 p1 = 0;
		if (y) {
			p1 = (bound-x) * y;
			rotateddiff.push_back( std::pair<uint32,double>
				( ((x<<rc)|y) - bound2, double(p1)*pinv ));			
		}

		uint32 p2 = ((bound-x)<<rc) - p1;
		rotateddiff.push_back( std::pair<uint32,double>
				( (x<<rc)|y, double(p2)*pinv ));

		uint32 p3 = 0, p4 = 0;
		if (x != 0) {
			y += 1; 
			y &= ~bound2; // erase possible carry at rc-th bit

			if (y) {
				p3 = x * y;
				rotateddiff.push_back( std::pair<uint32,double>
					( ((x<<rc)|y) - bound2, double(p3)*pinv ));			
			}

			p4 = (x<<rc) - p3;
			rotateddiff.push_back( std::pair<uint32,double>
					( (x<<rc)|y, double(p4)*pinv ));
		}
	}

	unsigned hw_table[0x800];
	struct sdr_carry {
		std::vector< std::vector<sdr> > positive, negative;
	};
	std::vector<sdr_carry> hashclash_sc(0);
	std::vector< std::vector< std::pair<unsigned,unsigned> > > hashclash_scn(0);
	void hashclash_init_scn();

	unsigned count_sdrs(uint32 n, unsigned maxw)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		uint32 l = n & 0xFFFF;
		uint32 h1 = n - l;
		uint32 h2 = h1 + 0x10000;
		unsigned h1w = hwnaf(h1);
		unsigned h2w = hwnaf(h2);
		unsigned hw = h1w;
		if (h2w<hw) hw = h2w;

		std::vector< std::pair<unsigned,unsigned> >& vl = hashclash_scn[l];
		std::vector< std::pair<unsigned,unsigned> >& vh1 = hashclash_scn[h1>>16];
		std::vector< std::pair<unsigned,unsigned> >& vh2 = hashclash_scn[h2>>16];
		unsigned tot = 0;

		for (unsigned i = 0; i < vl.size() && i+hw <= maxw; ++i)
		{
			if (vl[i].first == 0 && vl[i].second == 0) continue;
			unsigned m = 0;
			for (unsigned j = 0; j < vh1.size() && i+j<=maxw; ++j)
				m += vh1[j].first + vh1[j].second;
			tot += vl[i].first * m;
			m = 0;
			for (unsigned j = 0; j < vh2.size() && i+j<=maxw; ++j)
				m += vh2[j].first + vh2[j].second;
			tot += vl[i].second * m;
		}
		return tot;
	}

	typedef std::vector< std::vector<sdr> > vec_vec_sdr_t;
	void table_sdrs(std::vector<sdr>& result, uint32 n, unsigned maxw)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		result.clear();
		if (maxw < hwnaf(n)) return;
		std::vector< triple<uint32,uint32,uint32> > breakn(8);
		std::vector< triple<vec_vec_sdr_t*, vec_vec_sdr_t*, vec_vec_sdr_t*> > breaknsdr(8);
		uint32 m0 = n & 0x7FF;
		uint32 n0 = n - m0;
		uint32 m1 = n0 & 0x3FFFFF;
		uint32 m2 = n0 - m1;
		breakn[0]=make_triple(m0,m1,m2);
		breakn[4]=make_triple(m0,m1,m2);
		breaknsdr[0]=make_triple(
			&hashclash_sc[m0].positive,
			&hashclash_sc[m1>>11].positive,
			&hashclash_sc[m2>>21].positive);
		breaknsdr[4]=make_triple(
			&hashclash_sc[m0].positive,
			&hashclash_sc[m1>>11].positive,
			&hashclash_sc[m2>>21].negative);

		m2 += 0x400000;
		breakn[2]=make_triple(m0,m1,m2);
		breakn[6]=make_triple(m0,m1,m2);
		breaknsdr[2]=make_triple(
			&hashclash_sc[m0].positive,
			&hashclash_sc[m1>>11].negative,
			&hashclash_sc[m2>>21].positive);
		breaknsdr[6]=make_triple(
			&hashclash_sc[m0].positive,
			&hashclash_sc[m1>>11].negative,
			&hashclash_sc[m2>>21].negative);

		n0 += 0x800;
		m1 = n0 & 0x3FFFFF;
		m2 = n0 - m1;
		breakn[1]=make_triple(m0,m1,m2);
		breakn[5]=make_triple(m0,m1,m2);
		breaknsdr[1]=make_triple(
			&hashclash_sc[m0].negative,
			&hashclash_sc[m1>>11].positive,
			&hashclash_sc[m2>>21].positive);
		breaknsdr[5]=make_triple(
			&hashclash_sc[m0].negative,
			&hashclash_sc[m1>>11].positive,
			&hashclash_sc[m2>>21].negative);

		m2 += 0x400000;
		breakn[3]=make_triple(m0,m1,m2);
		breakn[7]=make_triple(m0,m1,m2);
		breaknsdr[3]=make_triple(
			&hashclash_sc[m0].negative,
			&hashclash_sc[m1>>11].negative,
			&hashclash_sc[m2>>21].positive);
		breaknsdr[7]=make_triple(
			&hashclash_sc[m0].negative,
			&hashclash_sc[m1>>11].negative,
			&hashclash_sc[m2>>21].negative);
		std::vector<sdr>::const_iterator cit, citend;
		for (unsigned l = 0; l < breakn.size(); ++l)
		{
			m0 = breakn[l].first;
			m1 = breakn[l].second;
			m2 = breakn[l].third;
			vec_vec_sdr_t& v0 = *breaknsdr[l].first;
			vec_vec_sdr_t& v1 = *breaknsdr[l].second;
			vec_vec_sdr_t& v2 = *breaknsdr[l].third;
			unsigned w0 = 0; while (w0 < v0.size() && 0 == v0[w0].size()) ++w0;
			unsigned w1 = 0; while (w1 < v1.size() && 0 == v1[w1].size()) ++w1;
			unsigned w2 = 0; while (w2 < v2.size() && 0 == v2[w2].size()) ++w2;
			unsigned w0max = maxw - w1 - w2;
			unsigned w1max = maxw - w2;
			for (unsigned i0 = w0; i0 <= w0max && i0 < v0.size(); ++i0)
			for (unsigned j0 = 0; j0 < v0[i0].size(); ++j0)
			{
				sdr temp0 = v0[i0][j0];
				for (unsigned i1 = w1; i0+i1 <= w1max && i1 < v1.size(); ++i1)
				for (unsigned j1 = 0; j1 < v1[i1].size(); ++j1)
				{
					sdr temp1 = temp0;
					temp1.mask ^= v1[i1][j1].mask << 11;
					temp1.sign ^= v1[i1][j1].sign << 11;
					for (unsigned i2 = w2; i0+i1+i2 <= maxw && i2 < v2.size(); ++i2)
					{
						cit = v2[i2].begin();
						citend = v2[i2].end();
						for (; cit != citend; ++cit)
						{
							sdr temp2 = temp1;
							temp2.mask ^= cit->mask << 21;
							temp2.sign ^= cit->sign << 21;
							result.push_back(temp2);
						}
					}
				}
			}
		}

	}

	unsigned count_sdrs(uint32 n, unsigned w, bool signpos)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		uint32 l = n & 0xFFFF;
		uint32 h1 = n - l;
		uint32 h2 = h1 + 0x10000;
		unsigned h1w = hwnaf(h1);
		unsigned h2w = hwnaf(h2);
		unsigned hw = h1w;
		if (h2w<hw) hw = h2w;

		std::vector< std::pair<unsigned,unsigned> >& vl = hashclash_scn[l];
		std::vector< std::pair<unsigned,unsigned> >& vh1 = hashclash_scn[h1>>16];
		std::vector< std::pair<unsigned,unsigned> >& vh2 = hashclash_scn[h2>>16];
		unsigned tot = 0;

		for (unsigned i = 0; i < vl.size() && i+hw <= w; ++i)
		{
			if (signpos) {
				if (w-i < vh1.size())
					tot += vl[i].first * vh1[w-i].first;
				if (w-i < vh2.size())
					tot += vl[i].second * vh2[w-i].first;
			} else {
				if (w-i < vh1.size())
					tot += vl[i].first * vh1[w-i].second;
				if (w-i < vh2.size())
					tot += vl[i].second * vh2[w-i].second;
			}
		}
		return tot;
	}

	unsigned count_sdrs(sdr n, unsigned maxw, unsigned rot)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		sdr low = n, high = n;
		high.mask >>= 32-rot; high.sign >>= 32-rot;
		low.mask &= (~uint32(0))>>rot; low.sign &= low.mask;

		uint32 highdiff = high.adddiff();
		uint32 lowdiff = low.adddiff();
		bool lowpos = (lowdiff & 0x80000000)==0;
		bool highpos = (highdiff & 0x80000000)==0;

		unsigned lowcnt[33];
		unsigned highcnt[33];
		unsigned lowhw = hwnaf(lowdiff);
		unsigned highhw = hwnaf(highdiff);
		lowdiff <<= rot;
		highdiff <<= 32-rot;
		for (unsigned i = lowhw; i+highhw <= maxw; ++i)
			lowcnt[i] = count_sdrs(lowdiff, i, lowpos);
		for (unsigned i = highhw; i+lowhw <= maxw; ++i)
			highcnt[i] = count_sdrs(highdiff, i, highpos);
		unsigned tot = 0;
		for (unsigned wl = lowhw; wl+highhw <= maxw; ++wl)
		{
			unsigned m = 0;
			for (unsigned wh = highhw; wh+wl <= maxw; ++wh)
				m += highcnt[wh];
			tot += lowcnt[wl] * m;
		}
		return tot;
	}

	void table_sdrs(std::vector<sdr>& result, uint32 n, unsigned w, bool signpos)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		result.clear();
		if (w < hwnaf(n)) return;
		std::vector< triple<uint32,uint32,uint32> > breakn(4);
		std::vector< triple<vec_vec_sdr_t*, vec_vec_sdr_t*, vec_vec_sdr_t*> > breaknsdr(4);
		uint32 m0 = n & 0x7FF;
		uint32 n0 = n - m0;
		uint32 m1 = n0 & 0x3FFFFF;
		uint32 m2 = n0 - m1;
		breakn[0]=make_triple(m0,m1,m2);
		if (signpos)
			breaknsdr[0]=make_triple(
				&hashclash_sc[m0].positive,
				&hashclash_sc[m1>>11].positive,
				&hashclash_sc[m2>>21].positive);
		else
			breaknsdr[0]=make_triple(
				&hashclash_sc[m0].positive,
				&hashclash_sc[m1>>11].positive,
				&hashclash_sc[m2>>21].negative);

		m2 += 0x400000;
		breakn[2]=make_triple(m0,m1,m2);
		if (signpos)
			breaknsdr[2]=make_triple(
				&hashclash_sc[m0].positive,
				&hashclash_sc[m1>>11].negative,
				&hashclash_sc[m2>>21].positive);
		else
			breaknsdr[2]=make_triple(
				&hashclash_sc[m0].positive,
				&hashclash_sc[m1>>11].negative,
				&hashclash_sc[m2>>21].negative);

		n0 += 0x800;
		m1 = n0 & 0x3FFFFF;
		m2 = n0 - m1;
		breakn[1]=make_triple(m0,m1,m2);
		if (signpos)
			breaknsdr[1]=make_triple(
				&hashclash_sc[m0].negative,
				&hashclash_sc[m1>>11].positive,
				&hashclash_sc[m2>>21].positive);
		else
			breaknsdr[1]=make_triple(
				&hashclash_sc[m0].negative,
				&hashclash_sc[m1>>11].positive,
				&hashclash_sc[m2>>21].negative);

		m2 += 0x400000;
		breakn[3]=make_triple(m0,m1,m2);
		if (signpos)
			breaknsdr[3]=make_triple(
				&hashclash_sc[m0].negative,
				&hashclash_sc[m1>>11].negative,
				&hashclash_sc[m2>>21].positive);
		else
			breaknsdr[3]=make_triple(
				&hashclash_sc[m0].negative,
				&hashclash_sc[m1>>11].negative,
				&hashclash_sc[m2>>21].negative);
		std::vector<sdr>::const_iterator cit, citend;
		for (unsigned l = 0; l < breakn.size(); ++l)
		{
			m0 = breakn[l].first;
			m1 = breakn[l].second;
			m2 = breakn[l].third;
			vec_vec_sdr_t& v0 = *breaknsdr[l].first;
			vec_vec_sdr_t& v1 = *breaknsdr[l].second;
			vec_vec_sdr_t& v2 = *breaknsdr[l].third;
			unsigned w0 = 0; while (w0 < v0.size() && 0 == v0[w0].size()) ++w0;
			unsigned w1 = 0; while (w1 < v1.size() && 0 == v1[w1].size()) ++w1;
			unsigned w2 = 0; while (w2 < v2.size() && 0 == v2[w2].size()) ++w2;
			unsigned w0max = w - w1 - w2;
			unsigned w1max = w - w2;
			for (unsigned i0 = w0; i0 <= w0max && i0 < v0.size(); ++i0)
			for (unsigned j0 = 0; j0 < v0[i0].size(); ++j0)
			{
				sdr temp0 = v0[i0][j0];
				for (unsigned i1 = w1; i0+i1 <= w1max && i1 < v1.size(); ++i1)
				if (w-i0-i1 < v2.size())
				for (unsigned j1 = 0; j1 < v1[i1].size(); ++j1)
				{
					sdr temp1 = temp0;
					temp1.mask ^= v1[i1][j1].mask << 11;
					temp1.sign ^= v1[i1][j1].sign << 11;
					cit = v2[w - i0 - i1].begin();
					citend = v2[w - i0 - i1].end();
					for (; cit != citend; ++cit)
					{
						sdr temp2 = temp1;
						temp2.mask ^= cit->mask << 21;
						temp2.sign ^= cit->sign << 21;
						result.push_back(temp2);
					}
				}
			}
		}
	}

	void table_sdrs(std::vector<sdr>& result, sdr n, unsigned maxw, unsigned rot)
	{
		if (hashclash_scn.size() == 0) hashclash_init_scn();

		result.clear();
		if (maxw < hwnaf(n)) return;
		sdr low = n, high = n;
		high.mask >>= 32-rot; high.sign >>= 32-rot;
		low.mask &= (~uint32(0))>>rot; low.sign &= low.mask;

		uint32 highdiff = high.adddiff();
		uint32 lowdiff = low.adddiff();
		bool lowpos = (lowdiff & 0x80000000)==0;
		bool highpos = (highdiff & 0x80000000)==0;

		std::vector< std::vector<sdr> > lowsdrs(33);
		std::vector< std::vector<sdr> > highsdrs(33);
		unsigned lowhw = hwnaf(lowdiff);
		unsigned highhw = hwnaf(highdiff);
		lowdiff <<= rot;
		highdiff <<= 32-rot;
		for (unsigned i = lowhw; i <= maxw-highhw; ++i)
			table_sdrs(lowsdrs[i], lowdiff, i, lowpos);			
		for (unsigned i = highhw; i <= maxw-lowhw; ++i)
			table_sdrs(highsdrs[i], highdiff, i, highpos);
		for (unsigned wl = lowhw; wl <= maxw-highhw; ++wl)
			for (unsigned l0 = 0; l0 < lowsdrs[wl].size(); ++l0)
			{
				sdr temp0 = lowsdrs[wl][l0];
				if (temp0.hw() != wl) throw;
				temp0.mask >>= rot;
				temp0.sign >>= rot;
				for (unsigned wh = highhw; wh+wl <= maxw; ++wh)
					for (unsigned h0 = 0; h0 < highsdrs[wh].size(); ++h0)
					{
						sdr temp1 = highsdrs[wh][h0];
						if (temp1.hw() != wh) throw;
						temp1.mask ^= temp0.mask;
						temp1.sign ^= temp0.sign;
						result.push_back(temp1);
					}
			}
	}


	boost::mutex hashclash_init_scn_mutex;
	void hashclash_init_scn()
	{
		boost::lock_guard<boost::mutex> lock(hashclash_init_scn_mutex);
		if (hashclash_scn.size()) return;
		std::vector< std::vector< std::pair<unsigned,unsigned> > > hashclash_scn2(0);
/*		try {
			load_gz(hashclash_scn2, "hashclash_scn", binary_archive);
			hashclash_scn.swap(hashclash_scn2);
		} catch(...) {
			hashclash_scn2.clear();
		}
		if (hashclash_scn.size()) return;
*/
		sdr temp(0);
		hashclash_scn2.resize(1<<16);	
		for (temp.mask = 0; temp.mask < 0x10000; ++temp.mask)
			for (temp.sign = 0; temp.sign <= temp.mask; ++temp.sign)
			{
				if (temp.sign & (~temp.mask)) continue;
				unsigned w = temp.hw();
				uint32 n = temp.adddiff();
				if (n & 0x80000000) // negative
				{
					n += 0x10000;
					if (hashclash_scn2[n].size() < w+1)
						hashclash_scn2[n].resize(w+1);
					++hashclash_scn2[n][w].second;
				} else // positive
				{
					if (hashclash_scn2[n].size() < w+1)
						hashclash_scn2[n].resize(w+1);
					++hashclash_scn2[n][w].first;
				}
			}
/*
		try {
			save_gz(hashclash_scn2, "hashclash_scn", binary_archive);
		} catch (...)
		{}
*/
		hashclash_scn.swap(hashclash_scn2);
	}

	struct hashclash_sdr__init {
		hashclash_sdr__init()
		{
			init_hwtable();
			init_sc();
		}

		void init_hwtable()
		{
			for (uint32 n = 0; n < 0x800; ++n)
			{
				unsigned w = 0;
				uint32 k = n;
				while (k) {
					w += k & 1;
					k >>= 1;
				}
				hw_table[n] = w;
			}
		}

		void init_sc()
		{
			hashclash_sc.resize(0x800);
			sdr temp;
			for (temp.mask = 0; temp.mask < 0x800; ++temp.mask)
				for (temp.sign = 0; temp.sign < 0x800; ++temp.sign)
				{
					if (temp.mask != (temp.sign | temp.mask)) continue;
					uint32 n = temp.adddiff();
					unsigned w = temp.hw();
					if (n & 0x80000000) // negative
					{
						n += 0x800;
						if (hashclash_sc[n].negative.size() < w+1)
							hashclash_sc[n].negative.resize(w+1);
						hashclash_sc[n].negative[w].push_back(temp);
					} else // positive
					{
						if (hashclash_sc[n].positive.size() < w+1)
							hashclash_sc[n].positive.resize(w+1);
						hashclash_sc[n].positive[w].push_back(temp);
					}
				}
		}
	};
	hashclash_sdr__init hashclash_sdr__init__now;

	void hashclash_sdr_hpp_init()
	{
		hashclash_sdr__init here;
	}


} // namespace
