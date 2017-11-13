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

#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <vector>
#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/thread.hpp>

#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>

using namespace hashclash;
using namespace std;

extern boost::mutex mut;
extern std::string workdir;
class path_container_autobalance;
void dostep(path_container_autobalance& container, bool savetocache = false);

struct md5_forward_thread {
	void md5_forward_differential_step(const hashclash::differentialpath& path, path_container_autobalance& container);
	differentialpath newpath;
	vector<sdr> sdrs;
	bitcondition Qtb[32], Qtm1b[32], Qtm2b[32];
	vector<unsigned> bval;
	uint32 fdiv[32];
	bf_outcome foutcomes[32];	
};


class path_container_autobalance {
public:
	path_container_autobalance()
		: modn(1), modi(0), inputfile(), showinputpaths(false)
		, t(0), trange(0), maxweight(1), maxsdrs(1), maxcond(2176), tbegin(-3), ubound(0)
		, estimatefactor(0), noverify(false), normalt01(false), includenaf(false)
		, nafestweight(0), halfnafweight(false), minQ456tunnel(0), minQ91011tunnel(0), minQ314tunnel(0)
		, pathsout(0), size(0), count(0), count_balanced(0)
		, condcount(0), verified(0), verifiedbad(0), minweight(0), fillfraction(0), threads(1)
		, uct(-4), ucb(-1), ucc('.')
	{
		if (uct < -3 || uct > 64 || ucb < 0 || ucb > 32
			|| (ucc != '0' && ucc != '1' && ucc != '^' && ucc != '!' && ucc != 'm' && ucc != '#')
			)
		{
			uct = -4;
			ucb = -1;
			ucc = '.';
		}
		for (unsigned k = 0; k < 4; ++k)
			IHV1[k] = IHV2[k] = 0;
		for (unsigned k = 0; k < 16; ++k)
			m_diff[k] = 0;
	}

	~path_container_autobalance() 
	{
		cerr << "Autobalance parameters: maxcond=" << maxcond << endl;
		if (!noverify)
			cerr << "Verified: " << verifiedbad << " bad out of " << verified << endl;
	}

	void set_parameters()
	{
		pathsout.clear();
		pathsout.resize(maxcond+1);
		condcount.clear();
		condcount.resize(maxcond+1);
	}

	void estimate(unsigned cond, unsigned amount)
	{
		if (cond > maxcond) return;
		if (ubound == 0) return;
		boost::lock_guard<boost::mutex> lock(mut);
		condcount[cond] += amount;
		size += amount;
		unsigned uboundf = unsigned(double(ubound) * estimatefactor);
		if (size > uboundf) {
			unsigned newsize = 0;
			unsigned k = 0;
			while (k < condcount.size() && newsize+condcount[k] <= uboundf)
			{
				newsize += condcount[k];
				++k;
			}
			for (unsigned j = k+1; j <= maxcond && j < condcount.size(); ++j)
				condcount[j] = 0;
			if (newsize != 0) {
				size = newsize;
				maxcond = k-1;
				condcount[k] = 0;
			} else {
				size = newsize + condcount[k];
				maxcond = k;
			}
		}
	}

	void finish_estimate()
	{
		estimatefactor = 0;
		size = 0;
		count = 0;
		count_balanced = 0;
		if (includenaf)
			if (halfnafweight)
				maxcond += (nafestweight>>1);
			else
				maxcond += nafestweight;
		pathsout.clear();
		pathsout.resize(maxcond + 1);
		condcount.clear();
		condcount.resize(maxcond + 1);
	}

	void push_back(const differentialpath& path, unsigned cond = 0)
	{ 
		if (!test_uc(path))
			return;

		if (cond == 0)
		{
			for (int k = max(tbegin,path.tbegin()); k <= int(t); ++k)
				cond += path[k].hw();
			if (includenaf)
				if (halfnafweight)
					cond += (path[t+1].hw()>>1);
				else
					cond += path[t+1].hw();
		}
		if (cond > maxcond) return;
		if (minQ314tunnel > 0 && t >= 3) {
			uint32 q314tunnel = ~path[3].mask() & ~path[2].next();
			if (t >= 4)
				q314tunnel &= ~path[4].mask() | ~path[4].set0();
			if (t >= 5)
				q314tunnel &= ~path[5].mask() | path[5].next() | path[5].set1();
			if (t >= 14)
				q314tunnel &= ~path[14].mask() & ~path[14].next();
			if (t >= 15)
				q314tunnel &= ~path[15].mask() | path[15].next() | ~path[15].set0();
			if (t >= 16)
				q314tunnel &= ~path[16].mask() | ~path[16].set0() | path[16].next();
			if (hw(q314tunnel)<minQ314tunnel) 
				return;			
		}
		if (minQ456tunnel > 0 && t >= 5 && t <= 8) {
			uint32 q4tunnel = ~path[4].mask() & ~path[3].next();
			uint32 q5tunnel = ~path[5].mask() & ~path[4].next();
			if (t == 5) {
				if (hw(q4tunnel|q5tunnel)<minQ456tunnel) return;
			} else {
				uint32 q6free = path[6].next() | ~path[6].mask();
				q4tunnel &= q6free | path[6].set1();
				q5tunnel &= q6free | ~path[6].set0();
				if (hw(q4tunnel|q5tunnel)<minQ456tunnel) 
					return;
			}
		}
		if (minQ91011tunnel > 0 && t >= 10 && t <= 13) {
			uint32 q9tunnel = ~path[9].mask() & ~path[8].next();
			uint32 q10tunnel = ~path[10].mask() & ~path[9].next();
			if (t == 10) {
				if (hw(q9tunnel|q10tunnel)<minQ91011tunnel) return;
			} else {
				uint32 q11free = path[11].next() | ~path[11].mask();
				q9tunnel &= q11free | path[11].set1();
				q10tunnel &= q11free | ~path[11].set0();
				if (hw(q9tunnel|q10tunnel)<minQ91011tunnel) return;
			}
		}
		if (!noverify)
		{
			if (!test_path_fast(path, m_diff))
			{
				mut.lock();
				++verifiedbad;
				++verified;
				mut.unlock();
				return;
			}
		}
		mut.lock();
		if (!noverify) ++verified;
		if (ubound == 0) {
			pathsout[cond].push_back(path);
			mut.unlock();
			return;
		}
		if (pathsout[cond].size() < ubound) {
			pathsout[cond].push_back(path); 

			++size, ++count;
			autobalance();
		}
		mut.unlock();
	}

	void autobalance()
	{
		if (size > ubound) {
			++count_balanced;
			unsigned newsize = 0;
			unsigned k = 0;
			while (k < pathsout.size() && k <= maxcond && newsize+pathsout[k].size() <= ubound)
			{
				newsize += unsigned(pathsout[k].size());
				++k;
			}
			if (newsize != 0) {
				size = newsize;
				maxcond = k-1;
			} else {
				size = newsize + pathsout[k].size();
				maxcond = k;
			}
			for (unsigned j = maxcond+2; j < pathsout.size(); ++j)
				pathsout[j].clear();
		}
	}

	void export_results(vector< differentialpath >& outpaths) const
	{
		outpaths.clear();
		outpaths.reserve(ubound);
		for (unsigned k = 0; k < pathsout.size() && k <= maxcond; ++k)
		{
			unsigned index = unsigned(outpaths.size());
			outpaths.resize(index + pathsout[k].size());
			std::copy(pathsout[k].begin(), pathsout[k].end(), outpaths.begin()+index);
		}
                unsigned hbound = unsigned(double(ubound)*fillfraction);
                unsigned k = maxcond+1;
                while (outpaths.size() < hbound && k < pathsout.size()) {
                        unsigned index = unsigned(outpaths.size());
                        unsigned length = hbound-index;
                        if (length > pathsout[k].size())
                                length = pathsout[k].size();
			outpaths.resize(index + length);
                        std::copy(pathsout[k].begin(), pathsout[k].begin()+length, outpaths.begin()+index);
                        ++k;
                }
	}

	bool test_uc(const differentialpath& path)
	{
		if (ucb != -1 && uct >= path.tbegin() + 1 && uct < path.tend() - 1)
		{
			switch (ucc)
			{
			default: throw std::runtime_error("test_uc: value ucc not allowed");
			case '0':
				return !(path(uct, ucb) == bc_zero || path(uct, ucb) == bc_plus);
			case '1':
				return !(path(uct, ucb) == bc_one || path(uct, ucb) == bc_minus);
			case '^':
				if (uct - 1 < path.tbegin()+1) return true;
				if (path(uct, ucb) == bc_prev) return false;
				if (path(uct - 1, ucb) == bc_next) return false;
				if ((path(uct - 1, ucb) == bc_zero || path(uct - 1, ucb) == bc_plus) && (path(uct, ucb) == bc_zero || path(uct, ucb) == bc_plus)) return false;
				if ((path(uct - 1, ucb) == bc_one || path(uct - 1, ucb) == bc_minus) && (path(uct, ucb) == bc_one || path(uct, ucb) == bc_minus)) return false;
				return true;
			case '!':
				if (uct - 1 < path.tbegin()+1) return true;
				if (path(uct, ucb) == bc_prevn) return false;
				if (path(uct - 1, ucb) == bc_nextn) return false;
				if ((path(uct - 1, ucb) == bc_zero || path(uct - 1, ucb) == bc_plus) && (path(uct, ucb) == bc_one || path(uct, ucb) == bc_minus)) return false;
				if ((path(uct - 1, ucb) == bc_one || path(uct - 1, ucb) == bc_minus) && (path(uct, ucb) == bc_zero || path(uct, ucb) == bc_plus)) return false;
				return true;
			case 'm':
				if (uct - 2 < path.tbegin()+1) return true;
				if (path(uct, ucb) == bc_prev2) return false;
				if (path(uct - 2, ucb) == bc_next2) return false;
				if ((path(uct, ucb) == bc_prev || path(uct - 1, ucb) == bc_next) && (path(uct - 2, ucb) == bc_next || path(uct - 1, ucb) == bc_prev)) return false;
				if ((path(uct - 2, ucb) == bc_zero || path(uct - 2, ucb) == bc_plus) && (path(uct, ucb) == bc_zero || path(uct, ucb) == bc_plus)) return false;
				if ((path(uct - 2, ucb) == bc_one || path(uct - 2, ucb) == bc_minus) && (path(uct, ucb) == bc_one || path(uct, ucb) == bc_minus)) return false;
				return true;
			case '#':
				if (uct - 2 < path.tbegin()+1) return true;
				if (path(uct, ucb) == bc_prev2n) return false;
				if (path(uct - 2, ucb) == bc_next2n) return false;
				if ((path(uct, ucb) == bc_prev || path(uct - 1, ucb) == bc_next) && (path(uct - 2, ucb) == bc_nextn || path(uct - 1, ucb) == bc_prevn)) return false;
				if ((path(uct, ucb) == bc_prevn || path(uct - 1, ucb) == bc_nextn) && (path(uct - 2, ucb) == bc_nextn || path(uct - 1, ucb) == bc_prevn)) return false;
				if ((path(uct - 2, ucb) == bc_zero || path(uct - 2, ucb) == bc_plus) && (path(uct, ucb) == bc_one || path(uct, ucb) == bc_minus)) return false;
				if ((path(uct - 2, ucb) == bc_one || path(uct - 2, ucb) == bc_minus) && (path(uct, ucb) == bc_zero || path(uct, ucb) == bc_plus)) return false;
				return true;
			}
		}
		return true;
	}

	uint32 IHV1[4];
	uint32 IHV2[4];
	uint32 m_diff[16];

	unsigned t, trange;
	int tbegin;
	unsigned maxcond;
	unsigned maxsdrs;
	unsigned maxweight;
	unsigned minweight;
	
	bool includenaf;
	unsigned nafestweight;
	bool halfnafweight;
	
	unsigned minQ456tunnel;
	unsigned minQ91011tunnel;
	unsigned minQ314tunnel;
	
	bool noverify;
	unsigned verified, verifiedbad;

	unsigned modn;
	unsigned modi;
	std::string inputfile;
	bool showinputpaths;
	bool normalt01;

	double estimatefactor;
	unsigned ubound;
	unsigned size;
	unsigned count, count_balanced;
	vector< unsigned > condcount;
	vector< vector< differentialpath > > pathsout;

	double fillfraction;

	int threads;

	int uct, ucb;
	char ucc;

};



#endif // MAIN_HPP
