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

#include <hashclash/booleanfunction.hpp>
#include <hashclash/sha1differentialpath.hpp>

using namespace hashclash;
using namespace std;

extern std::string workdir;
class path_container_autobalance;
void dostep(path_container_autobalance& container, bool savetocache = false);
void sha1_forward_differential_step(const hashclash::sha1differentialpath& path
								   , path_container_autobalance& container);



class path_container_autobalance {
public:
	path_container_autobalance()
		: modn(1), modi(0), inputfile(), showinputpaths(false), newinputpath(false)
		, t(0), trange(0), maxweight(1), maxsdrs(1), maxcond(2176), tbegin(-4), ubound(0)
		, estimatefactor(0), includenaf(false)
		, nafestweight(0), halfnafweight(false)
		, pathsout(0), size(0), count(0), count_balanced(0)
		, condcount(0), minweight(0), fillfraction(0)
		, onemessagediff(false), expandprevmessagediff(false)
	{
		for (unsigned k = 0; k < 80; ++k)
			m_mask[k] = 0;
	}

	~path_container_autobalance() 
	{
		cerr << "Autobalance parameters: maxcond=" << maxcond << endl;
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
		condcount[cond] += amount;
		size += amount;
		static unsigned uboundf = unsigned(double(ubound) * estimatefactor);
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

	void push_back(const sha1differentialpath& path, unsigned cond = 0)
	{ 
		if (cond == 0)
		{
			for (int k = tbegin; k <= int(t); ++k)
				if (k >= path.tbegin())
					cond += path[k].hw();
			if (includenaf)
				if (halfnafweight)
					cond += (path[t+1].hw()>>1);
				else
					cond += path[t+1].hw();
		}
		if (cond > maxcond) return;

		for (int k = tunnelconditions.tbegin(); k < tunnelconditions.tend(); ++k) {
		    if (k >= path.tbegin() && k < path.tend())
			for (unsigned b = 0; b < 32; ++b) {
				bitcondition bc = tunnelconditions(k,b), bcpath = path(k,b);
				if (bcpath == bc_constant && bc != bc_prev && bc != bc_prevn) continue;
				if (bcpath == bc_plus) bcpath = bc_zero;
				else if (bcpath == bc_minus) bcpath = bc_one;
				switch (bc) {
					case bc_constant: break;
					case bc_plus: return;
					case bc_minus: return;
					case bc_one:
					case bc_zero:
						if (bcpath == bc) break;
						if (bcpath == bc_one || bcpath == bc_zero) return;
						if (bcpath == bc_plus || bcpath == bc_minus) return;
						break;
					case bc_prev:{
						bitcondition bcpathm1 = path(k-1,b);
						if (bcpathm1 == bc_plus) bcpathm1 = bc_zero;
						else if (bcpathm1 == bc_minus) bcpathm1 = bc_one;
						if (bcpathm1 == bc_next) break;
//						if (bcpath == bc_plus || bcpath == bc_minus || bcpathm1 == bc_plus || bcpathm1 == bc_minus) return;
						if (bcpathm1 == bc_constant) break;
						if (bcpathm1 == bc_one && bcpath == bc_zero) return;
						if (bcpathm1 == bc_zero && bcpath == bc_one) return;
						break;}
					case bc_prevn:{
						bitcondition bcpathm1 = path(k-1,b);
						if (bcpathm1 == bc_plus) bcpathm1 = bc_zero;
						else if (bcpathm1 == bc_minus) bcpathm1 = bc_one;
						if (bcpathm1 == bc_nextn) break;
						if (bcpath == bc_plus || bcpath == bc_minus || bcpathm1 == bc_plus || bcpathm1 == bc_minus) return;
						if (bcpathm1 == bc_constant) break;
						if (bcpathm1 == bc_one && bcpath == bc_one) return;
						if (bcpathm1 == bc_zero && bcpath == bc_zero) return;
						break;}
				}
			}
		}

		if (ubound == 0) {
			pathsout[cond].push_back(path);
			return;
		}
		if (pathsout[cond].size() < ubound) {
			pathsout[cond].push_back(path); 
			++size, ++count;
			autobalance();
		}
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
		if (ubound != 0 && maxcond+1 < pathsout.size() && size + pathsout[maxcond+1].size() > ubound)
			pathsout[maxcond+1].pop_back();
	}

	void export_results(vector< sha1differentialpath >& outpaths)
	{
		outpaths.clear();
		outpaths.reserve(ubound);
		for (unsigned k = 0; k < pathsout.size() && k <= maxcond; ++k)
		{
			unsigned index = unsigned(outpaths.size());
			outpaths.resize(index + pathsout[k].size());
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
			std::copy(pathsout[k].begin(), pathsout[k].end(), outpaths.begin()+index);
#else
			for (std::size_t i = 0; i < pathsout[k].size(); ++i)
				outpaths[index+i] = std::move(pathsout[k][i]);
#endif
			pathsout[k].clear();
		}
		unsigned hbound = unsigned(double(ubound)*fillfraction);
		unsigned k = maxcond+1;
		while (outpaths.size() < hbound && k < pathsout.size()) {
			unsigned index = unsigned(outpaths.size());
			unsigned length = hbound-index;
			if (length > pathsout[k].size())
				length = pathsout[k].size();
			outpaths.resize(index + length);
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
			std::copy(pathsout[k].begin(), pathsout[k].begin()+length, outpaths.begin()+index);
#else
			for (std::size_t i = 0; i < length; ++i)
				outpaths[index+i] = std::move(pathsout[k][i]);
#endif
			pathsout[k].clear();
			++k;
		}
		pathsout.clear();
	}

	uint32 m_mask[80];

	unsigned t, trange;
	int tbegin;
	unsigned maxcond;
	unsigned maxsdrs;
	unsigned maxweight;
	unsigned minweight;
	
	bool includenaf;
	unsigned nafestweight;
	bool halfnafweight;
	
	unsigned modn;
	unsigned modi;
	std::string inputfile;
	bool showinputpaths;
	bool newinputpath;

	double estimatefactor;
	unsigned ubound;
	unsigned size;
	unsigned count, count_balanced;
	vector< unsigned > condcount;
	vector< vector< sha1differentialpath > > pathsout;

	bool onemessagediff;
	bool expandprevmessagediff;

	double fillfraction;
	sha1differentialpath tunnelconditions;
};



#endif // MAIN_HPP
