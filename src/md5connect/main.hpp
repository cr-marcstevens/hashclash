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
#include <stdexcept>

#include <boost/filesystem/operations.hpp>
#include <boost/thread.hpp>

#include <hashclash/saveload_gz.hpp>
#include <hashclash/sdr.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/timer.hpp>

using namespace hashclash;
using namespace std;

extern boost::mutex mut;
extern std::string workdir;
class path_container;
void dostep(path_container& container);
bool check_path_collfind(const differentialpath& diffpath, const uint32 mdiff[16]);

struct connect_bitdata {
	uint32 dQt;
	uint32 dQtp1;
	uint32 dFt;
	uint32 dFtp1;
	uint32 dFtp2;
	uint32 dFtp3;

	inline bool operator== (const connect_bitdata& r) const {
		return dQt==r.dQt && dQtp1==r.dQtp1 && dFt==r.dFt && dFtp1==r.dFtp1 && dFtp2==r.dFtp2 && dFtp3==r.dFtp3;
	}
	inline bool operator!= (const connect_bitdata& r) const {
		return !(*this == r);
	}
	inline bool operator< (const connect_bitdata& r) const {
		if (dQt < r.dQt) return true;
		if (dQt > r.dQt) return false;
		if (dQtp1 < r.dQtp1) return true;
		if (dQtp1 > r.dQtp1) return false;
		if (dFt < r.dFt) return true;
		if (dFt > r.dFt) return false;
		if (dFtp1 < r.dFtp1) return true;
		if (dFtp1 > r.dFtp1) return false;
		if (dFtp2 < r.dFtp2) return true;
		if (dFtp2 > r.dFtp2) return false;
		if (dFtp3 < r.dFtp3) return true;
		return false;
	}
};

struct md5_connect_thread {
	md5_connect_thread(): sw(true), countb(33,0), countbaborted(33,0), countbdepth(33,0)
	{}
	void md5_connect(const vector<differentialpath>& lowerpaths, const differentialpath& upperpath, path_container& container);
	timer sw/*(true)*/;
	vector<unsigned> countb/*(33, 0)*/;
	vector<unsigned> countbaborted/*(33, 0)*/;
	vector<unsigned> countbdepth/*(33, 0)*/;
	uint64 count, countall;
	vector<unsigned char> isgood;
	vector<int> lowerpathsmaxtunnel;

	unsigned md5_connect_bits(const vector<differentialpath>& lowers, unsigned index, const differentialpath& upper, path_container& container);
	connect_bitdata startbit0;
	vector<connect_bitdata> bitdataresults[33];
	vector<connect_bitdata> bitdatastart[32];
	vector<connect_bitdata> bitdataend[32];
	vector<byteconditions>  bitdatanewcond[32];
	differentialpath newpath2, tmppath;
	unsigned bindex[32];

	void connectbits2(const connect_bitdata& in, unsigned b, const differentialpath& lower, const differentialpath& upper);
	bf_outcome bfo0, bfo1, bfo2, bfo3;
	bf_conditions bfc0, bfc1, bfc2, bfc3;

	void connectbits(const connect_bitdata& in, vector<connect_bitdata>& out, unsigned b, 
		const differentialpath& lower, const differentialpath& upper, vector<byteconditions>* newconds = 0);
	connect_bitdata result;
	//bf_outcome bfo0, bfo1, bfo2, bfo3;
	//bf_conditions bfc0, bfc1, bfc2, bfc3;
	bitcondition Qt, Qtp1;
	byteconditions newcond;
	bool lastdFp1, lastdFp2, lastdFp3;

	unsigned t;
	booleanfunction* Ft;
	booleanfunction* Ftp1;
	booleanfunction* Ftp2;
	booleanfunction* Ftp3;
	vector<uint32> dFt, dFtp1, dFtp2, dFtp3;
	uint32 dmt, dmtp1, dmtp2, dmtp3;
	uint32 dQtp1, dQtp2, dQtp3, dQtp4;
	differentialpath newpath;

	inline bool isinrange(const differentialpath& needle, const differentialpath& haystack, unsigned b)
	{
		for (unsigned j = t-2; j <= t; ++j)
		{
			const wordconditions& nt = needle[j];
			const wordconditions& ht = haystack[j];
			unsigned k = 0;
			while (k+7 <= b)
			{
				if (nt.bytes[k>>3] != ht.bytes[k>>3])
					return false;
				k += 8;
			}
			for (; k <= b; ++k)
				if (nt[k] != ht[k])
					return false;
		}
		return true;
	}
	inline unsigned binary_search_lower_paths(const vector<differentialpath>& lowerpaths,
									   unsigned i, unsigned b)
	{
		const differentialpath& needle = lowerpaths[i];
		unsigned l = i, u = lowerpaths.size()-1;
		while (l < u)
		{
			unsigned j = (l+u+1)>>1;
			if (isinrange(needle, lowerpaths[j], b)) {
				l = j;
			} else {
				u = j-1;
			}
		}
		return l;
	}
};


extern vector<uint32> lowdQt, lowdQtm1, lowdQtm2, lowdQtm3;

class path_container {
public:
	path_container()
		: modn(1), modi(0), inputfilelow(), inputfilehigh()
		, t(0), noverify(false), noenhancepath(false)
		, verified(0), verifiedbad(0), bestpathcond(1<<20)
		, bestmaxtunnel(0), showstats(false), bestmaxcomp(-1000), threads(1)
	{
		for (unsigned k = 0; k < 16; ++k)
			m_diff[k] = 0;
	}

	~path_container() 
	{
		cout << "Best path: totcompl=" << bestmaxcomp << " tottunnel=" << bestmaxtunnel << ", totcond=" << bestpathcond << endl;
		if (!noverify)
			cerr << "Verified: " << verifiedbad << " bad out of " << verified << endl;
	}

	void push_back(const differentialpath& fullpath)
	{
		differentialpath pathback = fullpath;
		//pathback = fullpath;
		try {
		cleanup(pathback);
		} catch (std::exception& e) {
			cerr << "hashclash::cleanup(differentialpath&): unknown exception!:" << endl << e.what() << endl;
			show_path(pathback, m_diff);
		}catch (...) {
			cerr << "hashclash::cleanup(differentialpath&): unknown exception!:" << endl;
			show_path(pathback, m_diff);
		}

//		++verified;
		if (!noverify && !test_path_fast(pathback, m_diff))
		{
			mut.lock();
			++verified;
			++verifiedbad;
			mut.unlock();
			return;
		}

		unsigned tunnel = totaltunnelstrength(pathback);
		int tuncompl = tunnel;
		for (int k = pathback.tbegin(); k < pathback.tend() && k < 64; ++k)
			if (k >= Qcondstart)
				tuncompl -= int(pathback[k].hw());
		if (tuncompl < bestmaxcomp)
			return;
		if (tunnel < bestmaxtunnel)
			return;
		unsigned cond = 0;
		for (int k = pathback.tbegin(); k < pathback.tend(); ++k)
			cond += pathback[k].hw();
		if (!(tuncompl >= bestmaxcomp || tunnel >= bestmaxtunnel || cond <= bestpathcond)) 
			return;
		try {
			if (!noenhancepath) {
				enhancepath(pathback, m_diff);
				tunnel = totaltunnelstrength(pathback);
				tuncompl = tunnel;
				for (int k = pathback.tbegin(); k < pathback.tend() && k < 64; ++k)
					if (k >= Qcondstart)
						tuncompl -= int(pathback[k].hw());
				if (tuncompl < bestmaxcomp)
					return;
				if (tunnel < bestmaxtunnel)
					return;
				cond = 0;
				for (int k = pathback.tbegin(); k < pathback.tend(); ++k)
					cond += pathback[k].hw();
				if (!(tuncompl >= bestmaxcomp || tunnel >= bestmaxtunnel || cond <= bestpathcond))
					return;
			}
		} catch (std::exception&) {
			return;
		} catch (...) {
			return;
		}

		if (!check_path_collfind(pathback, m_diff)) 
			return;

		if (tuncompl == bestmaxcomp && tunnel == bestmaxtunnel && cond == bestpathcond) {
			mut.lock();
			++verified;
			bestpaths.push_back(pathback);
			if (hw(uint32(bestpaths.size()))==1) {
				save_gz(bestpaths, workdir + "/bestpaths", binary_archive);
				cout << "Best paths: " << bestpaths.size() << endl;
			}
			mut.unlock();
			return;
		}
		if (tuncompl < bestmaxcomp) return;
		if (tunnel < bestmaxtunnel) return;
		if (cond > bestpathcond) return;
		mut.lock();
		++verified;
		bestpaths.clear();
		bestpaths.push_back(pathback);

		bestmaxcomp = tuncompl;
		bestpathcond = cond;
		bestpath = pathback;
		bestmaxtunnel = tunnel;
		show_path(pathback, m_diff);
		double p = test_path(pathback, m_diff);
		cout << "Best path: totcompl=" << bestmaxcomp << " tottunnel=" << bestmaxtunnel << ", totcond=" << bestpathcond << ", p=" << p << endl;
		save_gz(pathback, workdir + "/bestpath_t" + boost::lexical_cast<string>(tunnel) + "_c" + boost::lexical_cast<string>(cond), binary_archive);
		save_gz(pathback, workdir + "/bestpath_new", binary_archive);
		save_gz(bestpaths, workdir + "/bestpaths_new", binary_archive);
		try { boost::filesystem::rename(workdir + "/bestpath.bin.gz", workdir + "/bestpath_old.bin.gz"); } catch (...) {}
		try { boost::filesystem::rename(workdir + "/bestpaths.bin.gz", workdir + "/bestpaths_old.bin.gz"); } catch (...) {}
		try { boost::filesystem::rename(workdir + "/bestpath_new.bin.gz", workdir + "/bestpath.bin.gz"); } catch (...) {}
		try { boost::filesystem::rename(workdir + "/bestpaths_new.bin.gz", workdir + "/bestpaths.bin.gz"); } catch (...) {}
		mut.unlock();
	}

	uint32 m_diff[16];

	unsigned t;
	int Qcondstart;
	
	bool showstats;
	bool noenhancepath;
	bool noverify;
	unsigned verified, verifiedbad;

	unsigned modn;
	unsigned modi;
	std::string inputfilelow, inputfilehigh;
	bool showinputpaths;

	differentialpath bestpath;
	vector<differentialpath> bestpaths;
	volatile unsigned bestpathcond;
	volatile unsigned bestmaxtunnel;
	volatile int bestmaxcomp;
	int threads;
};

struct diffpathlower_less
	: public std::binary_function<differentialpath, differentialpath, bool>
{
	bool operator()(const differentialpath& _Left, const differentialpath& _Right) const
	{
		if (_Left.tend() < _Right.tend() || _Left.tbegin() > _Right.tbegin())
			return true;
		if (_Left.tend() > _Right.tend() || _Left.tbegin() < _Right.tbegin())
			return false;
		if (_Left.path.size() < 3)
			throw std::runtime_error("Lower differential paths of insufficient size");
		unsigned t = _Left.tend() - 1;
		uint32 LdQt   = _Left[t].diff();
		uint32 RdQt   = _Right[t].diff();		
		for (unsigned b = 0; b < 32; ++b)
		{
			if ((LdQt & (1<<b)) < (RdQt & (1<<b)))
				return true;
			if ((LdQt & (1<<b)) > (RdQt & (1<<b)))
				return false;
			if (_Left[t-2][b] < _Right[t-2][b])
				return true;
			if (_Left[t-2][b] > _Right[t-2][b])
				return false;
			if (_Left[t-1][b] < _Right[t-1][b])
				return true;
			if (_Left[t-1][b] > _Right[t-1][b])
				return false;
		}
		return false;
	}
};

struct diffpathupper_less
	: public std::binary_function<differentialpath, differentialpath, bool>
{
	bool operator()(const differentialpath& _Left, const differentialpath& _Right) const
	{
		if (_Left.tend() < _Right.tend() || _Left.tbegin() > _Right.tbegin())
			return true;
		if (_Left.tend() > _Right.tend() || _Left.tbegin() < _Right.tbegin())
			return false;
		if (_Left.path.size() < 3)
			throw std::runtime_error("Lower differential paths of insufficient size");
		unsigned t = _Left.tbegin() - 1;
		uint32 LdQtp1   = _Left[t+1].diff();
		uint32 RdQtp1   = _Right[t+1].diff();
		if (_Left[t+1].hw() > _Right[t+1].hw()) return true;
		if (_Left[t+1].hw() < _Right[t+1].hw()) return false;
		for (unsigned b = 0; b < 32; ++b)
		{
			if ((LdQtp1 & (1<<b)) < (RdQtp1 & (1<<b)))
				return true;
			if ((LdQtp1 & (1<<b)) > (RdQtp1 & (1<<b)))
				return false;
			if (_Left[t+2][b] < _Right[t+2][b])
				return true;
			if (_Left[t+2][b] > _Right[t+2][b])
				return false;
			if (_Left[t+3][b] < _Right[t+3][b])
				return true;
			if (_Left[t+3][b] > _Right[t+3][b])
				return false;
		}
		return false;
	}
};

#endif // MAIN_HPP
