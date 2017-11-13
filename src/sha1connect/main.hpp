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
#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sdr.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/timer.hpp>

using namespace hashclash;
using namespace std;

extern boost::mutex mut;
extern std::string workdir;
class path_container;
void dostep(path_container& container);

struct connect_bitdata {
	bitcondition rqtm2[2], rqtm1[2], rqt[2], rqtp1[2];
	bitcondition fqtm1[2], fqt[2], fqtp1[2], fqtp2[2];
	uint32 dQtm1;
	uint32 dQt;
	uint32 dQtp1;
	uint32 dFt;
	uint32 dFtp1;
	uint32 dFtp2;
	uint32 dFtp3;
	uint32 dFtp4;

	inline bool operator== (const connect_bitdata& r) const {
		for (unsigned i = 0; i < 2; ++i)
			if (rqtm2[i] != r.rqtm2[i] || rqtm1[i] != r.rqtm1[i] || rqt[i] != r.rqt[i] || rqtp1[i] != r.rqtp1[i]
				|| fqtm1[i] != r.fqtm1[i] || fqt[i] != r.fqt[i] || fqtp1[i] != r.fqtp1[i] || fqtp2[i] != r.fqtp2[i]) 
					return false;
		return dQtm1==r.dQtm1 && dQt==r.dQt && dQtp1==r.dQtp1
			&& dFt==r.dFt && dFtp1==r.dFtp1 && dFtp2==r.dFtp2 && dFtp3==r.dFtp3 && dFtp4==r.dFtp4
			;
	}
	inline bool operator!= (const connect_bitdata& r) const {
		return !(*this == r);
	}
	inline bool operator< (const connect_bitdata& r) const {
		if (dQtm1 < r.dQtm1) return true;
		if (dQtm1 > r.dQtm1) return false;
		if (dFt < r.dFt) return true;
		if (dFt > r.dFt) return false;
		if (dQt < r.dQt) return true;
		if (dQt > r.dQt) return false;
		if (dFtp1 < r.dFtp1) return true;
		if (dFtp1 > r.dFtp1) return false;
		if (dQtp1 < r.dQtp1) return true;
		if (dQtp1 > r.dQtp1) return false;
		if (dFtp2 < r.dFtp2) return true;
		if (dFtp2 > r.dFtp2) return false;
		if (dFtp3 < r.dFtp3) return true;
		if (dFtp3 > r.dFtp3) return false;
		if (dFtp4 < r.dFtp4) return true;
		if (dFtp4 > r.dFtp4) return false;
		for (unsigned i = 0; i < 2; ++i) {
			if (rqtm2[i] < r.rqtm2[i]) return true;
			if (rqtm2[i] > r.rqtm2[i]) return false;
			if (rqtm1[i] < r.rqtm1[i]) return true;
			if (rqtm1[i] > r.rqtm1[i]) return false;
			if (rqt[i] < r.rqt[i]) return true;
			if (rqt[i] > r.rqt[i]) return false;
			if (rqtp1[i] < r.rqtp1[i]) return true;
			if (rqtp1[i] > r.rqtp1[i]) return false;
			if (fqtm1[i] < r.fqtm1[i]) return true;
			if (fqtm1[i] > r.fqtm1[i]) return false;
			if (fqt[i] < r.fqt[i]) return true;
			if (fqt[i] > r.fqt[i]) return false;
			if (fqtp1[i] < r.fqtp1[i]) return true;
			if (fqtp1[i] > r.fqtp1[i]) return false;
			if (fqtp2[i] < r.fqtp2[i]) return true;
			if (fqtp2[i] > r.fqtp2[i]) return false;
		}
		return false;
	}
};

struct sha1_connect_thread {
	sha1_connect_thread() {
		testokcnt = 0; testbadcnt = 0;
		sw.start(); sw2.start();
		bcnt.assign(41,0);
	}
	void sha1_connect(const sha1differentialpath& lowerpath, const vector<sha1differentialpath>& upperpaths, path_container& container, unsigned prev_eq_b);
	vector<uint64> bcnt;
	timer sw, sw2;
	void sha1_connect_clear_uppercache() {
		dpFt.clear();
		dpFtp1.clear();
		dpFtp2.clear();
		dpFtp3.clear();
		dpFtp4.clear();
		dQtp1.clear();
	}
	unsigned int lowpathindex;
	int t;
	uint32 mmaskt, mmasktp1, mmasktp2, mmasktp3, mmasktp4;
	booleanfunction* Ft;
	booleanfunction* Ftp1;
	booleanfunction* Ftp2;
	booleanfunction* Ftp3;
	booleanfunction* Ftp4;
	vector<uint32> dFt, dFtp1, dFtp2, dFtp3, dFtp4;
	sha1differentialpath newpath;

	vector<uint32> dpFt, dpFtp1, dpFtp2, dpFtp3, dpFtp4, dQtp1;
	bitcondition Qtm1b31not, Qtb31not, Qtp1b31not;
	bitcondition Qtm1b31not2, Qtb31not2, Qtp1b31not2;
	uint32 dQtm1b2b31, dQtm1b27b31, dQtb2b31, dQtb27b31, dQtp1b2b31, dQtp1b27b31;

	uint32 Qtm3freemask, Qtm2freemask, Qtm1freemask, Qtfreemask, Qtp1freemask, Qtp2freemask, Qtp3freemask;
	uint64 testokcnt, testbadcnt;
	sha1differentialpath tbc;
	vector<connect_bitdata> connectbits_vec, connectbits_vec2;
	bool usedFtp1, usedFtp2, usedFtp3, usedFtp4;
	vector<unsigned> prevb;

	template<int steps, bool storefdata, bool storenewconds, bool compmincond = false>
	void connectbits_01234(const pair<connect_bitdata,unsigned>& in, map<connect_bitdata,unsigned>& out, unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& dataend, vector<unsigned>& mincond);
	template<bool storenewconds>
	void connectbits_1234(const pair<connect_bitdata,unsigned>& in, map<connect_bitdata,unsigned>& out, unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& dataend, vector<unsigned>& mincond);
	template<bool storenewconds>
	void connectbits_234(const pair<connect_bitdata,unsigned>& in, map<connect_bitdata,unsigned>& out, unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& dataend, vector<unsigned>& mincond);
	template<bool storenewconds>
	void connectbits_34(const pair<connect_bitdata,unsigned>& in, map<connect_bitdata,unsigned>& out, unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& dataend, vector<unsigned>& mincond);
	template<bool storenewconds>
	void connectbits_4(const pair<connect_bitdata,unsigned>& in, map<connect_bitdata,unsigned>& out, unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& dataend, vector<unsigned>& mincond);
	 connect_bitdata result;
	 bf_outcome bfo0, bfo1, bfo2, bfo3, bfo4;
	 bf_conditions bfc0, bfc1, bfc2, bfc3, bfc4;
	 bitcondition Qtm3b, Qtm2b, Qtm1b, Qtb, Qtp1b, Qtp2b, Qtp3b;
	 uint32 dft, dftp1, dftp2, dftp3, dftp4;
	 pair<byteconditions,byteconditions> newcond;
	 bitcondition Ttm3b, Ttm2b, Ttm1b, Ttb, Ttp1b, Ttp2b, Ttp3b;


	unsigned connect_paths(const sha1differentialpath& lowerpath, const sha1differentialpath& upperpath, connect_bitdata& startdata, path_container& container, bool returnb40);
	map<connect_bitdata, unsigned> bitdataresults[41];
	vector<connect_bitdata> bitdataresults_vec[41];
	vector<connect_bitdata> bitdatastart[40];
	vector<connect_bitdata> bitdataend[40];
	vector< unsigned > bitdatamincond[40];
	vector< pair<byteconditions,byteconditions> >  bitdatanewcond[40];
	vector< unsigned > bitdatanewcondhw[40];
	vector< unsigned > mincond[41];	

};

extern uint64 pbcount;// = 0;
extern int lastlowpathindex;// = -1;
extern vector< vector<unsigned> > badcounts;//(85, vector<unsigned>(32,0));

class path_container {
public:
	path_container()
		: modn(1), modi(0), inputfilelow(), inputfilehigh()
		, t(0), bestpathcond(1<<20), lastsavecnt(0), determinelowestcond(false), splitmode(0), threads(1), showstats(false)
		, keepall(false)
	{
		for (unsigned k = 0; k < 80; ++k)
			m_mask[k] = 0;
		shownbadpath = false;
	}

	~path_container() 
	{
		cout << "Best path: totcond=" << bestpathcond << " count=" << bestpaths.size() << endl;
		save_bestpaths();
	}

	void push_back(const sha1differentialpath& fullpath)
	{
		boost::lock_guard<boost::mutex> lock(mut);
		if (hw(++pbcount)==1) cout << "{" << pbcount << "}" << flush;
		unsigned cond = fullpath.nrcond();
		if (cond > bestpathcond) return;
		if (!test_path(fullpath)) {
			if (!shownbadpath) {
				shownbadpath = true;
				cout << "bad path " /*<< lowpathindex*/ << endl;
				show_path(fullpath);
			}
			return;
		}
//		for (int k = fullpath.tbegin(); k < fullpath.tend(); ++k)
//			cond += fullpath[k].hw();
		bool show = false, bad = false;
                for (int k = tunnelconditions.tbegin(); k < tunnelconditions.tend() && !bad; ++k) {
                    if (k >= fullpath.tbegin() && k < fullpath.tend())
                        for (unsigned b = 0; b < 32 && !bad; ++b) {   
                                bitcondition bc = tunnelconditions(k,b), bcpath = fullpath(k,b);
                                if (bcpath == bc_constant && bc != bc_prev && bc != bc_prevn) continue;
                                if (bcpath == bc_plus) bcpath = bc_zero;
                                else if (bcpath == bc_minus) bcpath = bc_one;
                                switch (bc) {
                                        case bc_constant: break;
                                        case bc_plus: { bad = true; ++badcounts[4+k][b]; break; }
                                        case bc_minus: { bad = true; ++badcounts[4+k][b]; break; } 
                                        case bc_one:
                                        case bc_zero:
                                                if (bcpath == bc) break;
                                                if (bcpath == bc_one || bcpath == bc_zero) { bad = true; ++badcounts[4+k][b]; break; }
                                                if (bcpath == bc_plus || bcpath == bc_minus) { bad = true; ++badcounts[4+k][b]; break; }
                                                break;
                                        case bc_prev:{
                                                bitcondition bcpathm1 = fullpath(k-1,b);
                                                if (bcpathm1 == bc_plus) bcpathm1 = bc_zero;
                                                else if (bcpathm1 == bc_minus) bcpathm1 = bc_one;
                                                if (bcpathm1 == bc_next) break;
						if (bcpath == bc_prev) break;
						if (bcpathm1 == bc_nextn) { bad = true; ++badcounts[4+k][b]; break; }
						if (bcpath == bc_prevn) { bad = true; ++badcounts[4+k][b]; break; }
//                                                if (bcpath == bc_plus || bcpath == bc_minus || bcpathm1 == bc_plus || bcpathm1 == bc_minus) return;
                                                if (bcpathm1 == bc_constant) break;
                                                if (bcpathm1 == bc_one && bcpath == bc_zero) { bad = true; ++badcounts[4+k][b]; break; }
                                                if (bcpathm1 == bc_zero && bcpath == bc_one) { bad = true; ++badcounts[4+k][b]; break; }
                                                break;}
                                        case bc_prevn:{
                                                bitcondition bcpathm1 = fullpath(k-1,b);
                                                if (bcpathm1 == bc_plus) bcpathm1 = bc_zero;
                                                else if (bcpathm1 == bc_minus) bcpathm1 = bc_one;
                                                if (bcpathm1 == bc_nextn) break;
                                                if (bcpath == bc_prevn) break;
                                                if (bcpathm1 == bc_next) { bad = true; ++badcounts[4+k][b]; break; }
                                                if (bcpath == bc_prev) { bad = true; ++badcounts[4+k][b]; break; }
  //                                              if (bcpath == bc_plus || bcpath == bc_minus || bcpathm1 == bc_plus || bcpathm1 == bc_minus) return;
                                                if (bcpathm1 == bc_constant) break;
                                                if (bcpathm1 == bc_one && bcpath == bc_one) { bad = true; ++badcounts[4+k][b]; break; }
                                                if (bcpathm1 == bc_zero && bcpath == bc_zero) { bad = true; ++badcounts[4+k][b]; break; }
                                                break;}
                                }
                                if (bad && hw(uint32(badcounts[4+k][b])) == 1) show = true;
                        }
                }
                if (show) {
                	for (unsigned i = 0; i < badcounts.size(); ++i)
                		for (unsigned b = 0; b < badcounts[i].size(); ++b)
                			if (badcounts[i][b])
                				cout << "badcounts[" << int(i)-4 << "," << b << "] = \t" << badcounts[i][b] << endl;
			cout << endl;
                }
		if (bad) {
			static bool badshown = false;
			if (!badshown) {
				badshown = true;
				show_path(fullpath);
			}
			return;
		}
		
		if (!keepall && (cond < bestpathcond || bestpaths.size() == 0))
		{
			bestpath = fullpath;
			bestpaths.clear();
			bestpaths.push_back(fullpath);
			bestpathcond = cond;
			show_path(fullpath);
			lastsavecnt = 0;
			cout << "Best path: totcond=" << bestpathcond << endl;
			save_bz2(fullpath, workdir + "/bestpath_c" + boost::lexical_cast<string>(bestpathcond), binary_archive);
			save_bz2(fullpath, workdir + "/bestpath", binary_archive);
		} else if (keepall || cond == bestpathcond) {
			bestpaths.push_back(fullpath);
			if (hw(uint32(bestpaths.size()))==1)
				cout << "Best path: totcond=" << bestpathcond << " count=" << bestpaths.size() << endl;
		}
	}

	void save_bestpaths(bool silent = false) {
		mut.lock();
		if (lastsavecnt != bestpaths.size()) {
			save_bz2(bestpaths, workdir + "/bestpaths_c" + boost::lexical_cast<string>(bestpathcond), binary_archive);
			save_bz2(bestpaths, workdir + "/bestpaths", binary_archive);
			lastsavecnt = bestpaths.size();
			if (!silent)
				cout << "Best path: totcond=" << bestpathcond << " count=" << bestpaths.size() << endl;
		}
		mut.unlock();
	}
	bool shownbadpath;
	
	uint32 m_mask[80];

	unsigned t;
	int Qcondstart;
	
	bool showstats;
	bool determinelowestcond;

	unsigned modn;
	unsigned modi;
	std::string inputfilelow, inputfilehigh, inputfileredo;
	bool showinputpaths;

	sha1differentialpath bestpath;
	vector<sha1differentialpath> bestpaths;
	unsigned bestpathcond;
	unsigned bestmaxtunnel;
	int bestmaxcomp;
	uint32 lastsavecnt;
	unsigned splitmode;
	unsigned loworder;

	bool keepall;

	sha1differentialpath tunnelconditions;
	int threads;
};

struct diffpathupper_less
	: public std::binary_function<sha1differentialpath, sha1differentialpath, bool>
{
	bool operator()(const sha1differentialpath& _Left, const sha1differentialpath& _Right) const
	{
		if (_Left.tend() < _Right.tend() || _Left.tbegin() > _Right.tbegin())
			return true;
		if (_Left.tend() > _Right.tend() || _Left.tbegin() < _Right.tbegin())
			return false;
		if (_Left.path.size() < 4)
			throw std::runtime_error("Upper differential paths of insufficient size");
		unsigned t = _Left.tbegin() - 1;
		uint32 LdFt = _Left[t+1].diff();
		uint32 LdFtp1 = _Left[t+2].diff() - _Left[t+1].getsdr().rotate_left(5).adddiff();
		uint32 LdFtp2 = _Left[t+3].diff() - _Left[t+2].getsdr().rotate_left(5).adddiff();
		uint32 LdFtp3 = _Left[t+4].diff() - _Left[t+3].getsdr().rotate_left(5).adddiff();
		uint32 LdFtp4 = _Left[t+5].diff() - _Left[t+4].getsdr().rotate_left(5).adddiff();
		uint32 RdFt = _Right[t+1].diff();
		uint32 RdFtp1 = _Right[t+2].diff() - _Right[t+1].getsdr().rotate_left(5).adddiff();
		uint32 RdFtp2 = _Right[t+3].diff() - _Right[t+2].getsdr().rotate_left(5).adddiff();
		uint32 RdFtp3 = _Right[t+4].diff() - _Right[t+3].getsdr().rotate_left(5).adddiff();
		uint32 RdFtp4 = _Right[t+5].diff() - _Right[t+4].getsdr().rotate_left(5).adddiff();
		uint32 LdQtp1 = LdFt;
		uint32 RdQtp1 = RdFt;
		unsigned b = 0;
		while (b < 32+8) {
			if (b < 32 && (LdFt^RdFt)&(1<<b)) break;
			if (b >= 2 && b-2 < 32 && ((LdFtp1^RdFtp1)&(1<<(b-2)))) break;
			if (b >= 4 && b-4 < 32 && ((LdQtp1^RdQtp1)&(1<<(b-4)))) break;
			if (b >= 4 && b-4 < 32 && ((LdFtp2^RdFtp2)&(1<<(b-4)))) break;
			if (b >= 6 && b-6 < 32 && _Left(t+2,b-6) != _Right(t+2,b-6)) break;
			if (b >= 6 && b-6 < 32 && ((LdFtp3^RdFtp3)&(1<<(b-6)))) break;
			if (b >= 8 && b-8 < 32 && _Left(t+3,b-8) != _Right(t+3,b-8)) break;
			if (b >= 8 && b-8 < 32 && ((LdFtp4^RdFtp4)&(1<<(b-8)))) break;
			++b;
		}
		if (b == 32+8)
			return false;
		if (b < 32) {
			LdFt >>= b; LdFt &= 1; RdFt >>= b; RdFt &= 1;
			if (LdFt < RdFt) return true;
			if (LdFt > RdFt) return false;
		}
		if (b >= 2 && b-2 < 32) {
			LdFtp1 >>= b-2; LdFtp1 &= 1; RdFtp1 >>= b-2; RdFtp1 &= 1;
			if (LdFtp1 < RdFtp1) return true;
			if (LdFtp1 > RdFtp1) return false;
		}
		if (b >= 4 && b-4 < 32) {
			LdQtp1 >>= b-4; LdQtp1 &= 1; RdQtp1 >>= b-4; RdQtp1 &= 1;
			if (LdQtp1 < RdQtp1) return true;
			if (LdQtp1 > RdQtp1) return false;		
			LdFtp2 >>= b-4; LdFtp2 &= 1; RdFtp2 >>= b-4; RdFtp2 &= 1;
			if (LdFtp2 < RdFtp2) return true;
			if (LdFtp2 > RdFtp2) return false;
		}
		if (b >= 6 && b-6 < 32) {
			if (_Left(t+2,b-6) < _Right(t+2,b-6)) return true;
			if (_Left(t+2,b-6) > _Right(t+2,b-6)) return false;		
			LdFtp3 >>= b-6; LdFtp3 &= 1; RdFtp3 >>= b-6; RdFtp3 &= 1;
			if (LdFtp3 < RdFtp3) return true;
			if (LdFtp3 > RdFtp3) return false;
		}
		if (b >= 8 && b-8 < 32) {
			if (_Left(t+2,b-8) < _Right(t+2,b-8)) return true;
			if (_Left(t+2,b-8) > _Right(t+2,b-8)) return false;		
			LdFtp4 >>= b-8; LdFtp4 &= 1; RdFtp4 >>= b-8; RdFtp4 &= 1;
			if (LdFtp4 < RdFtp4) return true;
			if (LdFtp4 > RdFtp4) return false;
		}
		return false;
	}
};

struct diffpathlower_less
	: public std::binary_function<sha1differentialpath, sha1differentialpath, bool>
{
	bool operator()(const sha1differentialpath& _Left, const sha1differentialpath& _Right) const
	{
		if (_Left.tend() < _Right.tend() || _Left.tbegin() > _Right.tbegin())
			return true;
		if (_Left.tend() > _Right.tend() || _Left.tbegin() < _Right.tbegin())
			return false;
		if (_Left.path.size() < 4)
			throw std::runtime_error("Lower differential paths of insufficient size");
		unsigned t = _Left.tend() - 1;
		uint32 LdFt = 0 - _Left[t].getsdr().rotate_left(5).adddiff() - _Left[t-4].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp1 = 0 - _Left[t-3].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp2 = 0 - _Left[t-2].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp3 = 0 - _Left[t-1].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp4 = 0 - _Left[t].getsdr().rotate_left(30).adddiff();
		uint32 LdQtm1 = _Left[t-1].diff();
		uint32 LdQt = _Left[t].diff();
		uint32 RdFt = 0 - _Right[t].getsdr().rotate_left(5).adddiff() - _Right[t-4].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp1 = 0 - _Right[t-3].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp2 = 0 - _Right[t-2].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp3 = 0 - _Right[t-1].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp4 = 0 - _Right[t].getsdr().rotate_left(30).adddiff();
		uint32 RdQtm1 = _Right[t-1].diff();
		uint32 RdQt = _Right[t].diff();
		unsigned b = 0;
		while (b < 32+8) {
			if (b < 32) {
				if (_Left(t-3,(b+2)&31)!=_Right(t-3,(b+2)&31)) break;
				if (_Left(t-2,(b+2)&31)!=_Right(t-2,(b+2)&31)) break;
				if ((LdFt^RdFt)&(1<<b)) break;
				if ((LdQtm1^RdQtm1)&(1<<b)) break;
			}
			if (b >= 2 && b-2 < 32) {
				if ((LdFtp1^RdFtp1)&(1<<(b-2))) break;
				if ((LdQt^RdQt)&(1<<(b-2))) break;
			}
			if (b >= 4 && b-4 < 32 && ((LdFtp2^RdFtp2)&(1<<(b-4)))) break;
			if (b >= 6 && b-6 < 32 && ((LdFtp3^RdFtp3)&(1<<(b-6)))) break;
			if (b >= 8 && b-8 < 32 && ((LdFtp4^RdFtp4)&(1<<(b-8)))) break;
			++b;
		}
		if (b == 32+8)
			return false;
		if (b < 32) {
			LdFt >>= b; LdFt &= 1; RdFt >>= b; RdFt &= 1;
			if (LdFt < RdFt) return true;
			if (LdFt > RdFt) return false;
			LdQtm1 >>= b; LdQtm1 &= 1; RdQtm1 >>= b; RdQtm1 &= 1;
			if (LdQtm1 < RdQtm1) return true;
			if (LdQtm1 > RdQtm1) return false;
			if (_Left(t-3,(b+2)&31) < _Right(t-3,(b+2)&31)) return true;
			if (_Left(t-3,(b+2)&31) > _Right(t-3,(b+2)&31)) return false;
			if (_Left(t-2,(b+2)&31) < _Right(t-2,(b+2)&31)) return true;
			if (_Left(t-2,(b+2)&31) > _Right(t-2,(b+2)&31)) return false;
		}
		if (b >= 2 && b-2 < 32) {
			LdFtp1 >>= b-2; LdFtp1 &= 1; RdFtp1 >>= b-2; RdFtp1 &= 1;
			if (LdFtp1 < RdFtp1) return true;
			if (LdFtp1 > RdFtp1) return false;
			LdQt >>= b-2; LdQt &= 1; RdQt >>= b-2; RdQt &= 1;
			if (LdQt < RdQt) return true;
			if (LdQt > RdQt) return false;
		}
		if (b >= 4 && b-4 < 32) {
			LdFtp2 >>= b-4; LdFtp2 &= 1; RdFtp2 >>= b-4; RdFtp2 &= 1;
			if (LdFtp2 < RdFtp2) return true;
			if (LdFtp2 > RdFtp2) return false;
		}
		if (b >= 6 && b-6 < 32) {
			LdFtp3 >>= b-6; LdFtp3 &= 1; RdFtp3 >>= b-6; RdFtp3 &= 1;
			if (LdFtp3 < RdFtp3) return true;
			if (LdFtp3 > RdFtp3) return false;
		}
		if (b >= 8 && b-8 < 32) {
			LdFtp4 >>= b-8; LdFtp4 &= 1; RdFtp4 >>= b-8; RdFtp4 &= 1;
			if (LdFtp4 < RdFtp4) return true;
			if (LdFtp4 > RdFtp4) return false;
		}
		return false;
	}
};

struct diffpathlower_less_weight
	: public std::binary_function<sha1differentialpath, sha1differentialpath, bool>
{
	bool operator()(const sha1differentialpath& _Left, const sha1differentialpath& _Right) const
	{
		if (_Left.tend() < _Right.tend() || _Left.tbegin() > _Right.tbegin())
			return true;
		if (_Left.tend() > _Right.tend() || _Left.tbegin() < _Right.tbegin())
			return false;
		if (_Left.path.size() < 4)
			throw std::runtime_error("Lower differential paths of insufficient size");
		unsigned t = _Left.tend() - 1;
		unsigned lw = _Left[t+1].getsdr().hw() + _Left[t].getsdr().hw();
		unsigned rw = _Right[t+1].getsdr().hw() + _Right[t].getsdr().hw();
		if (lw != rw)
			return lw > rw;
		uint32 LdFt = 0 - _Left[t].getsdr().rotate_left(5).adddiff() - _Left[t-4].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp1 = 0 - _Left[t-3].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp2 = 0 - _Left[t-2].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp3 = 0 - _Left[t-1].getsdr().rotate_left(30).adddiff();
		uint32 LdFtp4 = 0 - _Left[t].getsdr().rotate_left(30).adddiff();
		uint32 LdQtm1 = _Left[t-1].diff();
		uint32 LdQt = _Left[t].diff();
		uint32 RdFt = 0 - _Right[t].getsdr().rotate_left(5).adddiff() - _Right[t-4].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp1 = 0 - _Right[t-3].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp2 = 0 - _Right[t-2].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp3 = 0 - _Right[t-1].getsdr().rotate_left(30).adddiff();
		uint32 RdFtp4 = 0 - _Right[t].getsdr().rotate_left(30).adddiff();
		uint32 RdQtm1 = _Right[t-1].diff();
		uint32 RdQt = _Right[t].diff();
		unsigned b = 0;
		while (b < 32+8) {
			if (b < 32) {
				if (_Left(t-3,(b+2)&31)!=_Right(t-3,(b+2)&31)) break;
				if (_Left(t-2,(b+2)&31)!=_Right(t-2,(b+2)&31)) break;
				if ((LdFt^RdFt)&(1<<b)) break;
				if ((LdQtm1^RdQtm1)&(1<<b)) break;
			}
			if (b >= 2 && b-2 < 32) {
				if ((LdFtp1^RdFtp1)&(1<<(b-2))) break;
				if ((LdQt^RdQt)&(1<<(b-2))) break;
			}
			if (b >= 4 && b-4 < 32 && ((LdFtp2^RdFtp2)&(1<<(b-4)))) break;
			if (b >= 6 && b-6 < 32 && ((LdFtp3^RdFtp3)&(1<<(b-6)))) break;
			if (b >= 8 && b-8 < 32 && ((LdFtp4^RdFtp4)&(1<<(b-8)))) break;
			++b;
		}
		if (b == 32+8)
			return false;
		if (b < 32) {
			LdFt >>= b; LdFt &= 1; RdFt >>= b; RdFt &= 1;
			if (LdFt < RdFt) return true;
			if (LdFt > RdFt) return false;
			LdQtm1 >>= b; LdQtm1 &= 1; RdQtm1 >>= b; RdQtm1 &= 1;
			if (LdQtm1 < RdQtm1) return true;
			if (LdQtm1 > RdQtm1) return false;
			if (_Left(t-3,(b+2)&31) < _Right(t-3,(b+2)&31)) return true;
			if (_Left(t-3,(b+2)&31) > _Right(t-3,(b+2)&31)) return false;
			if (_Left(t-2,(b+2)&31) < _Right(t-2,(b+2)&31)) return true;
			if (_Left(t-2,(b+2)&31) > _Right(t-2,(b+2)&31)) return false;
		}
		if (b >= 2 && b-2 < 32) {
			LdFtp1 >>= b-2; LdFtp1 &= 1; RdFtp1 >>= b-2; RdFtp1 &= 1;
			if (LdFtp1 < RdFtp1) return true;
			if (LdFtp1 > RdFtp1) return false;
			LdQt >>= b-2; LdQt &= 1; RdQt >>= b-2; RdQt &= 1;
			if (LdQt < RdQt) return true;
			if (LdQt > RdQt) return false;
		}
		if (b >= 4 && b-4 < 32) {
			LdFtp2 >>= b-4; LdFtp2 &= 1; RdFtp2 >>= b-4; RdFtp2 &= 1;
			if (LdFtp2 < RdFtp2) return true;
			if (LdFtp2 > RdFtp2) return false;
		}
		if (b >= 6 && b-6 < 32) {
			LdFtp3 >>= b-6; LdFtp3 &= 1; RdFtp3 >>= b-6; RdFtp3 &= 1;
			if (LdFtp3 < RdFtp3) return true;
			if (LdFtp3 > RdFtp3) return false;
		}
		if (b >= 8 && b-8 < 32) {
			LdFtp4 >>= b-8; LdFtp4 &= 1; RdFtp4 >>= b-8; RdFtp4 &= 1;
			if (LdFtp4 < RdFtp4) return true;
			if (LdFtp4 > RdFtp4) return false;
		}
		return false;
	}
};


#endif // MAIN_HPP
