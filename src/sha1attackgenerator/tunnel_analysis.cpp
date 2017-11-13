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

#include "main.hpp"
#include <hashclash/rng.hpp>
#include <hashclash/progress_display.hpp>
#include <hashclash/timer.hpp>
#include <boost/lexical_cast.hpp>

uint32 m[80];
uint32 m2[80];
uint32 Q[85];
uint32 Q2[85];
uint32 Qr30[85];
uint32 Q2r30[85];
uint32 dQ[85];
uint32 dT[85];
uint32 dR[85];
uint32 dF[80];

uint32 Qvaluemask[85];
uint32 Qvalue[85];
uint32 Qprev[85];
uint32 Qprev2[85];

void fill_tables(const sha1differentialpath& diffpath)
{
	for (int t = -4; t <= 80; ++t)
	{
		dF[t] = dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
		Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
		Qprev[offset+t] = Qprev2[offset+t] = 0;
	}

	// build tables
	for (int t = diffpath.tbegin(); t < diffpath.tend(); ++t)
	{
		dQ[offset+t] = diffpath[t].diff();
		Qprev[offset+t] = diffpath[t].prev() | diffpath[t].prevn();
		Qprev2[offset+t] = diffpath[t].prev2() | diffpath[t].prev2n();
		Qvaluemask[offset+t] = (~diffpath[t].set0()) | diffpath[t].set1()
							| diffpath[t].prev() | diffpath[t].prevn()
							| diffpath[t].prev2() | diffpath[t].prev2n();
		Qvalue[offset+t] = diffpath[t].set1() | diffpath[t].prevn() | diffpath[t].prev2n();
	}
	booleanfunction* F[4] = { &SHA1_F1_data, &SHA1_F2_data, &SHA1_F3_data, &SHA1_F4_data };
	for (int t = diffpath.tbegin()+4; t+1 < diffpath.tend(); ++t) {
		booleanfunction& FF = *F[t/20];
		uint32 dFt = 0;
		for (unsigned b = 0; b < 32; ++b)
			dFt += FF.outcome(diffpath[t-1][b],diffpath[t-2][(b+2)&31],diffpath[t-3][(b+2)&31])(0,b);
		dF[t] = dFt;
	}
	
	Q[0] = Qvalue[0]; Q2[0] = Q[0] + dQ[0];
	Q[1] = Qvalue[1]; Q2[1] = Q[1] + dQ[1];
	Q[2] = Qvalue[2]; Q2[2] = Q[2] + dQ[2];
	Q[3] = Qvalue[3]; Q2[3] = Q[3] + dQ[3];
	Q[4] = Qvalue[4]; Q2[4] = Q[4] + dQ[4];
}


unsigned analyze_tunnel(const sha1differentialpath& tunnel, bool verbose, bool deepanalyze)
{
#define AT_TEST_COUNT (1<<18) // 16
#define AT_MIN_PROB 0
#define AT_MIN_COUNT (AT_TEST_COUNT>>AT_MIN_PROB)

	vector< vector<unsigned> > qtbcnts(80, vector<unsigned>(32,0));
	static uint32 tcnts[80], qmask[80];
	static uint32 mmask[80], mset1[80];
	for (unsigned t = 0; t < 40; ++t) {
		tcnts[t] = 0;
		qmask[t] = maindiffpath[t].mask();
	}
	sha1differentialpath diffpath = tunnel;
	diffpath.get(-4);
	diffpath.get(80);
	fill_tables(diffpath);
	for (int t = 16; t < 80; ++t) {
		if (t-4 >= tunnel.tbegin() && t+1 < tunnel.tend()) {
			mmask[t] = ~tunnel.getme(t).mask;
			mset1[t] = tunnel.getme(t).set1conditions();
		} else {
			mmask[t] = ~uint32(0);
			mset1[t] = 0;
		}
	}
	

	for (unsigned k = 0; k < AT_TEST_COUNT; ++k) {
//		if (k == (AT_TEST_COUNT>>6) && double(tcnts[22])/double(k) >= 0.8) return 17;
		int t = -4;
		while (t < 17) {
			bool ok = false;
			for (unsigned j = 0; j < 3; ++j) {
				Q[offset+t] = (xrng128()&~Qvaluemask[offset+t]) ^ Qvalue[offset+t];
				if (t > -4)
					Q[offset+t] ^= Qprev[offset+t]&Q[offset+t-1];
				if (t > -3)
					Q[offset+t] ^= Qprev2[offset+t]&Q[offset+t-2];
				Q2[offset+t] = Q[offset+t] + dQ[offset+t];
				if (t >= 1) {
					m[t-1] = Q[offset+t] - (
						sha1_f1(Q[offset+t-2], rotate_left(Q[offset+t-3],30), rotate_left(Q[offset+t-4],30)) + sha1_ac[0] + rotate_left(Q[offset+t-1],5) + rotate_left(Q[offset+t-5],30)
						);
					m2[t-1] = Q2[offset+t] - (
						sha1_f1(Q2[offset+t-2], rotate_left(Q2[offset+t-3],30), rotate_left(Q2[offset+t-4],30)) + sha1_ac[0] + rotate_left(Q2[offset+t-1],5) + rotate_left(Q2[offset+t-5],30)
						);
					if ((m2[t-1]^m[t-1]) == diffpath.getme(t-1).mask) {
						ok = true;
						break;
					} else {
						if (verbose)
							cout << "(" << t << ":" << naf(m2[t-1]^m[t-1]) << "!=" << diffpath.getme(t-1) << ") " << flush;
						continue;
					}
				} else
					ok = true;
				break;
			}
			if (ok)
				++t;
			else
				t -= 4;
			if (t < -4) throw std::runtime_error("t < -4");
		}
		for (unsigned i = 16; i < 26; ++i) {
			m[i]=rotate_left(m[i-3] ^ m[i-8] ^ m[i-14] ^ m[i-16], 1);
			m2[i]=rotate_left(m2[i-3] ^ m2[i-8] ^ m2[i-14] ^ m2[i-16], 1);
			uint32 md = m[i]^m2[i];
			m[i] = (m[i]&mmask[i])^mset1[i];
			m2[i] = m[i] ^ md;
		}
		for (int t = 16; t < 20; ++t) {
			sha1_step_round1(t, Q, m);
			sha1_step_round1(t, Q2, m2);
		}
		for (int t = 20; t < 26; ++t) {
			sha1_step_round2(t, Q, m);
			sha1_step_round2(t, Q2, m2);
		}
		for (int t = 17; t < 27; ++t) {
			for (unsigned b = 0; b < 32; ++b)
				if ((Q2[offset+t]^Q[offset+t])&(1<<b))
					++qtbcnts[t][b];
		}
		int breakt = 80;
		for (int t = 17; t < 27; ++t) {
			if ((Q2[offset+t]^Q[offset+t])&(qmask[t]|0x00000000)) {
				breakt = t;
				for (int t1 = t; t1 < 27; ++t1)
					++tcnts[t1];
				break;
			}
		}
	}
	int tret = 80;
	for (int t = 17; t < 27; ++t)
	{
		if (tcnts[t] >= AT_MIN_COUNT) { tret = t; break; }
		if (double(tcnts[t]) >= double(AT_TEST_COUNT)*0.1) { tret = t; break; }
	}
	if (!deepanalyze || tret < mainparameters.mintunnelpov || tret == 80)
		return tret;
	cout << "Tunnel: " << endl;
	show_path(tunnel);
	cout << "  << ";
	for (int t1 = 17; t1 <= 25; ++t1)
		cout << t1 << ":" << ((double(100)*double(tcnts[t1]))/double(AT_TEST_COUNT)) << "% ";
	cout << ">>  " << flush;
	cout << tret << "b[";
	for (unsigned b = 0; b < 32; ++b)
		cout << b << ":" << unsigned((double(100)*double(qtbcnts[tret][b]))/double(AT_TEST_COUNT)) << "%" << (b==31?"]":",");
	cout << endl;
	return tret;
}










#define GOODBIT(i,t,b) ((goodqtb[i][t]>>b)&1)
#define BADBIT(i,t,b) ((badqtb[i][t]>>b)&1)
#define STAT_TEST(ba,go) ((double(ba)/double(ba+go))<0.6*(double(badqtb.size())/double(badqtb.size()+goodqtb.size())))
#define STAT_ANA(ba,go,s) if (STAT_TEST(ba,go)) { cout << "poss. bc: q_" << t << "[" << b << "]=" << s << ": \t" << (double(100*ba)/double(ba+go)) << "%" << endl; }
void stat_ana(vector< vector<uint32> >& goodqtb, vector< vector<uint32> >& badqtb
	, vector< vector<uint32> >& goodq2tb, vector< vector<uint32> >& badq2tb
	, int t, int b)
{
	unsigned g_zero = 0, g_one = 0, g_prev = 0, g_prevn = 0, g_prevr = 0, g_prevrn = 0, g_prev2 = 0, g_prev2n = 0;
	unsigned b_zero = 0, b_one = 0, b_prev = 0, b_prevn = 0, b_prevr = 0, b_prevrn = 0, b_prev2 = 0, b_prev2n = 0;
	for (unsigned i = 0; i < goodqtb.size(); ++i) {
		if (GOODBIT(i,t,b)==0) ++g_zero; else ++g_one;
		if (GOODBIT(i,t,b)==GOODBIT(i,t-1,b)) ++g_prev; else ++g_prevn;
		if (GOODBIT(i,t,b)==GOODBIT(i,t-1,((b+2)&31))) ++g_prevr; else ++g_prevrn;
		if (GOODBIT(i,t,b)==GOODBIT(i,t-2,((b+2)&31))) ++g_prev2; else ++g_prev2n;
	}
	for (unsigned i = 0; i < badqtb.size(); ++i) {
		if (BADBIT(i,t,b)==0) ++b_zero; else ++b_one;
		if (BADBIT(i,t,b)==BADBIT(i,t-1,b)) ++b_prev; else ++b_prevn;
		if (BADBIT(i,t,b)==BADBIT(i,t-1,((b+2)&31))) ++b_prevr; else ++b_prevrn;
		if (BADBIT(i,t,b)==BADBIT(i,t-2,((b+2)&31))) ++b_prev2; else ++b_prev2n;
	}
	double baseprob = double(badqtb.size())/double(badqtb.size()+goodqtb.size());
	STAT_ANA(b_zero,g_zero,"0");
	STAT_ANA(b_one,g_one,"1");
	STAT_ANA(b_prev,g_prev,"^");
	STAT_ANA(b_prevn,g_prevn,"!");
	STAT_ANA(b_prevr,g_prevr,"^r");
	STAT_ANA(b_prevrn,g_prevrn,"!r");
	STAT_ANA(b_prev2,g_prev2,"M");
	STAT_ANA(b_prev2n,g_prev2n,"#");
	if (b != 0) return;
	if (1) {
		map<sdr, pair<unsigned,unsigned> > dqtstat;
		for (unsigned i = 0; i < goodqtb.size(); ++i)
			++dqtstat[sdr(goodqtb[i][t],goodq2tb[i][t])].first;
		for (unsigned i = 0; i < badqtb.size(); ++i)
			++dqtstat[sdr(badqtb[i][t],badq2tb[i][t])].second;
		typedef map<sdr, pair<unsigned,unsigned> >::const_iterator dqtstatcit;
		vector<dqtstatcit> goodcits;
		unsigned goodcnt = 0, badcnt = 0;	
		for (dqtstatcit cit = dqtstat.begin(); cit != dqtstat.end(); ++cit)
			if ((double(cit->second.second)/double(cit->second.first+cit->second.second))<0.6*baseprob) {
				goodcits.push_back(cit);
				goodcnt += cit->second.first;
				badcnt += cit->second.second;
			}
		if (goodcits.size() && double(goodcnt) > 0.25*double(goodqtb.size()) && goodcits.size()*2 < goodqtb.size()) {
			cout << "good dQ_" << t << " " << double(badcnt*100)/double(badcnt+goodcnt) << "% bad #=" << goodcits.size() << ":\t";
			for (unsigned i = 0; i < goodcits.size(); ++i)
				if (double(goodcits[i]->second.first) > 0.01*double(goodqtb.size()))
				cout << (goodcits[i]->first) << "(" << 100.0*double(goodcits[i]->second.first)/double(goodcnt) << "%) ";
			cout << endl;
		}
	}
	if (1) {
		map<sdr, pair<unsigned,unsigned> > dftstat;
		if (t < 20) {
			for (unsigned i = 0; i < goodqtb.size(); ++i) {
				uint32 F1 = sha1_f1(goodqtb[i][t-1], rotate_left(goodqtb[i][t-2],30), rotate_left(goodqtb[i][t-3],30));
				uint32 F2 = sha1_f1(goodq2tb[i][t-1], rotate_left(goodq2tb[i][t-2],30), rotate_left(goodq2tb[i][t-3],30));
				++dftstat[sdr(F1,F2)].first;
			}
			for (unsigned i = 0; i < badqtb.size(); ++i) {
				uint32 F1 = sha1_f1(badqtb[i][t-1], rotate_left(badqtb[i][t-2],30), rotate_left(badqtb[i][t-3],30));
				uint32 F2 = sha1_f1(badq2tb[i][t-1], rotate_left(badq2tb[i][t-2],30), rotate_left(badq2tb[i][t-3],30));
				++dftstat[sdr(F1,F2)].second;
			}
		} else {
			for (unsigned i = 0; i < goodqtb.size(); ++i) {
				uint32 F1 = sha1_f2(goodqtb[i][t-1], rotate_left(goodqtb[i][t-2],30), rotate_left(goodqtb[i][t-3],30));
				uint32 F2 = sha1_f2(goodq2tb[i][t-1], rotate_left(goodq2tb[i][t-2],30), rotate_left(goodq2tb[i][t-3],30));
				++dftstat[sdr(F1,F2)].first;
			}
			for (unsigned i = 0; i < badqtb.size(); ++i) {
				uint32 F1 = sha1_f2(badqtb[i][t-1], rotate_left(badqtb[i][t-2],30), rotate_left(badqtb[i][t-3],30));
				uint32 F2 = sha1_f2(badq2tb[i][t-1], rotate_left(badq2tb[i][t-2],30), rotate_left(badq2tb[i][t-3],30));
				++dftstat[sdr(F1,F2)].second;
			}
		}
		typedef map<sdr, pair<unsigned,unsigned> >::const_iterator dftstatcit;
		vector<dftstatcit> goodcits;
		unsigned goodcnt = 0, badcnt = 0;	
		for (dftstatcit cit = dftstat.begin(); cit != dftstat.end(); ++cit)
			if ((double(cit->second.second)/double(cit->second.first+cit->second.second))<0.6*baseprob) {
				goodcits.push_back(cit);
				goodcnt += cit->second.first;
				badcnt += cit->second.second;
			}
		if (goodcits.size() && double(goodcnt) > 0.25*double(goodqtb.size()) && goodcits.size()*2 < goodqtb.size()) {
			cout << "good dF_" << t << " " << double(badcnt*100)/double(badcnt+goodcnt) << "% bad #=" << goodcits.size() << ":\t";
			for (unsigned i = 0; i < goodcits.size(); ++i)
				if (double(goodcits[i]->second.first) > 0.01*double(goodqtb.size()))
				cout << (goodcits[i]->first) << "(" << 100.0*double(goodcits[i]->second.first)/double(goodcnt) << "%) ";
			cout << endl;
		}
	}
}

unsigned analyze_bc_tunnel(const sha1differentialpath& tunnel)
{
#define ATLC_TEST_COUNT (1<<18) // 16
#define ATLC_MIN_PROB 0
#define ATLC_MIN_COUNT (ATLC_TEST_COUNT>>ATLC_MIN_PROB)
	unsigned tend = analyze_tunnel(tunnel, false, mainparameters.changebitstat);
	memset(m, 0, sizeof(uint32)*16);
	for (int i = 0; i < 16; ++i)
		if (i >= tunnel.tbegin() && i < tunnel.tend())
			m[i] = tunnel.getme(i).mask;
	for (unsigned i = 16; i < 32; ++i)
		m[i]=rotate_left(m[i-3] ^ m[i-8] ^ m[i-14] ^ m[i-16], 1);
//	if (tend < 21)
//		return tend;
	show_path(tunnel);
	cout << tend << endl;
	vector< vector<uint32> > dmes(80);
	unsigned tend2 = tend;
	uint64 mcnt = 1;
	for (unsigned i = 16; i < tend; ++i) {
		sdr tmp;
		tmp.mask = m[i];
		tmp.sign = 0;
		do {
			dmes[i].push_back(tmp.adddiff());
			tmp.sign += ~tmp.mask + 1;
			tmp.sign &= tmp.mask;
		} while (tmp.sign != 0);
		sort(dmes[i].begin(), dmes[i].end());
		dmes[i].erase( unique(dmes[i].begin(),dmes[i].end()), dmes[i].end());
		if ((mcnt * dmes[i].size()) > (1<<6)) {
			// don't let mcnt grow too big
			tend2 = i;
			break;
		}
		mcnt *= uint64(dmes[i].size());
	}	
	cout << mcnt << endl;
	uint32 qmask[80];
	for (unsigned t = 0; t < 40; ++t)
		qmask[t] = maindiffpath[t].mask();

	for (unsigned mi = 0; mi < mcnt; ++mi) {
		vector< vector<uint32> > goodqtb, badqtb, goodq2tb, badq2tb;
		vector<unsigned> dmesindex(80,0);
		unsigned mi2 = mi;
		cout << mi << ": ";
		for (unsigned i = 16; i < tend2; ++i)		
			if (dmes[i].size() > 1) {
				dmesindex[i] = mi2 % dmes[i].size();
				mi2 /= dmes[i].size();
				cout << i << naf(dmes[i][dmesindex[i]]);
			}
		cout << ": ";
		sha1differentialpath diffpath = tunnel;
		diffpath.get(-4);
		diffpath.get(80);
		fill_tables(diffpath);
		for (unsigned k = 0; k < ATLC_TEST_COUNT; ++k) {
			int t = -4;
			while (t < 17) {
				bool ok = false;
				for (unsigned j = 0; j < 3; ++j) {
					Q[offset+t] = (xrng128()&~Qvaluemask[offset+t]) ^ Qvalue[offset+t];
					if (t > -4)
						Q[offset+t] ^= Qprev[offset+t]&Q[offset+t-1];
					if (t > -3)
						Q[offset+t] ^= Qprev2[offset+t]&Q[offset+t-2];
					Q2[offset+t] = Q[offset+t] + dQ[offset+t];
					if (t >= 1) {
						m[t-1] = Q[offset+t] - (
							sha1_f1(Q[offset+t-2], rotate_left(Q[offset+t-3],30), rotate_left(Q[offset+t-4],30)) + sha1_ac[0] + rotate_left(Q[offset+t-1],5) + rotate_left(Q[offset+t-5],30)
							);
						m2[t-1] = Q2[offset+t] - (
							sha1_f1(Q2[offset+t-2], rotate_left(Q2[offset+t-3],30), rotate_left(Q2[offset+t-4],30)) + sha1_ac[0] + rotate_left(Q2[offset+t-1],5) + rotate_left(Q2[offset+t-5],30)
							);
						if ((m2[t-1]^m[t-1]) == diffpath.getme(t-1).mask) {
							ok = true;
							break;
						} else {
							continue;
						}
					} else
						ok = true;
					break;
				}
				if (ok)
					++t;
				else
					t -= 4;
				if (t < -4) throw std::runtime_error("t < -4");
			}
			for (unsigned i = 16; i < 26; ++i) {
				m[i]=rotate_left(m[i-3] ^ m[i-8] ^ m[i-14] ^ m[i-16], 1);
				if (i < tend2)
					m2[i] = m[i] + dmes[i][dmesindex[i]];
				else
					m2[i]=rotate_left(m2[i-3] ^ m2[i-8] ^ m2[i-14] ^ m2[i-16], 1);
			}
			for (int t = 16; t < 20; ++t) {
				sha1_step_round1(t, Q, m);
				sha1_step_round1(t, Q2, m2);
			}
			for (int t = 20; t < 26; ++t) {
				sha1_step_round2(t, Q, m);
				sha1_step_round2(t, Q2, m2);
			}
			for (int t = 17; t < 27; ++t) {
				if ((Q2[offset+t]^Q[offset+t])&(qmask[t]|0x00000000)) {
					static vector<uint32> tmpqtb(27);
					for (int t1 = 0; t1 < tmpqtb.size(); ++t1)
						tmpqtb[t1] = Q[4+t1];
					if (t <= tend)
						badqtb.push_back(tmpqtb);
					else
						goodqtb.push_back(tmpqtb);
					for (int t1 = 0; t1 < tmpqtb.size(); ++t1)
						tmpqtb[t1] = Q2[4+t1];
					if (t <= tend)
						badq2tb.push_back(tmpqtb);
					else
						goodq2tb.push_back(tmpqtb);
					break;
				}
			}
		}
		cout << goodqtb.size() << " vs " << badqtb.size() << ": " << (double(100)*double(badqtb.size()))/double(goodqtb.size()+badqtb.size()) << "% bad" << endl;
		for (unsigned t = 16-4; t < tend; ++t)
			for (unsigned b = 0; b < 32; ++b) {
				stat_ana(goodqtb, badqtb, goodq2tb, badq2tb, t, b);
			}
	}
	return tend;
}





void find_best_tunnel(vector< sha1differentialpath >& tunnels, const sha1differentialpath& mypath, const vector< vector<uint32> >& bitrels, bool force_lc = false) 
{
//#define ADD_ME_CARRIES 0
	unsigned ADD_ME_CARRIES = mainparameters.mecarry;
	if (mainparameters.tunnelfilterbitcondition)
		filter_tunnels_bitconditions(tunnels, mypath);
	static vector<sdr> sdrsvec[16];
	sdr metbu;
	static vector<unsigned> iQvec, ncvec;
	iQvec.assign(tunnels.size(),0);
	ncvec.assign(tunnels.size(),1<<20);
//	cout << tunnels.size() << endl;
	for (unsigned i = 0; i < tunnels.size(); ++i) {
		uint64 cnt = 1;
		for (int t = 0; t < 16; ++t) {
			sdrsvec[t].clear();
			if (t < tunnels[i].tbegin() || t >= tunnels[i].tend()) continue;
			metbu = tunnels[i].getme(t);
			table_sdrs(sdrsvec[t], metbu.adddiff(), hw(metbu.mask)+ADD_ME_CARRIES);
			unsigned j = 0;
			while (j < sdrsvec[t].size()) {
				if (sdrsvec[t][j].get(31) == -1) {
					swap(sdrsvec[t][j], sdrsvec[t][sdrsvec[t].size()-1]);
					sdrsvec[t].pop_back();
				} else
					++j;
			}
			cnt *= sdrsvec[t].size();
		}
		if (cnt == 0)
			throw std::runtime_error("cnt == 0 ??");
		unsigned bestiQ = 0, bestnc = 1<<20;
		uint64 bestc = cnt;
		unsigned basenc = tunnels[i].nrcond();
		for (uint64 c = 0; c < cnt; ++c) {
			uint64 cc = c;
			for (int t = 0; t < 16; ++t) {
				if (sdrsvec[t].size()>1) {
					uint64 k = cc % sdrsvec[t].size();
					cc /= sdrsvec[t].size();
					tunnels[i].getme(t) = sdrsvec[t][k];
				}
			}
			if (!filter_tunnel_bitrelations(tunnels[i],bitrels)) continue;
			unsigned iQ;
			if (force_lc)
				iQ = analyze_bc_tunnel(tunnels[i]);
			else
				iQ = analyze_tunnel(tunnels[i], false, mainparameters.changebitstat);
			if (iQ < bestiQ) continue;
			if (iQ > bestiQ) {
				bestiQ = iQ;
				bestnc = 1<<20;
				bestc = c;
			}
			unsigned nc = basenc;
			for (int t = tunnels[i].tbegin(); t < tunnels[i].tend(); ++t)
				nc += hw(tunnels[i].getme(t).mask&0x7FFFFFFF);
//			if (nc >= 15) continue;
			if (nc < bestnc) {
				bestnc = nc;
				bestc = c;
			}
		}
		if (bestc == cnt) {
			tunnels[i].clear();
			continue;
		}
		uint64 cc = bestc;
		for (int t = 0; t < 16; ++t) {
			if (sdrsvec[t].size()) {
				uint64 k = cc % sdrsvec[t].size();
				cc /= sdrsvec[t].size();
				tunnels[i].getme(t) = sdrsvec[t][k];
			}
		}
		iQvec[i] = bestiQ;
		ncvec[i] = bestnc;
//		if (bestiQ >= 23) break;
	}
	unsigned bestiQ = 0, bestnc = 1<<20, besti = tunnels.size(), minnc = 1<<20;
	for (unsigned i = 0; i < tunnels.size(); ++i) {
		if (tunnels[i].path.size()) {
			if (iQvec[i] > bestiQ) {
				bestiQ = iQvec[i];
				bestnc = ncvec[i];
				besti = i;
			} else if (iQvec[i] == bestiQ && ncvec[i] < bestnc) {
				bestnc = ncvec[i];
				besti = i;
			}
			if (ncvec[i] < minnc)
				minnc = ncvec[i];
		}
	}
	if (besti < tunnels.size()) {
//		if (bestiQ >= 21)
			cout << endl << bestiQ << "(" << bestnc << "," << minnc << ") " << flush;
		if (bestiQ >= 16) { //21
			//show_path(tunnels[besti]);
#if 0
			analyze_tunnel(tunnels[besti], false, true);
#else
			for (unsigned j = 0; j < tunnels.size(); ++j)
				if (iQvec[j] == iQvec[besti])
					analyze_tunnel(tunnels[j], false, true);
#endif
		} 
		sha1differentialpath bestpath;
		bestpath.swap(tunnels[besti]);
		tunnels.resize(1);
		bestpath.swap(tunnels[0]);
		
		if (bestiQ < 19) tunnels.clear();
		
	} else
		tunnels.clear();
}

void create_tunnel_stept(vector< sha1differentialpath >& tunnels, sha1differentialpath& tunnel, int t, bool force_lc = false)
{
	if (t >= 0 && t+1 >= tunnel.tend()) {
		create_tunnel_stept(tunnels, tunnel, t-1, force_lc);
		return;
	}

	if (t-4 < tunnel.tbegin() || t < 0) {
		tunnels.push_back(tunnel);
		return;
	}
	if (tunnel[t-1].diff() == 0 && tunnel[t-2].diff() == 0 && tunnel[t-3].diff() == 0) {
		// dFt = 0
		uint32 dme = tunnel[t+1].diff() - tunnel[t].getsdr().rotate_left(5).adddiff() - tunnel[t-4].getsdr().rotate_left(30).adddiff();
		tunnel.getme(t) = naf(dme);
		create_tunnel_stept(tunnels, tunnel, t-1, force_lc);
		return;
	}

	bf_conditions newconditions[32][2];
	bool bitactive[32];
	unsigned cnt = 1;
	for (unsigned b = 0; b < 32; ++b) {
		bitactive[b] = false;
		//sha1_f1(Q[offset+t-1], rotate_left(Q[offset+t-2],30), rotate_left(Q[offset+t-3],30))
		bf_outcome bo = SHA1_F1_data.outcome(tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31));
		if (bo.size() == 0) throw std::runtime_error("bo.size() == 0");
		if (bo.size() == 1) continue;
		if (bo.size() == 2) {
			if (/*b == 31 &&*/ bo[0] != bc_constant && bo[1] != bc_constant)
				continue;
			if (force_lc) {
				bf_conditions bf;
				if (bo[0] != bc_constant)
					bf = SHA1_F1_data.backwardconditions( tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31), bo[0] );
				else
					bf = SHA1_F1_data.backwardconditions( tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31), bo[1] );
				tunnel.setbitcondition(t-1,b,bf.first);
				tunnel.setbitcondition(t-2,(b+2)&31,bf.second);
				tunnel.setbitcondition(t-3,(b+2)&31,bf.third);
				continue;
			}
			newconditions[b][0] = SHA1_F1_data.backwardconditions( tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31), bo[0]);
			newconditions[b][1] = SHA1_F1_data.backwardconditions( tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31), bo[1]);
			cnt *= 2;
			bitactive[b] = true;
		} else { // bo.size() == 3
			bf_conditions bf = SHA1_F1_data.backwardconditions( tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31), bc_constant );
#if 1
			if (mainparameters.simpletunnels) {
				tunnel.setbitcondition(t-1,b,bf.first);
				tunnel.setbitcondition(t-2,(b+2)&31,bf.second);
				tunnel.setbitcondition(t-3,(b+2)&31,bf.third);
				continue;
			}
#endif
			if (force_lc) {
				tunnel.setbitcondition(t-1,b,bf.first);
				tunnel.setbitcondition(t-2,(b+2)&31,bc_prevn);
				tunnel.setbitcondition(t-3,(b+2)&31,bf.third);
				continue;
			}
			if (bf.second != bc_prev) throw std::runtime_error("bf.second != bc_prev");
			newconditions[b][0] = bf;
			bf.second = bc_prevn;
			newconditions[b][1] = bf;
			cnt *= 2;
			bitactive[b] = true;
		}
	}
	for (unsigned c = 0; c < cnt; ++c) {
		unsigned cc = c;
		uint32 dFt = 0;
		for (unsigned b = 0; b < 32; ++b) {
			if (bitactive[b]) {
				unsigned k = cc % 2;
				cc /= 2;
				tunnel.setbitcondition(t-1,b,newconditions[b][k].first);
				tunnel.setbitcondition(t-2,((b+2)&31),newconditions[b][k].second);
				tunnel.setbitcondition(t-3,((b+2)&31),newconditions[b][k].third);
			}
			dFt += SHA1_F1_data.outcome(tunnel(t-1,b), tunnel(t-2,(b+2)&31), tunnel(t-3,(b+2)&31))(0,b);
		}

		uint32 dme = tunnel[t+1].diff() - dFt - tunnel[t].getsdr().rotate_left(5).adddiff() - tunnel[t-4].getsdr().rotate_left(30).adddiff();
		tunnel.getme(t) = naf(dme);
		create_tunnel_stept(tunnels, tunnel, t-1, force_lc);
	}
}

void analyze_tunnels_diffpath(const sha1differentialpath& mypath, const vector< vector<uint32> >& bitrels, vector<sha1differentialpath> tunnelsin)
{
	cout << "Analyzing tunnels for differential path:" << endl;
	sha1differentialpath path = mypath;
	cleanup_path(path);
	show_path(path);
	maindiffpath = path;
	cout << "Message space determined by " << bitrels.size() << " bitrelations." << endl;

	bool deep_analysis = mainparameters.tunneldeepanalysis;
	
	if (tunnelsin.size()) {
		vector<sha1differentialpath> tunnelstmp;
		for (unsigned i = 0; i < tunnelsin.size(); ++i)
			create_tunnel_stept(tunnelstmp, tunnelsin[i], 15);
		find_best_tunnel(tunnelstmp, path, bitrels, deep_analysis);	
		return;
	}

	sha1differentialpath tunnelpath;
	vector< sha1differentialpath > tunnels, tunnelstmp, tunnelstmp2;

	cout << "Generating tunnels..." << endl;
	set<int> goodbits;
	set< pair<int,int> > good2bits;
#if 1
	for (int bit1 = 0; bit1 < 16*32; ++bit1) {
		int t1 = bit1 / 32;
		int b1 = bit1 % 32;
		tunnelpath.clear();
		tunnelpath.get(t1-4);
		if (t1 + 6 < 16)
			tunnelpath.get(t1+6);
		else
			tunnelpath.get(16);
		tunnelpath.setbitcondition(t1+1,b1,bc_plus);
		tunnelstmp.clear();
		create_tunnel_stept(tunnelstmp, tunnelpath, 15);
		find_best_tunnel(tunnelstmp, path, bitrels, deep_analysis);
		if (tunnelstmp.size()) {
			goodbits.insert(bit1);
			continue;
		}
	}
#if 1
	for (int bit1 = 0; bit1 < 16*32; ++bit1) {
		for (int bit2 = bit1+1; bit2 < 16*32; ++bit2) {
			if (goodbits.find(bit1) != goodbits.end()) continue;
			if (goodbits.find(bit2) != goodbits.end()) continue;
			int t1 = bit1 / 32;
			int b1 = bit1 % 32;
			int t2 = bit2/32;
			int b2 = bit2%32;
			tunnelpath.clear();
			tunnelpath.get(t1-4);
			if (t2 + 6 < 16)
				tunnelpath.get(t2+6);
			else
				tunnelpath.get(16);
			tunnelpath.setbitcondition(t1+1,b1,bc_plus);
			tunnelpath.setbitcondition(t2+1,b2,bc_plus);
			tunnelstmp.clear();
			create_tunnel_stept(tunnelstmp, tunnelpath, 15);
			find_best_tunnel(tunnelstmp, path, bitrels, deep_analysis);
			if (tunnelstmp.size())
				good2bits.insert(pair<int,int>(bit1,bit2));
		}
	}
#endif
//	exit(0);
#endif
	for (int bit1 = 16*32-1; bit1 >= 0; --bit1) {
//	for (int bit1 = 169; bit1 < 16*32; ++bit1) {
		cout << "[" << bit1 << "] " << flush;
		for (int bit2 = bit1; bit2 < 16*32; ++bit2) {
			for (int bit2sign = 0; bit2sign < 1; ++bit2sign)
			for (int bit3 = bit2; bit3 < 16*32; ++bit3)
			for (int bit3sign = 0; bit3sign < 1; ++bit3sign)
			{
				//if (bit2 != bit1 || bit3 != bit1) continue; // force single bit difference
				if (!mainparameters.threebittunnels && bit3 != bit2) continue; // force max 2 bit difference
				if (bit1 == bit2 && bit2sign == 1) continue;
				if (bit3 == bit2 && bit3sign != bit2sign) continue;
				if (bit1 == bit2 && bit3 > bit2) continue;
				if (goodbits.find(bit1) != goodbits.end()) continue;
				if (goodbits.find(bit2) != goodbits.end()) continue;
				if (goodbits.find(bit3) != goodbits.end()) continue;
				if (good2bits.find(pair<int,int>(bit1,bit2)) != good2bits.end()) continue;
				if (good2bits.find(pair<int,int>(bit1,bit3)) != good2bits.end()) continue;
				if (good2bits.find(pair<int,int>(bit2,bit3)) != good2bits.end()) continue;
				
				//cout << "[" << bit1 << "," << bit2 << "] " << flush;
				int t1 = bit1 / 32;
				int b1 = bit1 % 32;
				int t2 = bit2 / 32;
				int b2 = bit2 % 32;
				int t3 = bit3 / 32;
				int b3 = bit3 % 32;
				tunnelpath.clear();
				tunnelpath.get(t1-4);
				if (t3 + 6 < 16)
					tunnelpath.get(t3+6);
				else
					tunnelpath.get(16);
				tunnelpath.setbitcondition(t1+1,b1,bc_plus);
				if (bit2sign == 0)
					tunnelpath.setbitcondition(t2+1,b2,bc_plus);
				else
					tunnelpath.setbitcondition(t2+1,b2,bc_minus);
				if (bit3sign == 0)
					tunnelpath.setbitcondition(t3+1,b3,bc_plus);
				else
					tunnelpath.setbitcondition(t3+1,b3,bc_minus);
				tunnelstmp.clear();
				create_tunnel_stept(tunnelstmp, tunnelpath, 15);
				find_best_tunnel(tunnelstmp, path, bitrels, deep_analysis);
			}
		}
	}
	exit(0);
}














































struct ihv_type {
	uint32 ihv[5];
	bool operator<(const ihv_type& r) const {
		for (int i = 0; i < 4; ++i) {
			if (ihv[i] < r.ihv[i]) return true;
			if (ihv[i] > r.ihv[i]) return false;
		}
		if (ihv[4] < r.ihv[4]) return true;
		return false;
	}
	bool operator==(const ihv_type& r) const {
		for (int i = 0; i < 5; ++i)
			if (ihv[i] != r.ihv[i]) return false;
		return true;
	}
};

void analyze_indepsection_prob()
{
	vector< vector< vector<uint32> > > pathbitrelationsmatrix;
	mespace_to_pathbitrelationsmatrix(mainmespace, pathbitrelationsmatrix);
	sha1differentialpath path = maindiffpath;
	int tbegin = 16;
	for (; tbegin < 80; ++tbegin)
		if (0 == (path[tbegin].diff() | path[tbegin-1].diff() | path[tbegin-2].diff() | path[tbegin-3].diff() | path[tbegin-4].diff()))
			break;
	while (tbegin < 80) {
		int tend = tbegin+1;
		for (; tend < 80; ++tend)
			if (0 == (path[tend].diff() | path[tend-1].diff() | path[tend-2].diff() | path[tend-3].diff() | path[tend-4].diff()))
				break;
		if (0 != (path[tend].diff() | path[tend-1].diff() | path[tend-2].diff() | path[tend-3].diff() | path[tend-4].diff()))
			break;
		if (tend == tbegin+1) {
			tbegin = tend;
			continue;
		}
		cout << "Prob t=[" << tbegin << "-" << tend << "): \t" << flush;
#if 0
		tbegin = tend; cout << endl; continue;
#endif
		uint64 cnt = 0;
		uint32 okcnt = 0;
		while (okcnt < (1<<15)) {
			++cnt; if (hw(cnt)==1) cout << cnt << " " << flush;
			random_me(m, pathbitrelationsmatrix);
			for (int t = 0; t < 16; ++t)
				m2[t] = m[t] ^ path.getme(t).mask;
			for (int t = 16; t < 80; ++t)
				m2[t]=rotate_left(m2[t-3] ^ m2[t-8] ^ m2[t-14] ^ m2[t-16], 1);
			for (int t = tbegin-4; t <= tbegin; ++t)
				Q2[offset+t] = Q[offset+t] = xrng64();
			for (int t = tbegin; t < tend; ++t) {
				sha1_step(t, Q, m);
				sha1_step(t, Q2, m2);
			}
			bool ok = true;
			for (int t = tend-4; t <= tend; ++t)
				if (Q2[offset+t] != Q[offset+t]) {
					ok = false;
					break;
				}
			if (ok)
				if (hw(++okcnt)==1) cout << "[" << okcnt << "] " << flush;
		}
		cout << " prob=" << log(double(okcnt)/double(cnt))/log(2.0) << endl;
		
		tbegin = tend;
	}		
	cout << "Prob t=[" << tbegin << "-" << 80 << "):" << endl;
	vector<ihv_type> target_dihvs(20);
	vector<ihv_type> dihvs(1<<28);
	progress_display pd(dihvs.size());
	for (unsigned i = 0; i < dihvs.size(); ++i,++pd) {
		random_me(m, pathbitrelationsmatrix);
		for (int t = 0; t < 16; ++t)
			m2[t] = m[t] ^ path.getme(t).mask;
		for (int t = 16; t < 80; ++t)
			m2[t]=rotate_left(m2[t-3] ^ m2[t-8] ^ m2[t-14] ^ m2[t-16], 1);
		for (int t = tbegin-4; t <= tbegin; ++t)
			Q2[offset+t] = Q[offset+t] = xrng64();
		for (int t = -4; t <= 0; ++t)
			Q2[offset+t] = Q[offset+t] = xrng64();
		for (int t = tbegin; t < 80; ++t) {
			sha1_step(t, Q, m);
			sha1_step(t, Q2, m2);
		}
		dihvs[i].ihv[0] = rotate_left(Q2[offset+80-4],30)-rotate_left(Q[offset+80-4],30);
		dihvs[i].ihv[1] = rotate_left(Q2[offset+80-3],30)-rotate_left(Q[offset+80-3],30);
		dihvs[i].ihv[2] = rotate_left(Q2[offset+80-2],30)-rotate_left(Q[offset+80-2],30);
		dihvs[i].ihv[3] = Q2[offset+80-1]-Q[offset+80-1];
		dihvs[i].ihv[4] = Q2[offset+80]-Q[offset+80];
	}
	sort(dihvs.begin(), dihvs.end());
	unsigned i = 0;
	vector<unsigned> bestcnts;
	while (i < dihvs.size()) {
		unsigned j = i+1;
		while (j < dihvs.size() && dihvs[j] == dihvs[i]) ++j;
		bestcnts.push_back(j-i);
		i = j;
	}
	cout << "Prob t=[" << tbegin << "-" << 80 << "): " << flush;
	sort(bestcnts.begin(),bestcnts.end());
	unsigned bestcnt = bestcnts[bestcnts.size()-1];
#if 1
	target_dihvs.clear();
	i = 0;
	while (i < dihvs.size()) {
		unsigned j = i+1;
		while (j < dihvs.size() && dihvs[j] == dihvs[i]) ++j;
		if ((j-i)*2 >= bestcnt) {
			cout << "\tprob=" << log(double(j-i)/double(dihvs.size()))/log(2.0) << ":\t";
			for (unsigned k = 0; k < 5; ++k)
				cout << naf(dihvs[i].ihv[k]) << " ";
			cout << endl;
			target_dihvs.push_back(dihvs[i]);
		}
		i = j;
	}
	{ vector<ihv_type> tmptmp(1); dihvs.swap(tmptmp); } // free memory
	uint64 cnt = 0, okcnt = 0;
	ihv_type dihv;
	while (true) {
		random_me(m, pathbitrelationsmatrix);
		for (int t = 0; t < 16; ++t)
			m2[t] = m[t] ^ path.getme(t).mask;
		for (int t = 16; t < 80; ++t)
			m2[t]=rotate_left(m2[t-3] ^ m2[t-8] ^ m2[t-14] ^ m2[t-16], 1);
		for (int t = tbegin-4; t <= tbegin; ++t)
			Q2[offset+t] = Q[offset+t] = xrng64();
		for (int t = -4; t <= 0; ++t)
			Q2[offset+t] = Q[offset+t] = xrng64();
		for (int t = tbegin; t < 80; ++t) {
			sha1_step(t, Q, m);
			sha1_step(t, Q2, m2);
		}
		dihv.ihv[0] = rotate_left(Q2[offset+80-4],30)-rotate_left(Q[offset+80-4],30);
		dihv.ihv[1] = rotate_left(Q2[offset+80-3],30)-rotate_left(Q[offset+80-3],30);
		dihv.ihv[2] = rotate_left(Q2[offset+80-2],30)-rotate_left(Q[offset+80-2],30);
		dihv.ihv[3] = Q2[offset+80-1]-Q[offset+80-1];
		dihv.ihv[4] = Q2[offset+80]-Q[offset+80];
		++cnt;
		for (unsigned k = 0; k < target_dihvs.size(); ++k)
			if (dihv == target_dihvs[k]) {
				++okcnt;
				break;
			}
		if (hw(cnt)+hw(cnt>>32)==1 && okcnt > 1) {
			cout << "[" << cnt << ":p=" << log(double(okcnt)/double(cnt))/log(2.0) << "]" << flush;
		}
	}
#else	
	for (unsigned i = 0; i < bestcnts.size(); ++i)
		if (bestcnts[i]*2 >= bestcnt)
			cout << " prob=" << log(double(bestcnts[i])/double(dihvs.size()))/log(2.0) << flush;
#endif
	cout << endl;
	exit(0);
}
