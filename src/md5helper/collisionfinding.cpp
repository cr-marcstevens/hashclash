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

#include <cmath>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <unordered_map>

#define MD5DETAIL_INLINE_IMPL
#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/progress_display.hpp>

#include <boost/lexical_cast.hpp>

#include "main.hpp"

//#define CPUPERFORMANCE
#ifdef CPUPERFORMANCE
#include <hashclash/cpuperformance.hpp>
//uint64 cpu_step_t[64];
#define UPDATE(s) update_performance_counter __update(cpu_step_t[s]);
#ifdef __GNUC__
#include <sched.h>
#endif
#else
#define UPDATE(s)
#endif

uint64 tendcount = 0, t61count = 0;
//vector<uint64> testcounts(1<<20,0);
//#define DOTESTCOUNTS
#ifdef DOTESTCOUNTS
#define TESTCOUNT(s) ++testcounts[s];
#else
#define TESTCOUNT(s)
#endif

/*
differentialpath diffpath;
uint32 m_diff[16];

vector< triple<uint32,uint32,uint32> > Q1Q2m0withm1ok;
vector< triple<uint32,uint32,uint32> >::const_iterator Q1Q2m0it, Q1Q2m0itend;

const int offset = 3;
uint32 m[16];
uint32 m2[16];
uint32 Q[68];
uint32 Q2[68];

uint32 dQ[68];
uint32 dT[68];
uint32 dR[68];

uint32 Qvaluemask[68];
uint32 Qvalue[68];
uint32 Qprev[68];
uint32 Qprev2[68];

uint32 Q5m5tunnel=0, Q4m5tunnel=0, Q4m4tunnel=0;
uint32 Q10m10tunnel=0, Q9m10tunnel=0, Q9m9tunnel=0;
uint32 Q8Q12m15tunnel=0, Q14m6Q3m5tunnel = 0, Q14Q3m14tunnel = 0;
*/

inline bool checkrotation(uint32 Qt, uint32 Qtp1, uint32 dT, uint32 dR, unsigned rc)
{
	uint32 R1 = Qtp1-Qt;
	uint32 R2 = R1 + dR;
	uint32 T1 = rotate_right(R1, rc);
	uint32 T2 = rotate_right(R2, rc);
	return (T2-T1 == dT);
}

int collisionfinding_thread::verifyconds()
{
	int badt = 64;
	for (int t = -2; t < 64; ++t)
	{
		switch (t) {
			case 3:
			case 14:
				if ((Qvalue[offset+t] ^ (Q[offset+t]&Qvaluemask[offset+t]) ^ (Qprev[offset+t]&Q[offset+t-1])) & ~(Q14m6Q3m5tunnel|Q14Q3m14tunnel))
					if (t < badt) badt = t;
				break;
			case 8:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&~Q8Q12m15tunnel)
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;
			case 12:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&rotate_left(~Q8Q12m15tunnel,22))
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;
			case 4:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&~Q4m4tunnel&~Q4m5tunnel)
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;

			case 5:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&~Q5m5tunnel)
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;

			case 9:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&~Q9m9tunnel&~Q9m10tunnel)
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;

			case 10:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t]&~Q10m10tunnel)
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;

			default:
				if (Qvalue[offset+t] != 
					((Q[offset+t]&Qvaluemask[offset+t])
						^(Qprev[offset+t]&Q[offset+t-1])))
					if (t < badt) badt = t;
				break;
		}
	}
	uint32 block[16];
	uint32 block2[16];
	for (int k = 0; k < 16; ++k)
	{
		uint32 F = md5_ff(Q[offset+k], Q[offset+k-1], Q[offset+k-2]);
		uint32 T = F + Q[offset+k-3] + md5_ac[k];
		uint32 R = Q[offset+k+1] - Q[offset+k];
		block[k] = rotate_right(R, md5_rc[k]) - T;
		block2[k] = block[k] + m_diff[k];
	}
	for (int k = 16; k < 32; ++k)
	{
		uint32 F = md5_gg(Q[offset+k], Q[offset+k-1], Q[offset+k-2]);
		uint32 T = F + Q[offset+k-3] + md5_ac[k] + block[md5_wt[k]];
		uint32 R = rotate_left(T, md5_rc[k]);
		if (Q[offset+k+1] != Q[offset+k] + R)
			if (k+1 < badt) badt = (k+1);
	}
	for (int k = 0; k < 16; ++k)
	{
		uint32 F = md5_ff(Q2[offset+k], Q2[offset+k-1], Q2[offset+k-2]);
		uint32 T = F + Q2[offset+k-3] + md5_ac[k] + block2[md5_wt[k]];
		uint32 R = rotate_left(T, md5_rc[k]);
		Q2[offset+k+1] = Q2[offset+k] + R;
		if (Q2[offset+k+1]-Q[offset+k+1] != dQ[offset+k+1])
			if (k+1 < badt) badt = k+1;
	}
	for (int k = 16; k < 32; ++k)
	{
		uint32 F = md5_gg(Q2[offset+k], Q2[offset+k-1], Q2[offset+k-2]);
		uint32 T = F + Q2[offset+k-3] + md5_ac[k] + block2[md5_wt[k]];
		uint32 R = rotate_left(T, md5_rc[k]);
		Q2[offset+k+1] = Q2[offset+k] + R;
		if (Q2[offset+k+1]-Q[offset+k+1] != dQ[offset+k+1])
			if (k+1 < badt) badt = k+1;
	}
	
	return badt;
}

void collisionfinding_thread::do_step26()
{
	UPDATE(26);
	uint32 block[16], block2[16];
	for (int k = 0; k < 16; ++k)
	{
		uint32 F = md5_ff(Q[offset+k], Q[offset+k-1], Q[offset+k-2]);
		uint32 T = F + Q[offset+k-3] + md5_ac[k];
		uint32 R = Q[offset+k+1] - Q[offset+k];
		block[k] = rotate_right(R, md5_rc[k]) - T;
		block2[k] = block[k] + m_diff[k];
	}
	uint32 ihv[4], ihv2[4];
	ihv[0] = Q[0]; ihv[1] = Q[3]; ihv[2] = Q[2]; ihv[3] = Q[1];
	ihv2[0] = Q2[0]; ihv2[1] = Q2[3]; ihv2[2] = Q2[2]; ihv2[3] = Q2[1];
	md5compress(ihv, block);
	md5compress(ihv2, block2);
	{
	boost::lock_guard<boost::mutex> lock(mut);
	if (ihv2[0]-ihv[0] == dQ[0] + dQ[offset+61]) {
		++t61count;
		if (hw(t61count)==1)
			cout << tendcount << " " << t61count << endl;
	} else return;

	if (ihv2[1]-ihv[1] != dQ[3] + dQ[offset+64]) return;
	if (ihv2[2]-ihv[2] != dQ[2] + dQ[offset+63]) return;
	if (ihv2[3]-ihv[3] != dQ[1] + dQ[offset+62]) return;

	uint32 x = xrng128();
	string filename1 = workdir + "/coll1_" + boost::lexical_cast<string>(x);
	string filename2 = workdir + "/coll2_" + boost::lexical_cast<string>(x);
	ofstream of1(filename1.c_str(), ios::binary);
	ofstream of2(filename2.c_str(), ios::binary);
	if ((!of1) || (!of2)) {
		cerr << "Cannot open output files!" << endl;
	}
	save_block(of1, block);
	save_block(of2, block2);
	of1.close();
	of2.close();
	cout << "Block 1: " << filename1 << endl;
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout << hex;
			cout.width(2); cout.fill('0');
			cout << ((block[k]>>(c*8))&0xFF) << " ";
			if ((k&3)==3 && c == 3) cout << endl;
		}
	cout << dec << "Block 2: " << filename2 << endl;
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout << hex;
			cout.width(2); cout.fill('0');
			cout << ((block2[k]>>(c*8))&0xFF) << " ";
			if ((k&3)==3 && c == 3) cout << endl;
		}
	cout << dec;
	cout << "Found collision!" << endl;
	}
	exit(0);
	//throw std::runtime_error("Found collision!");
}

void collisionfinding_thread::do_step25()
{
	UPDATE(25);

	uint32 Q14tmask = Q14Q3m14tunnel;
//	uint32 Q14tadd = ~Q14tmask + 1;
	uint32 Q14value = Q[offset+14];
	uint32 Q3value = Q[offset+3];

	uint32 pT25 = md5_gg(Q[offset+25],Q[offset+24],Q[offset+23]) + Q[offset+22] + md5_ac[25];

	uint32 Q14tcur = 0;
	do {
		Q14tcur -= 1; Q14tcur &= Q14tmask;
		Q[offset+14] = Q14value + Q14tcur;
		Q[offset+3] = Q3value + Q14tcur;

		uint32 R14 = Q[offset+15]-Q[offset+14];
		m[14] = rotate_right(R14,md5_rc[14]) - md5_ff(Q[offset+14],Q[offset+13],Q[offset+12]) - Q[offset+11] - md5_ac[14];
		uint32 T25 = pT25 + m[14];
		uint32 R25 = rotate_left(T25, md5_rc[25]);

		TESTCOUNT(16);
		Q[offset+26] = Q[offset+25] + R25;
		if (Qvalue[offset+26] != 
			((Q[offset+26]&Qvaluemask[offset+26])
				^(Qprev[offset+26]&Q[offset+25])))
			continue;
		TESTCOUNT(17);
		uint32 T25b = T25 + dT[offset+25];
		uint32 R25b = rotate_left(T25b, md5_rc[25]);
		if (R25b-R25 != dR[offset+25]) 
			continue;
		TESTCOUNT(18);

		mut.lock();
		++tendcount;
		if (hw(uint32(tendcount))==1)
			cout << tendcount << " " << t61count << endl;
		mut.unlock();

		uint32 R3 = Q[offset+4]-Q[offset+3];
		m[3] = rotate_right(R3,md5_rc[3]) - md5_ff(Q[offset+3],Q[offset+2],Q[offset+1])
			- Q[offset+0] - md5_ac[3];
		uint32 T26 = md5_gg(Q[offset+26],Q[offset+25],Q[offset+24]) + Q[offset+23] + md5_ac[26] + m[3];
		uint32 R26 = rotate_left(T26, 14);
		Q[offset+27] = Q[offset+26] + R26;
		if (Qvalue[offset+27] != 
			((Q[offset+27]&Qvaluemask[offset+27])
				^(Qprev[offset+27]&Q[offset+26])))
			continue;

		uint32 R8 = Q[offset+9]-Q[offset+8];
		m[8] = rotate_right(R8,md5_rc[8]) - md5_ff(Q[offset+8],Q[offset+7],Q[offset+6])
			- Q[offset+5] - md5_ac[8];
		uint32 T27 = md5_gg(Q[offset+27],Q[offset+26],Q[offset+25]) + Q[offset+24] + md5_ac[27] + m[8];
		uint32 R27 = rotate_left(T27, 20);
		Q[offset+28] = Q[offset+27] + R27;
		if (Qvalue[offset+28] != 
			((Q[offset+28]&Qvaluemask[offset+28])
				^(Qprev[offset+28]&Q[offset+27])))
			continue;

#if 0
		if (verifyconds() <= 28)
			cerr << verifyconds() << endl;
#endif
		do_step26();

	} while (Q14tcur != 0);
}

void collisionfinding_thread::do_step24()
{
	UPDATE(24);

	uint32 Q9tmask = (Q9m10tunnel|Q9m9tunnel) & ~Q[offset+10];
//	uint32 Q9tadd = ~Q9tmask + 1;
	uint32 Q9value = Q[offset+9];

	uint32 pT24 = md5_gg(Q[offset+24],Q[offset+23],Q[offset+22]) + Q[offset+21] + md5_ac[24];

	uint32 Q9tcur = 0;
	do {
		Q9tcur -= 1; Q9tcur &= Q9tmask;
		Q[offset+9] = Q9tcur ^ Q9value;

		uint32 R9 = Q[offset+10]-Q[offset+9];
		m[9] = rotate_right(R9,md5_rc[9]) - md5_ff(Q[offset+9],Q[offset+8],Q[offset+7])
			- Q[offset+6] - md5_ac[9];
		
		uint32 T24 = pT24 + m[9];
		uint32 R24 = rotate_left(T24, md5_rc[24]);

		TESTCOUNT(13);
		Q[offset+25] = Q[offset+24] + R24;
		if (Qvalue[offset+25] != 
			((Q[offset+25]&Qvaluemask[offset+25])
				^(Qprev[offset+25]&Q[offset+24])))
			continue;
		TESTCOUNT(14);

		uint32 T24b = T24 + dT[offset+24];
		uint32 R24b = rotate_left(T24b, md5_rc[24]);
		if (R24b-R24 != dR[offset+24]) 
			continue;
		TESTCOUNT(15);

#if 0
		if (verifyconds() <= 25)
			cerr << verifyconds() << endl;
#else
		do_step25();
#endif
	} while (Q9tcur != 0);
}

void collisionfinding_thread::do_step23()
{
	UPDATE(23);

	//uint32 Q4tmask = (Q4m5tunnel|Q4m4tunnel) & ~Q[offset+5];
	uint32 Q4tmask = Q4m4tunnel;
//	uint32 Q4tadd = ~Q4tmask + 1;
	uint32 Q4value = Q[offset+4];

	uint32 pT23 = md5_gg(Q[offset+23],Q[offset+22],Q[offset+21]) + Q[offset+20] + md5_ac[23];

	uint32 Q4tcur = 0;
	do {
		Q4tcur -= 1; Q4tcur &= Q4tmask;
		Q[offset+4] = Q4tcur ^ Q4value;

		uint32 R4 = Q[offset+5]-Q[offset+4];
		m[4] = rotate_right(R4,md5_rc[4]) - md5_ff(Q[offset+4],Q[offset+3],Q[offset+2])
			- Q[offset+1] - md5_ac[4];
		
		uint32 T23 = pT23 + m[4];
		uint32 R23 = rotate_left(T23, md5_rc[23]);

		TESTCOUNT(10);
		Q[offset+24] = Q[offset+23] + R23;
		if (Qvalue[offset+24] != 
			((Q[offset+24]&Qvaluemask[offset+24])
				^(Qprev[offset+24]&Q[offset+23])))
			continue;
		TESTCOUNT(11);

		uint32 T23b = T23 + dT[offset+23];						
		uint32 R23b = rotate_left(T23b, md5_rc[23]);
		if (R23b-R23 != dR[offset+23]) 
			continue;
		TESTCOUNT(12);

#if 0
		if (verifyconds() <= 24)
			cerr << verifyconds() << endl;
#else
		do_step24();
#endif

	} while (Q4tcur != 0);
}


void collisionfinding_thread::do_step21()
{
	UPDATE(21);
	uint32 Q12val = Q[offset+12];

	uint32 Q10tmask = Q10m10tunnel;
//	uint32 Q10tadd = ~Q10tmask + 1;
	uint32 Q10value = Q[offset+10];
	uint32 Q9tmask = (Q9m10tunnel|Q9m9tunnel) & Q[offset+10];
//	uint32 Q9tadd = ~Q9tmask + 1;
	uint32 Q9value = Q[offset+9];

	uint32 Q8tmask = Q8Q12m15tunnel;
//	uint32 Q8tadd = ~Q8tmask + 1;
	uint32 Q8value = Q[offset+8];

	uint32 R15 = Q[offset+16]-Q[offset+15];
	uint32 T15 = rotate_right(R15,md5_rc[15]);
	m[15] = T15 - md5_ff(Q[offset+15],Q[offset+14],Q[offset+13]) - Q[offset+12] - md5_ac[15];

	uint32 pT21 = md5_gg(Q[offset+21],Q[offset+20],Q[offset+19]) + Q[offset+18] + md5_ac[21];
	uint32 pT22 = Q[offset+19] + md5_ac[22] + m[15];

	uint32 Q10tcur = 0;
	do {
		Q10tcur -= 1; Q10tcur &= Q10tmask;
		Q[offset+10] = Q10tcur ^ Q10value;

		uint32 R10 = Q[offset+11]-Q[offset+10];
		uint32 pT10 = rotate_right(R10,md5_rc[10]) - md5_ac[10] - Q[offset+7];

		uint32 Q9tcur = 0;
		do {
			Q9tcur -= 1; Q9tcur &= Q9tmask;
			Q[offset+9] = Q9tcur ^ Q9value;

			m[10] = pT10 - md5_ff(Q[offset+10],Q[offset+9],Q[offset+8]);
			uint32 T21 = pT21 + m[10];
			uint32 R21 = rotate_left(T21, md5_rc[21]);

			TESTCOUNT(3);
			Q[offset+22] = Q[offset+21] + R21;
			if (Qvalue[offset+22] != 
				((Q[offset+22]&Qvaluemask[offset+22])
					^(Qprev[offset+22]&Q[offset+21])))
				continue;

//			TESTCOUNT(4);
			uint32 T21b = T21 + dT[offset+21];						
			uint32 R21b = rotate_left(T21b, md5_rc[21]);
			if (R21b-R21 != dR[offset+21]) 
				continue;
//			TESTCOUNT(5);

			UPDATE(22);
			uint32 Q8tcur = 0;
			do {
				Q8tcur -= 1; Q8tcur &= Q8tmask;
				Q[offset+8] = Q8tcur ^ Q8value;

				uint32 T11 = md5_ff(Q[offset+11],Q[offset+10],Q[offset+9])
								+Q[offset+8] + m[11] + md5_ac[11];
				uint32 R11 = rotate_left(T11, md5_rc[11]);
				uint32 Q12 = Q[offset+11] + R11;
//				TESTCOUNT(6);
				if (Qvalue[offset+12] != 
					((Q12&Qvaluemask[offset+12]&~rotate_left(Q8tmask,22))
						^(Qprev[offset+12]&Q[offset+11])))
					continue;
				if (Qvalue[offset+13] != 
					((Q[offset+13]&Qvaluemask[offset+13])
						^(Qprev[offset+13]&Q12)))
					continue;
//				TESTCOUNT(7);
				Q[offset+12] = Q12;
				uint32 R15 = Q[offset+16]-Q[offset+15];
				m[15] = rotate_right(R15, md5_rc[15])
							- md5_ff(Q[offset+15], Q[offset+14], Q[offset+13])
							- Q[offset+12]
							- md5_ac[15];

				//uint32 T22 = md5_gg(Q[offset+22], Q[offset+21], Q[offset+20]) + pT22;
				uint32 T22 = md5_gg(Q[offset+22], Q[offset+21], Q[offset+20])
								+ Q[offset+19]
								+ m[15]
								+ md5_ac[22];
				uint32 R22 = rotate_left(T22, md5_rc[22]);
				Q[offset+23] = Q[offset+22] + R22;
				if (Qvalue[offset+23] != 
					((Q[offset+23]&Qvaluemask[offset+23])
						^(Qprev[offset+23]&Q[offset+22])))
					continue;
//				TESTCOUNT(8);
				uint32 T22b = T22 + dT[offset+22];						
				uint32 R22b = rotate_left(T22b, md5_rc[22]);
				if (R22b-R22 != dR[offset+22]) 
					continue;
//				TESTCOUNT(9);

#if 0
				if (verifyconds() <= 23)
					cerr << verifyconds() << endl;
#else
				do_step23();
#endif
			} while (Q8tcur != 0);
		} while (Q9tcur != 0);
	} while (Q10tcur != 0);
	Q[offset+12] = Q12val;
}

void collisionfinding_thread::do_step20()
{
	UPDATE(20);

	uint32 Q5tmask = Q5m5tunnel;
//	uint32 Q5tadd = ~Q5tmask + 1;
	uint32 Q5value = Q[offset+5];
	//uint32 Q4tmask = (Q4m5tunnel|Q4m4tunnel) & Q[offset+5];
	uint32 Q4tmask = Q4m5tunnel;
//	uint32 Q4tadd = ~Q4tmask + 1;
	uint32 Q4value = Q[offset+4];
	uint32 Q14tmask = Q14m6Q3m5tunnel & ~(Q[offset+15]^Q[offset+16]);
//	uint32 Q14tadd = ~Q14tmask + 1;
	uint32 Q14value = Q[offset+14];
	uint32 Q3value = Q[offset+3];

	uint32 pT20 = md5_gg(Q[offset+20],Q[offset+19],Q[offset+18]) + Q[offset+17] + md5_ac[20];
	uint32 testQ21val = Qvalue[offset+21] ^ (Qprev[offset+21]&Q[offset+20]);
	uint32 Q14tcur = 0;
	do {
		Q14tcur -= 1; Q14tcur &= Q14tmask;
		Q[offset+14] = Q14tcur ^ Q14value;
		Q[offset+3] = Q14tcur ^ Q3value;

		uint32 Q5tcur = 0;
		do {
			Q5tcur -= 1; Q5tcur &= Q5tmask;
			Q[offset+5] = Q5tcur ^ Q5value;

			uint32 R5 = Q[offset+6]-Q[offset+5];
			uint32 pT5 = rotate_right(R5,md5_rc[5]) - md5_ac[5] - Q[offset+2];
			uint32 p2T20 = pT20 + pT5;

			uint32 Q4tcur = 0;
			do {
				Q[offset+4] = Q4tcur ^ Q4value;
				Q4tcur -= 1; Q4tcur &= Q4tmask;				

				uint32 T20 = p2T20 - md5_ff(Q[offset+5],Q[offset+4],Q[offset+3]);
				uint32 Q21 = Q[offset+20] + rotate_left(T20, 5);
				if (testQ21val != (Q21&Qvaluemask[offset+21]))
					continue;
				Q[offset+21] = Q21;
				TESTCOUNT(1);

				uint32 R20 = rotate_left(T20, 5);
				uint32 T20b = T20 + dT[offset+20];						
				uint32 R20b = rotate_left(T20b, md5_rc[20]);
				if (R20b-R20 != dR[offset+20]) 
					continue;
				TESTCOUNT(2);

				m[5] = pT5 - md5_ff(Q[offset+5],Q[offset+4],Q[offset+3]);

#if 0
				if (verifyconds() <= 21)
					cerr << verifyconds() << endl;
#else
				do_step21();
#endif

			} while (Q4tcur != 0);
		} while (Q5tcur != 0);
	} while (Q14tcur != 0);
}

void collisionfinding_thread::do_step18()
{
	UPDATE(18);

	uint32 Q8pmask = ~Qvaluemask[offset+8] & Qprev[offset+9];
//	uint32 Q8padd = ~Q8pmask + 1;
	uint32 Q8pcur = 0;
	uint32 Q9mask = ~Qvaluemask[offset+9];
//	uint32 Q9add = ~Q9mask + 1;
	uint32 Q9cur = 0;
	uint32 Q10mask = ~Qvaluemask[offset+10];
//	uint32 Q10add = ~Q10mask +1;
	uint32 Q10cur = 0;
	uint32 Q11mask = ~Qvaluemask[offset+11];
//	uint32 Q11add = ~Q11mask +1;
	uint32 Q11cur = 0;
	uint32 Q12mask = ~Qvaluemask[offset+12] & ~Qprev[offset+13];
//	uint32 Q12add = ~Q12mask + 1;
	uint32 Q12cur = 0;
	uint32 Q8mask = (~Qvaluemask[offset+8]) & (~Qprev[offset+9]);
//	uint32 Q8add = ~Q8mask + 1;
	uint32 Q8cur = 0;

	uint32 Q8corr = Q8mask & rotate_right(Qvaluemask[offset+19],14);
	uint32 Q8mask2 = Q8mask & ~Q8corr;
//	uint32 Q8add2 = ~Q8mask2 + 1;
	uint32 Q19uncorr = Qvaluemask[offset+19] & ~rotate_left(Q8corr,14);

	uint32 pT18 = md5_gg(Q[offset+18],Q[offset+17],Q[offset+16]) + Q[offset+15] + md5_ac[18];
	uint32 pT19 = Q[offset+16] + md5_ac[19];
	uint32 Q19testval = Qvalue[offset+19] ^ (Qprev[offset+19]&Q[offset+18]);
	unsigned testcount = 0;	

	uint32 Q8pvalue = Qvalue[offset+8] ^ (Qprev[offset+8]&Q[offset+7])
		^ (xrng128()&Q8pmask);
	Q8pcur = 0;
	do {
		Q8pcur -= 1; Q8pcur &= Q8pmask;
		Q[offset+8] = Q8pcur ^ Q8pvalue;

		uint32 Q9value = Qvalue[offset+9] ^ (Qprev[offset+9]&Q[offset+8])
			^ (xrng128()&Q9mask);
		Q9cur = 0;
		do {
			Q9cur -= 1; Q9cur &= Q9mask;
			Q[offset+9] = Q9cur ^ Q9value;

			uint32 Q10value = Qvalue[offset+10] ^ (Qprev[offset+10]&Q[offset+9])
				^ (xrng128()&Q10mask);
			Q10cur = 0;
			do {
				Q10cur -= 1; Q10cur &= Q10mask;
				Q[offset+10] = Q10cur ^ Q10value;

				uint32 Q11value = Qvalue[offset+11] ^ (Qprev[offset+11]&Q[offset+10])
					^ (xrng128()&Q11mask);
				Q11cur = 0;
				do {
					Q11cur -= 1; Q11cur &= Q11mask;
					Q[offset+11] = Q11cur ^ Q11value;

					uint32 pT11 = md5_ff(Q[offset+11],Q[offset+10],Q[offset+9])
							+ md5_ac[11];

					uint32 Q12value = Qvalue[offset+12] 
							^ (Qprev[offset+12]&Q[offset+11]) 
							^ (Qprev[offset+13]&(Q[offset+13]^Qvalue[offset+13]))
							^ (xrng128()&Q12mask);
					Q12cur = 0;
					do {
						Q12cur -= 1; Q12cur &= Q12mask;
						Q[offset+12] = Q12cur ^ Q12value;

						uint32 R11 = Q[offset+12]-Q[offset+11];
						uint32 T11 = rotate_right(R11, md5_rc[11]) - pT11;						
						uint32 test = pT18 + T11;

						uint32 Q8value = Q[offset+8] ^ (rotate_right(Q19testval,14)&Q8corr);
						Q8cur = 0;

						/* fast rough test to see if we can correct bad bits */
/*						uint32 testT18 = test - Q8value;
						uint32 testQ19 = Q[offset+18] + rotate_left(testT18, 14);
						uint32 badmask = Q19testval ^ (testQ19&Qvaluemask[offset+19]);
						uint32 badmask2 = rotate_right(badmask, 14);
						if ((badmask2 & (Q8mask|(Q8mask>>1))) != badmask2)
						{
							if (testcount == 0)
								return;
							continue;
						}
*/						uint32 testT18 = test - ((xrng64() & Q8mask2) ^ Q8value);
						uint32 testQ19 = Q[offset+18] + rotate_left(testT18,14);
						if ((testQ19^Q19testval)&Q19uncorr) {
							if (testcount == 0)
								return;
							continue;
						}

//UPDATE(19);
/*** CRITICAL PART ***/
						do {
							Q8cur -= 1; Q8cur &= Q8mask2;

							uint32 T18 = test - (Q8cur ^ Q8value);
							uint32 Q19 = Q[offset+18] + rotate_left(T18, 14);
							if ((Q19^Q19testval)&Q19uncorr) 
								continue;
							uint32 corr = rotate_right(Q19^Q19testval,14) & Q8corr;
							T18 = test - (Q8cur ^ Q8value ^ corr);
							Q19 = Q[offset+18] + rotate_left(T18, 14);
							if (Q19testval != (Q19&Qvaluemask[offset+19]))
								continue;
//UPDATE(20);							

							UPDATE(19);
							uint32 p2T19 = pT19 + md5_gg(Q19, Q[offset+18], Q[offset+17]);
							uint32 Q20testval = Qvalue[offset+20] ^ (Qprev[offset+20]&Q19);
							for (Q1Q2m0it = Q1Q2m0withm1ok.begin(); Q1Q2m0it != Q1Q2m0withm1ok.end(); ++Q1Q2m0it)
							{
								uint32 T19 = p2T19 + Q1Q2m0it->third;
								uint32 Q20 = Q19 + rotate_left(T19, 20);
								if (Q20testval != (Q20&Qvaluemask[offset+20]))
									continue;
//UPDATE(21);
/*** END CRITICAL PART ***/
								Q[offset+20] = Q20;

								uint32 R18 = rotate_left(T18, 14);
								uint32 T18b = T18 + dT[offset+18];
								uint32 R18b = rotate_left(T18b, md5_rc[18]);
								if (R18b-R18 != dR[offset+18]) 
									continue;

								uint32 R19 = rotate_left(T19, md5_rc[19]);
								uint32 T19b = T19 + dT[offset+19];						
								uint32 R19b = rotate_left(T19b, md5_rc[19]);
								if (R19b-R19 != dR[offset+19]) 
									continue;

								Q[offset+19] = Q19;
								Q[offset+8] = Q8cur ^ Q8value ^ corr;
								m[11] = T11 - Q[offset+8];

								Q[offset+1] = Q1Q2m0it->first;
								Q[offset+2] = Q1Q2m0it->second;
								m[0] = Q1Q2m0it->third;
#if 0
								if (verifyconds() <= 20)
									cerr << verifyconds() << endl;
#else
								do_step20();
#endif
								++testcount;
							}
						} while (Q8cur != 0);
						if (testcount == 0)
							return;
					} while (Q12cur != 0);
				} while (Q11cur != 0);
			} while (Q10cur != 0);
		} while (Q9cur != 0);
	} while (Q8pcur != 0);
}


void collisionfinding_thread::do_step17()
{
	UPDATE(17);

	uint32 Q3mask = ~Qvaluemask[offset+3];
//	uint32 Q3add = Qvaluemask[offset+3]+1;
	uint32 Q3cur = 0;
	uint32 Q4mask = ~Qvaluemask[offset+4];
//	uint32 Q4add = Qvaluemask[offset+4]+1;
	uint32 Q4cur = 0;
	uint32 Q5mask = ~Qvaluemask[offset+5];
//	uint32 Q5add = Qvaluemask[offset+5]+1;
	uint32 Q5cur = 0;
	uint32 Q6mask = ~Qvaluemask[offset+6];
//	uint32 Q6add = Qvaluemask[offset+6]+1;
	uint32 Q6cur = 0;
	uint32 Q7mask = ~Qvaluemask[offset+7];
//	uint32 Q7add = Qvaluemask[offset+7]+1;
	uint32 Q7cur = 0;

	uint32 pT17 = md5_gg(Q[offset+17],Q[offset+16],Q[offset+15]) + Q[offset+14] + md5_ac[17];
	uint32 Q3value = Qvalue[offset+3] ^ (Qprev[offset+3]&Q[offset+2]);

	/* do not iterate over all values of Q3 */
	Q3value ^= xrng128() & Q3mask;
	Q3mask = 0;
	/**/

	do {
		Q3cur -= 1; Q3cur &= Q3mask;
		Q[offset+3] = Q3cur ^ Q3value;

		uint32 Q4value = Qvalue[offset+4] ^ (Qprev[offset+4]&Q[offset+3]);
		Q4cur = 0;
		do {
			Q4cur -= 1; Q4cur &= Q4mask;
			Q[offset+4] = Q4cur ^ Q4value;

			uint32 Q5value = Qvalue[offset+5] ^ (Qprev[offset+5]&Q[offset+4]);
			Q5cur = 0;
			do {
				Q5cur -= 1; Q5cur &= Q5mask;
				Q[offset+5] = Q5cur ^ Q5value;

				uint32 Q6value = Qvalue[offset+6] ^ (Qprev[offset+6]&Q[offset+5]);
				Q6cur = 0;
				do {
					Q6cur -= 1; Q6cur &= Q6mask;
					Q[offset+6] = Q6cur ^ Q6value;

					uint32 pT6 = md5_ff(Q[offset+6],Q[offset+5],Q[offset+4])
							+ Q[offset+3] + md5_ac[6];
					uint32 Q7value = Qvalue[offset+7] ^ (Qprev[offset+7]&Q[offset+6]);

					/* fast rough test to see if we can correct bad bits */
					uint32 testR6 = Q7value - Q[offset+6];
					uint32 testT6 = rotate_right(testR6, md5_rc[6]);
					uint32 testT17 = pT17 + testT6 - pT6;
					uint32 testR17 = rotate_left(testT17, md5_rc[17]);
					uint32 testQ18 = Q[offset+17] + testR17;
					uint32 badmask = Qvalue[offset+18] ^
							((testQ18&Qvaluemask[offset+18])
								^(Qprev[offset+18]&Q[offset+17]));
					uint32 badmask2 = rotate_left(badmask, 8);
					if ((badmask2 & (Q7mask|(Q7mask>>1))) != badmask2)
						continue;
					
					Q7cur = 0;
					do {
						Q7cur -= 1; Q7cur &= Q7mask;
						Q[offset+7] = Q7cur ^ Q7value;

						uint32 R6 = Q[offset+7]-Q[offset+6];
						uint32 T6 = rotate_right(R6, md5_rc[6]);
						m[6] = T6 - pT6;
						uint32 T17 = pT17 + m[6];
						uint32 R17 = rotate_left(T17, md5_rc[17]);
						Q[offset+18] = Q[offset+17] + R17;
						if (Qvalue[offset+18] != 
							((Q[offset+18]&Qvaluemask[offset+18])
								^(Qprev[offset+18]&Q[offset+17])))
							continue;

						uint32 T17b = T17 + dT[offset+17];						
						uint32 R17b = rotate_left(T17b, md5_rc[17]);
						if (R17b-R17 != dR[offset+17]) 
							continue;

						do_step18();

					} while (Q7cur != 0);
				} while (Q6cur != 0);
			} while (Q5cur != 0);
		} while (Q4cur != 0);
	} while (Q3cur != 0);
}

void collisionfinding_thread::do_step16()
{
	UPDATE(16);
	while (true) {
		Q[offset+1] = ((xrng128() & ~Qvaluemask[offset+1]) | Qvalue[offset+1]) ^ (Qprev[offset+1]&Q[offset+0]);
		if (!checkrotation(Q[offset+0], Q[offset+1], dT[offset+0], dR[offset+0], md5_rc[0]))
			continue;
		Q[offset+2] = ((xrng128() & ~Qvaluemask[offset+2]) | Qvalue[offset+2]) ^ (Qprev[offset+2]&Q[offset+1]);
		if (!checkrotation(Q[offset+1], Q[offset+2], dT[offset+1], dR[offset+1], md5_rc[1]))
			continue;
		// determine m0
		uint32 R0 = Q[offset+1] - Q[offset+0];
		uint32 T0 = rotate_right(R0, md5_rc[0]);
		uint32 F0 = md5_ff(Q[offset+0],Q[offset-1],Q[offset-2]);
		m[0] = T0 - F0 - Q[offset-3] - md5_ac[0];
		// determine m1
		uint32 R1 = Q[offset+2] - Q[offset+1];
		uint32 T1 = rotate_right(R1, md5_rc[1]);
		uint32 F1 = md5_ff(Q[offset+1],Q[offset+0],Q[offset-1]);
		m[1] = T1 - F1 - Q[offset-2] - md5_ac[1];

		Q[offset+13] = (xrng128() & (~Qvaluemask[offset+13] | Qprev[offset+13])) | Qvalue[offset+13];
		Q[offset+14] = ((xrng128() & ~Qvaluemask[offset+14]) | Qvalue[offset+14]) ^ (Qprev[offset+14]&Q[offset+13]);
		Q[offset+15] = ((xrng128() & ~Qvaluemask[offset+15]) | Qvalue[offset+15]) ^ (Qprev[offset+15]&Q[offset+14]);
		for (unsigned k = 0; k < 128; ++k)
		{
			Q[offset+16] = ((xrng128() & ~Qvaluemask[offset+16]) | Qvalue[offset+16]) ^ (Qprev[offset+16]&Q[offset+15]);
			if (!checkrotation(Q[offset+15], Q[offset+16], dT[offset+15], dR[offset+15], md5_rc[15]))
				continue;
			uint32 F16 = md5_gg(Q[offset+16],Q[offset+15],Q[offset+14]);
			uint32 T16 = F16 + m[1] + Q[offset+13] + md5_ac[16];
			uint32 R16 = rotate_left(T16, md5_rc[16]);
			Q[offset+17] = Q[offset+16] + R16;
			if (Qvalue[offset+17] != ((Q[offset+17]&Qvaluemask[offset+17])^(Qprev[offset+17]&Q[offset+16])))
				continue;
			uint32 T16b = T16 + dT[offset+16];
			uint32 R16b = rotate_left(T16b, md5_rc[16]);
			if (R16b - R16 != dR[offset+16]) 
				continue;

			Q1Q2m0withm1ok.clear();
			uint32 Q1mask = ~Qvaluemask[offset+1];
//			uint32 Q1add = ~Q1mask + 1;
			uint32 Q1val = Q[offset+1];
			uint32 Q1cur = 0;
			do {
				uint32 Q1 = Q1cur ^ Q1val;
				Q1cur -= 1; Q1cur &= Q1mask;

				uint32 F1 = md5_ff(Q1, Q[offset+0], Q[offset-1]);
				uint32 T1 = F1 + Q[offset-2] + md5_ac[1] + m[1];
				uint32 R1 = rotate_left(T1, md5_rc[1]);
				uint32 Q2 = Q1 + R1;
				if ((Q2^Q[offset+2])&Qprev[offset+3])
					continue;
				if (Qvalue[offset+2] != ((Q2&Qvaluemask[offset+2])^(Qprev[offset+2]&Q1)))
					continue;
				if (!checkrotation(Q[offset+0], Q1, dT[offset+0], dR[offset+0], md5_rc[0]))
					continue;
				if (!checkrotation(Q1, Q2, dT[offset+1], dR[offset+1], md5_rc[1]))
					continue;
				uint32 m0 = rotate_right(Q1 - Q[offset+0], md5_rc[0]) 
							- md5_ff(Q[offset+0],Q[offset-1],Q[offset-2])
							- Q[offset-3] 
							- md5_ac[0];
				
				Q1Q2m0withm1ok.push_back( make_triple(Q1,Q2,m0) );
			} while (Q1cur != 0 && Q1Q2m0withm1ok.size() < (1<<20));
//			cout << "Q1Q2m0withm1ok: " << Q1Q2m0withm1ok.size() << endl;

			if (isinfinite) { cout << "." << flush; isinfinite = false; }
			do_step17();
			return;
		}
	}	
}

void collisionfinding_thread::filltables()
{
	isinfinite = true;
	for (int t = -3; t <= 64; ++t)
	{
		dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
		Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
		Qprev[offset+t] = Qprev2[offset+t] = 0;
	}
	Q5m5tunnel=0; Q4m5tunnel=0; Q4m4tunnel=0;
	Q10m10tunnel=0; Q9m10tunnel=0; Q9m9tunnel=0;
	Q8Q12m15tunnel=0; Q14m6Q3m5tunnel=0; Q14Q3m14tunnel=0;

	// determine tunnels and update differential path
	// Q14Q3m14 tunnel
	if (1)
	for (unsigned b = 0; b < 32; ++b)
	{
		bitcondition Q3b = diffpath(3,b), Q4b = diffpath(4,b), Q5b = diffpath(5,b);
		bitcondition Q14b = diffpath(14,b), Q15b = diffpath(15,b), Q16b = diffpath(16,b);
		if (Q3b == bc_constant && Q14b == bc_constant
			&& (Q4b == bc_constant || Q4b == bc_zero || Q4b == bc_plus)
			&& (Q5b == bc_constant || Q5b == bc_one || Q5b == bc_minus || Q5b == bc_prevn)
			&& (Q15b == bc_constant || Q15b == bc_zero || Q15b == bc_plus)
			&& (Q16b == bc_constant || Q16b == bc_constant || Q16b == bc_plus || Q16b == bc_prev))
		{
			Q14Q3m14tunnel |= 1<<b;
			diffpath.setbitcondition(3,b,bc_zero);
			diffpath.setbitcondition(14,b,bc_zero);
			if (Q4b == bc_constant)
				diffpath.setbitcondition(4,b,bc_zero);
			if (Q5b == bc_constant || Q5b == bc_prevn)
				diffpath.setbitcondition(5,b,bc_one);
			if (Q15b == bc_constant)
				diffpath.setbitcondition(15,b,bc_zero);
			if (Q16b == bc_constant || Q16b == bc_prev)
				diffpath.setbitcondition(16,b,bc_zero);
		}
	}

	// Q8Q12m15 tunnel
	for (unsigned b = 0; b < 32; ++b)
	{
		bitcondition Q8b = diffpath(8,b), Q9b = diffpath(9,b), Q10b = diffpath(10,b);
		bitcondition Q12b = diffpath(12,(b+22)&31), Q13b = diffpath(13,(b+22)&31);
		if (Q8b == bc_constant && Q12b == bc_constant
				&& Q9b != bc_prev && Q9b != bc_prevn
				&& Q13b != bc_prev && Q13b != bc_prevn
				&& (Q10b == bc_constant || Q10b == bc_minus || Q10b == bc_one)
			)
		{
			Q8Q12m15tunnel |= 1<<b;
			diffpath.setbitcondition(8,b,bc_zero);
			diffpath.setbitcondition(12,(b+22)&31,bc_zero);
			if (Q10b == bc_constant)
				diffpath.setbitcondition(10,b,bc_one);
		}
	}
	// Q14m6Q3m5 tunnel
	for (unsigned b = 0; b < 32; ++b)
	{
		bitcondition Q14b = diffpath(14,b), Q3b = diffpath(3,b)
					, Q15b = diffpath(15,b), Q16b = diffpath(16,b)
					, Q4b = diffpath(4,b);
		if (Q14b == bc_constant && Q3b == bc_constant
				&& Q15b != bc_prev && Q15b != bc_prevn
				&& Q4b != bc_prev && Q4b != bc_prevn)
		{
			bool istunnel = false;
			if (Q16b == bc_prev)
				istunnel = true;
			else if (Q16b == bc_constant) {
				istunnel = true;
				if (Q15b == bc_constant)
					diffpath.setbitcondition(16,b,bc_prev);
				else if (Q15b == bc_zero || Q15b == bc_plus)
					diffpath.setbitcondition(16,b,bc_zero);
				else if (Q15b == bc_one || Q15b == bc_minus)
					diffpath.setbitcondition(16,b,bc_one);
				else
					istunnel = false;
			} else if (Q16b == bc_zero || Q16b == bc_plus) {
				if (Q15b == bc_constant)
					diffpath.setbitcondition(15,b,bc_zero);
				if (Q15b == bc_zero || Q15b == bc_plus)
					istunnel = true;
			} else if (Q16b == bc_one || Q16b == bc_minus) {
				if (Q15b == bc_constant)
					diffpath.setbitcondition(15,b,bc_one);
				if (Q15b == bc_one || Q15b == bc_minus)
					istunnel = true;
			}
			if (istunnel) {
				diffpath.setbitcondition(14,b,bc_zero);
				diffpath.setbitcondition(3,b,bc_zero);
				Q14m6Q3m5tunnel |= 1<<b;
			}
		}
	}
	// Q4m4, Q4m5, Q5m5 tunnels
	for (unsigned b = 0; b < 32; ++b)
	{
		bitcondition Q4b = diffpath(4,b), Q5b = diffpath(5,b), Q6b = diffpath(6,b);
		if ((Q4b == bc_constant) 
			&& (Q5b == bc_constant || Q5b == bc_zero || Q5b == bc_plus)
			&& (Q6b == bc_constant || Q6b == bc_one || Q6b == bc_minus || Q6b == bc_prevn))
		{
			// if possible we maintain a dynamic tunnel on Q4m4 Q4m5
			if (Q5b == bc_constant)// && Q6b == bc_prevn)
				diffpath.setbitcondition(5,b,bc_zero);
			if (Q6b == bc_constant || Q6b == bc_prevn)
				diffpath.setbitcondition(6,b,bc_one);
			diffpath.setbitcondition(4,b,bc_zero);
			Q4m4tunnel |= 1<<b;
		} 
		else if ((Q4b == bc_constant)
			&& (Q5b == bc_constant || Q5b == bc_one || Q5b == bc_minus)
			&& (Q6b == bc_constant || Q6b == bc_one || Q6b == bc_minus || Q6b == bc_prev))
		{
			// if possible we maintain a dynamic tunnel on Q4m4 Q4m5
			if (Q5b == bc_constant)// && Q6b == bc_prev) 
				diffpath.setbitcondition(5,b,bc_one);
			if (Q6b == bc_constant || Q6b == bc_prev)
				diffpath.setbitcondition(6,b,bc_one);
			diffpath.setbitcondition(4,b,bc_zero);
			Q4m5tunnel |= 1<<b;
		} 
		else if ((Q5b == bc_constant)
			&& (Q6b == bc_constant || Q6b == bc_zero || Q6b == bc_plus))
		{
			if (Q6b == bc_constant)
				diffpath.setbitcondition(6,b,bc_zero);
			diffpath.setbitcondition(5,b,bc_zero);
			Q5m5tunnel |= 1<<b;			
		}
	}
	// Q9m9, Q9m10, Q10m10 tunnels
	for (unsigned b = 0; b < 32; ++b)
	{
		bitcondition Q9b = diffpath(9,b), Q10b = diffpath(10,b), Q11b = diffpath(11,b);
		if ((Q9b == bc_constant) 
			&& (Q10b == bc_constant || Q10b == bc_zero || Q10b == bc_plus)
			&& (Q11b == bc_constant || Q11b == bc_one || Q11b == bc_minus || Q11b == bc_prevn))
		{
			// if possible we maintain a dynamic tunnel on Q9m9 Q9m10
			if (Q10b == bc_constant && Q11b == bc_prevn) 
				diffpath.setbitcondition(10,b,bc_zero);

			if (Q11b == bc_constant || Q11b == bc_prevn)
				diffpath.setbitcondition(11,b,bc_one);
			diffpath.setbitcondition(9,b,bc_zero);
			Q9m9tunnel |= 1<<b;			
		} 
		else if ((Q9b == bc_constant)
			&& (Q10b == bc_constant || Q10b == bc_one || Q10b == bc_minus)
			&& (Q11b == bc_constant || Q11b == bc_one || Q11b == bc_minus || Q11b == bc_prev))
		{
			// if possible we maintain a dynamic tunnel on Q9m9 Q9m10
			if (Q10b == bc_constant && Q11b == bc_prev) 
				diffpath.setbitcondition(10,b,bc_one);

			if (Q11b == bc_constant || Q11b == bc_prev)
				diffpath.setbitcondition(11,b,bc_one);
			diffpath.setbitcondition(9,b,bc_zero);
			Q9m10tunnel |= 1<<b;			
		} 
		else if ((Q10b == bc_constant)
			&& (Q11b == bc_constant || Q11b == bc_zero || Q11b == bc_plus))
		{
			if (Q11b == bc_constant)
				diffpath.setbitcondition(11,b,bc_zero);
			diffpath.setbitcondition(10,b,bc_zero);
			Q10m10tunnel |= 1<<b;			
		}
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
	for (int t = 0; t < 64; ++t)
	{
		booleanfunction* F = 0;
		if (t < 16) F = & MD5_F_data;
		else if (t < 32) F = & MD5_G_data;
		else if (t < 48) F = & MD5_H_data;
		else F = & MD5_I_data;
		uint32 dF = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			bf_outcome outcome = F->outcome(diffpath(t,b), diffpath(t-1,b), diffpath(t-2,b));
			if (outcome.size()) {
				if (outcome[0] == bc_plus) 			dF += 1<<b;
				else if (outcome[0] == bc_minus)	dF -= 1<<b;
			} else throw std::runtime_error("ambiguous path!!");
		}
		dT[offset+t] = dQ[offset+t-3] + dF + m_diff[md5_wt[t]];
		dR[offset+t] = dQ[offset+t+1] - dQ[offset+t];
	}

	Q[0] = Qvalue[0]; Q2[0] = Q[0] + dQ[0];
	Q[1] = Qvalue[1]; Q2[1] = Q[1] + dQ[1];
	Q[2] = Qvalue[2]; Q2[2] = Q[2] + dQ[2];
	Q[3] = Qvalue[3]; Q2[3] = Q[3] + dQ[3];
	
	
	cout << "20: Q5m5tunnel      = " << hw(Q5m5tunnel) << endl;
	cout << "20: Q4m5tunnel      = " << hw(Q4m5tunnel) << endl;
	cout << "20: Q14m6Q3m5tunnel = " << hw(Q14m6Q3m5tunnel) << endl;
	cout << "21: Q10m10tunnel    = " << hw(Q10m10tunnel) << endl;
	cout << "21: Q9m10tunnel     = " << hw(Q9m10tunnel) << endl;
	cout << "22: Q8Q12m15tunnel  = " << hw(Q8Q12m15tunnel) << endl;
	cout << "23: Q4m4tunnel      = " << hw(Q4m4tunnel) << endl;
	cout << "24: Q9m9tunnel      = " << hw(Q9m9tunnel) << endl;
	cout << "25: Q14Q3m14tunnel  = " << hw(Q14Q3m14tunnel) << endl;

}
static const uint64_t Q7m12m13m14size = 1ULL<<28;



void collisionfinding_alphabet_thread::do_step23()
{
	if (++testcounts[223] == 1)
		std::cout << "t22 t23 " << std::flush;
	std::cout << "T1Q7 " << std::flush;
	// known: m[15], m[1], m[6], m[11], m[0], m[5], m[10], set of { m[14], m[13] }
	// satisfied: Q12 -- Q23
	// precompute list of { m[12], m[13], m[14] => Q7 } for matching in dostep23
	Q7m12m13m14.clear();
	Q7m12m13m14.reserve(Q7m12m13m14size);
	unsigned l1 = 0, l2 = 0, l3 = 0, y = 0;
	for (auto& x : m14m13)
	{
		if (++y > 2 && Q7m12m13m14.empty())
			break;
		m[14] = x.first;
		m[13] = x.second;
		computeQtm3(14);
		computeQtm3(13);
		
		auto m12range = MA.word_range(12);
		auto Q9mv = masked_value_QtQtp1(9);
//		std::cout << "." << std::flush; // << m12range.count() << " " << hammingweight(masked_value_QtQtp1(8).mask) << " " << std::flush;
		for (uint32_t m12 : m12range)
		{
			m[12] = m12;
			computeQtm3(12);
//			std::cout << std::hex << Qt(9) << std::dec << " " << std::flush;
			if (!Q9mv.check(Qt(9)))
				continue;
//			if (hammingweight(++l1) == 1) std::cout << "L1=" << l1 << std::flush;
			computeQtm3(11);
			if (!masked_value_QtQtp1(8).check(Qt(8)))
				continue;
			if (hammingweight(++l2) == 1) std::cout << "L2=" << l2 << std::flush;
			computeQtm3(10);
			if (!masked_value_QtQtp1(7).check(Qt(7)))
				continue;
			if (hammingweight(++l3) == 1) std::cout << "L3=" << l3 << std::flush;
				
			Q7m12m13m14.emplace_back(Qt(7), m[12], m[13], m[14]);
			std::cout << "+" << std::flush;
			throw;
			if (Q7m12m13m14.size() >= Q7m12m13m14size)
				break;
		}
		if (Q7m12m13m14.size() >= Q7m12m13m14size)
			break;
	}
	std::cout << Q7m12m13m14.size() << std::endl;
}

void collisionfinding_alphabet_thread::do_step21()
{
	if (++testcounts[221] == 1)
		std::cout << "t21 " << std::flush;
	auto m10range = MA.word_range(10);
	auto Q22mv = masked_value_QtQtm1(22);
	if (m10range.count() < Q22mv.count())
	{
		std::cout << "t21m " << std::flush;
//		unsigned x = 0;
		for (uint32_t m10 : m10range)
		{
//			if (hammingweight(++x) == 1)
//				std::cout << x << " " << std::flush;
			m[10] = m10;
			computeQtp1(21);
			if (!Q22mv.check(Qt(22)))
				continue;
			auto Q23mv = masked_value_QtQtm1(23);
			computeQtp1(22); // using known m[15]
			if (Q23mv.check(Qt(23)))
				do_step23();
		}
	} else
	{
		std::cout << "t21Q " << std::flush;
		for (uint32_t Q22 : Q22mv.range_rnd())
		{
			Qt(22) = Q22;
			computeWt(21);
			if (!MA.checkword(10, m[10]))
				continue;
			auto Q23mv = masked_value_QtQtm1(23);
			computeQtp1(22); // using known m[15]
			if (Q23mv.check(Qt(23)))
				do_step23();
		}
	}
}

void collisionfinding_alphabet_thread::do_step20()
{
	if (++testcounts[220] == 1)
		std::cout << "t20 " << std::flush;
	auto m5range = MA.word_range(5);
	auto Q21mv = masked_value_QtQtm1(21);
	if (m5range.count() < Q21mv.count())
	{
		std::cout << "t20m " << std::flush;
//		unsigned x = 0;
		for (uint32_t m5 : m5range)
		{
//			if (hammingweight(++x) == 1)
//				std::cout << x << " " << std::flush;
			m[5] = m5;
			computeQtp1(20);
			if (!Q21mv.check(Qt(21)))
				continue;
			do_step21();
		}
	} else
	{
		std::cout << "t20Q " << std::flush;
		for (uint32_t Q21 : Q21mv.range_rnd())
		{
			Qt(21) = Q21;
			computeWt(20);
			if (!MA.checkword(5, m[5]))
				continue;
			do_step21();
		}
	}
}

void collisionfinding_alphabet_thread::do_step19()
{
	if (++testcounts[219] == 1)
		std::cout << "t19 " << std::flush;
	auto m0range = MA.word_range(0);
	auto Q20mv = masked_value_QtQtm1(20);
	if (m0range.count() < Q20mv.count())
	{
//		std::cout << "t19m " << std::flush;
		unsigned x = 0;
		for (uint32_t m0 : m0range)
		{
//			if (hammingweight(++x) == 1)
//				std::cout << x << " " << std::flush;
			m[0] = m0;
			computeQtp1(19);
			if (!Q20mv.check(Qt(20)))
				continue;
			do_step20();
		}
	} else
	{
		std::cout << "t19Q " << std::flush;
		for (uint32_t Q20 : Q20mv.range_rnd())
		{
			Qt(20) = Q20;
			computeWt(19);
			if (!MA.checkword(0, m[0]))
				continue;
			do_step20();
		}
	}
}

void collisionfinding_alphabet_thread::do_step18()
{
	if (++testcounts[218] == 1)
		std::cout << "t18 " << std::flush;
	auto m11range = MA.word_range(11);
	auto Q19mv = masked_value_QtQtm1(19);
	if (m11range.count() < Q19mv.count())
	{
		std::cout << "t18m " << std::flush;
//		unsigned x = 0;
		for (uint32_t m11 : m11range)
		{
//			if (hammingweight(++x) == 1)
//				std::cout << x << " " << std::flush;
			m[11] = m11;
			computeQtp1(18);
			if (!Q19mv.check(Qt(19)))
				continue;
			do_step19();
		}
	} else
	{
		std::cout << "t18Q " << std::flush;
		for (uint32_t Q19 : Q19mv.range_rnd())
		{
			Qt(19) = Q19;
			computeWt(18);
			if (!MA.checkword(11, m[11]))
				continue;
			do_step19();
		}
	}
}

void collisionfinding_alphabet_thread::do_step17()
{
//	if (++testcounts[217] == 1)
		std::cout << "t17 " << std::flush;
	auto m6range = MA.word_range(6);
	auto Q18mv = masked_value_QtQtm1(18);
	if (m6range.count() < Q18mv.count())
	{
		std::cout << "t17m " << std::flush;
//		unsigned x = 0;
		for (uint32_t m6 : m6range)
		{
//			if (hammingweight(++x) == 1)
//				std::cout << x << " " << std::flush;
			m[6] = m6;
			computeQtp1(17);
			if (!Q18mv.check(Qt(18)))
				continue;
			do_step18();
		}
	} else
	{
		std::cout << "t17Q " << std::flush;
		for (uint32_t Q18 : Q18mv.range_rnd())
		{
			Qt(18) = Q18;
			computeWt(17);
			if (!MA.checkword(6, m[6]))
				continue;
			do_step18();
		}
	}
}


void collisionfinding_alphabet_thread::extend_step_fw(int t, vec_state_t& out, vec_state_t& in, uint64_t N)
{
	int wt = md5_wt[t];
	out.clear(); out.reserve(N);
	
	std::cout << "\n======\nStep " << t << ": (Q" << (t+1) << ",m" << wt << "):";
	progress_display pd(N);
	
	auto wtrange = MA.word_range(wt);
	const size_t wtrangecount = wtrange.count(), Qtp1count = masked_value_QtQtm1(t+1).count();

	uint64_t attempts = 0;
	while (out.size() < N)
	{
		out.emplace_back();
		auto& cur = out.back();
		while (true)
		{
			++attempts;
			// sample random input value
			auto& incur = in[ xrng64() % in.size() ];
			// copy just the Q values needed for step t
			for (unsigned i = offset+t-3; i <= offset+t+1; ++i)
				Q[i] = incur[i];
			auto Qtp1mv = masked_value_QtQtm1(t+1);
			if (Qtp1count < wtrangecount)
			{
				Qt(t+1) = Qtp1mv.sample();
				computeWt(t);
				if (!MA.checkword(wt, m[wt]))
					continue;
			} else {
				m[wt] = MA.sampleword(wt);
				computeQtp1(t);
				if (!Qtp1mv.check(Qt(t+1)))
					continue;
			}
			if (!checkrotationQtQtp1(t))
				continue;
			// found solution
			cur = incur;
			cur[offset+t+1] = Qt(t+1);
			++pd;
			break;
		}
	}
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}

void collisionfinding_alphabet_thread::extend_step_bw(int t, vec_state_t& out, vec_state_t& in, uint64_t N)
{
	int wt = md5_wt[t];
	out.clear(); out.reserve(N);
	
	std::cout << "\n======\nStep " << t << ": (Q" << (t-3) << ",m" << wt << "):";
	progress_display pd(N);
	
	auto wtrange = MA.word_range(wt);
	const size_t wtrangecount = wtrange.count(), Qtm3count = masked_value_QtQtp1(t-3).count();

	uint64_t attempts = 0;
	while (out.size() < N)
	{
		out.emplace_back();
		auto& cur = out.back();
		while (true)
		{
			++attempts;
			// sample random input value
			auto& incur = in[ xrng64() % in.size() ];
			// copy just the Q values needed for step t
			for (unsigned i = offset+t-3; i <= offset+t+1; ++i)
				Q[i] = incur[i];
			auto Qtm3mv = masked_value_QtQtp1(t-3);
			if (Qtm3count < wtrangecount)
			{
				Qt(t-3) = Qtm3mv.sample();
				computeWt(t);
				if (!MA.checkword(wt, m[wt]))
					continue;
			} else {
				m[wt] = MA.sampleword(wt);
				computeQtm3(t);
				if (!Qtm3mv.check(Qt(t-3)))
					continue;
			}
			if (!checkrotationQtQtp1(t-3))
				continue;
			// found solution
			cur = incur;
			cur[offset+t-3] = Qt(t-3);
			++pd;
			break;
		}
	}
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}

void collisionfinding_alphabet_thread::extend_step_m11(vec_state_t& out, vec_state_t& in, uint64_t N)
{
	out.clear(); out.reserve(N);
	
	std::cout << "\n======\nStep: (m11,Q8,Q19):";
	progress_display pd(N);
	
	uint64_t attempts = 0;
	while (out.size() < N)
	{
		out.emplace_back();
		auto& cur = out.back();
		while (true)
		{
			++attempts;
			// sample random input value
			auto& incur = in[ xrng64() % in.size() ];
			// copy just the Q values needed for step t
			for (unsigned i = offset+8; i <= offset+19; ++i)
				Q[i] = incur[i];
			auto Q8mv = masked_value_QtQtp1(8);
			auto Q19mv = masked_value_QtQtm1(19);
			m[11] = MA.sampleword(11);
			computeQtp1(18);
			if (!Q19mv.check(Qt(19)))
				continue;
			computeQtm3(11);
			if (!Q8mv.check(Qt(8)))
				continue;
			if (!checkrotationQtQtp1(8))
				continue;
			if (!checkrotationQtQtp1(18))
				continue;
			// found solution
			cur = incur;
			cur[offset+8] = Qt(8);
			cur[offset+19] = Qt(19);
			++pd;
			break;
		}
	}
	std::cout << "Step m11: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}

void collisionfinding_alphabet_thread::extend_step_m10(vec_state_t& out, vec_state_t& in, uint64_t N)
{
	out.clear(); out.reserve(N);
	
	std::cout << "\n======\nStep: (m10,Q7,Q22):";
	progress_display pd(N);
	
	uint64_t attempts = 0;
	while (out.size() < N)
	{
		out.emplace_back();
		auto& cur = out.back();
		while (true)
		{
			++attempts;
			// sample random input value
			auto& incur = in[ xrng64() % in.size() ];
			// copy just the Q values needed for step t
			for (unsigned i = offset+7; i <= offset+22; ++i)
				Q[i] = incur[i];
			computeWt(15);
			
			auto Q7mv = masked_value_QtQtp1(7);
			auto Q22mv = masked_value_QtQtm1(22);
			m[10] = MA.sampleword(10);
			computeQtm3(10);
			if (!Q7mv.check(Qt(7)))
				continue;
			computeQtp1(21); // m10
			if (!Q22mv.check(Qt(22)))
				continue;
			computeQtp1(22); // m15
			if (!Q22mv.check(Qt(23)))
				continue;
			if (!checkrotationQtQtp1(7))
				continue;
			if (!checkrotationQtQtp1(21))
				continue;
			if (!checkrotationQtQtp1(22))
				continue;
			// found solution
			cur = incur;
			cur[offset+7] = Qt(7);
			cur[offset+22] = Qt(22);
			cur[offset+23] = Qt(23);
			++pd;
			break;
		}
	}
	std::cout << "Step m10: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}

void collisionfinding_alphabet_thread::do_step16_new()
{
	static const uint64_t N = 1ULL<<20;
	vector< array<uint32_t, 32> > in, out;

	std::cout << "Trying to load 'out.bin.gz'..." << std::flush;
	try {
		load_gz(out, "out", binary_archive);
		std::cout << "done." << out.size() << std::endl;
	} catch (std::exception& e)
	{
		out.clear();
		std::cout << "failed." << std::endl;
	}

if (out.empty())
{

	in.reserve(N);
	out.reserve(N);
	
	std::cout << "\n======\nGenerate Q13-Q17,m1:";
	progress_display genpd(N);
	uint64_t attempts = 0;
	auto m1range = MA.word_range(1);
	const size_t m1rangecount = m1range.count(), Q17count = masked_value_QtQtm1(17).count();
	while (out.size() < N)
	{
		++attempts;
		Qt(13) = masked_value_Qt(13).sample();
		Qt(14) = masked_value_QtQtm1(14).sample();
		Qt(15) = masked_value_QtQtm1(15).sample();
		Qt(16) = masked_value_QtQtm1(16).sample();
		
		auto Q17 = masked_value_QtQtm1(17);
		if (Q17count < m1rangecount)
		{
			Qt(17) = Q17.sample();
			computeWt(16);
			if (!MA.checkword(1, m[1]))
				continue;
		} else {
			m[1] = MA.sampleword(1);
			computeQtp1(16);
			if (!Q17.check(Qt(17)))
				continue;
		}
		++genpd;
		out.emplace_back();
		auto& cur = out.back();
		for (int t = 13; t <= 17; ++t)
			cur[offset+t] = Q[offset+t];
	}
	std::cout << "Generate Q13-Q17,m1: attempts: 2^" << log(double(attempts))/log(2.0) << " success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
	
	// m12 - m15
	for (int t = 15; t >= 12; --t)
	{
		std::swap(in, out); out.clear();
		extend_step_bw(t, out, in, N);
	}
	std::swap(in, out); out.clear();
	extend_step_fw(17, out, in, N); // m6
	std::swap(in, out); out.clear();
	extend_step_m11(out, in, N);    // m11
	std::swap(in, out); out.clear();
	extend_step_fw(19, out, in, N); // m0
	std::swap(in, out); out.clear();
	extend_step_fw(20, out, in, N); // m5
	std::swap(in, out); out.clear();
	extend_step_m10(out, in, N); // m10

	save_gz(out, "out", binary_archive);
} // if out not loaded successfully

	// sort on maximum Q9m9 tunnel	
	const uint32_t Q9m9 = (~Qvaluemask[offset+9] & ~Qprev[offset+10]);
	std::sort(out.begin(), out.end(), 
		[Q9m9](const state_t& l, const state_t& r) 
		{
			uint32_t lQ9m9 = Q9m9 & l[offset+11] & ~l[offset+10];
			uint32_t rQ9m9 = Q9m9 & r[offset+11] & ~r[offset+10];
			return hammingweight(lQ9m9) > hammingweight(rQ9m9);
		});
	std::cout << "Q9m9 tunnel best state: " << hammingweight(out.front()[offset+11] & ~(out.front()[offset+10]) & Q9m9) << std::endl;

	std::map< std::array<uint32_t, 3>, std::array< uint32_t, 2 > > Q7810m1213;
	std::unordered_multimap< uint32_t, std::array< uint32_t, 3 > > Q7m10m12m13;
	std::vector<uint32_t> m10good;
	for (auto& cur : out)
	{
		for (unsigned i = offset+7; i <= offset+23; ++i)
			Q[i] = cur[i];
		for (int t = 10; t < 16; ++t)
		{
			computeWt(t);
			if (!MA.checkword(t, m[t]))
				throw std::runtime_error("mt-alphabet violation");
		}
		const uint32_t m10org = m[10], m15org = m[15];
		for (int t = 16; t <= 22; ++t)
		{
			computeWt(t);
			if (!MA.checkword(md5_wt[t], m[md5_wt[t]]))
				throw std::runtime_error("mt-alphabet violation");
		}
		if (m[10] != m10org)
			throw std::runtime_error("m10-inconsistency");
		if (m[15] != m15org)
			throw std::runtime_error("m15-inconsistency");

		uint32_t m13org = m[13], Q10org = Qt(10);
		std::vector< std::pair<uint32_t,uint32_t> > m13good;
		for (uint32_t m13 : MA.word_range(13))
		{
			m[13] = m13;
			computeQtm3(13);
			if (!masked_value_QtQtp1(10).check(Qt(10)))
				continue;
			if (!checkrotationQtQtp1(10))
				continue;
			m13good.emplace_back(m[13], Qt(10));
		}
		std::cout << "m13goodsize: " << m13good.size() << std::endl;
		uint32 Q11 = Qt(11);
		std::sort(m13good.begin(), m13good.end(), 
			[Q11,Q9m9](auto& l, auto& r)
			{
				return hammingweight(Q9m9 & Q11 & ~l.second)
					> hammingweight(Q9m9 & Q11 & ~r.second);
			});
		
		// find m12 that still satisfies Q9 and Q8 conditions
		// forget m12 for which there is already another m12 with same Q8 and Q7
		// i.e. the difference is part of the Q9m9 tunnel
		progress_display pdQ7810(m13good.size());
		Q7810m1213.clear();
		auto m12range = MA.word_range(12);
		for (auto& m13Q10 : m13good)
		{
			++pdQ7810;
		for (uint32_t m12 : m12range)
		{
			m[13] = m13Q10.first;
			m[12] = m12;
			m[10] = m10org;
			computeQtm3(13);
			if (!masked_value_QtQtp1(10).check(Qt(10)))
				continue;
			computeQtm3(12);
			if (!masked_value_QtQtp1(9).check(Qt(9)))
				continue;
			computeQtm3(11);
			if (!masked_value_QtQtp1(8).check(Qt(8)))
				continue;
			computeQtm3(10);
			if (!checkrotationQtQtp1(10))
				continue;
			if (!checkrotationQtQtp1(9))
				continue;
			if (!checkrotationQtQtp1(8))
				continue;
			computeQtm3(10); // virtual Q7-under-m10org
			Q7810m1213[ {Qt(7),Qt(8),Qt(10)} ] = { m[12], m[13] };
			if (Q7810m1213.size() >= 1ULL<<18)
				break;
		}
			if (Q7810m1213.size() >= 1ULL<<18)
				break;
		}
		std::cout << "\nQ7810m1213size: " << Q7810m1213.size() << std::endl;

		// find all m10 that satisfy Q22-Q23
		m[4] = MA.sampleword(4);
		m10good.clear();
		auto m10range = MA.word_range(10);
		auto Q22mv = masked_value_QtQtm1(22);
		for (uint32_t m10 : m10range)
		{
			m[10] = m10;
			computeQtp1(21); // m10
			if (!Q22mv.check(Qt(22)))
				continue;
			computeQtp1(22); // m15
			if (!masked_value_QtQtm1(23).check(Qt(23)))
				continue;
			computeQtp1(23); // m4
			if (!masked_value_QtQtm1(24).check(Qt(24)))
				continue;
			if (!checkrotationQtQtp1(21))
				continue;
			if (!checkrotationQtQtp1(22))
				continue;
			if (!checkrotationQtQtp1(23))
				continue;
			m10good.emplace_back(m10);
		}
		std::cout << "m10goodsize: " << m10good.size() << std::endl;
		
		// iterate over those m12, and find all (m10, m12) satisfying Q7-Q9, Q22-Q23
		progress_display pdm12m13(Q7810m1213.size());
		Q7m10m12m13.clear();
		auto Q7mv = masked_value_QtQtp1(7);
		uint64_t Q7lookupstore = 0;
		for (auto& Qm12it : Q7810m1213)
		{
			++pdm12m13;
			m[12] = Qm12it.second[0];
			m[13] = Qm12it.second[1];
			computeQtm3(13);
			computeQtm3(12);
			computeQtm3(11);
			for (uint32_t m10 : m10good)
			{
				m[10] = m10;
				computeQtm3(10); // m10
				if (!Q7mv.check(Qt(7)))
					continue;
				computeQtp1(21); // m10
				if (!Q22mv.check(Qt(22)))
					continue;
				computeQtp1(22); // m15
				if (!masked_value_QtQtm1(23).check(Qt(23)))
					continue;
				computeQtp1(23); // m4
				if (!masked_value_QtQtm1(24).check(Qt(24)))
					continue;
				if (!checkrotationQtQtp1(7))
					continue;
				if (!checkrotationQtQtp1(21))
					continue;
				if (!checkrotationQtQtp1(22))
					continue;
				if (!checkrotationQtQtp1(23))
					continue;
				++Q7lookupstore;
				Q7m10m12m13.emplace(Qt(7), std::array<uint32_t,3>({m[10],m[12],m[13]}));//std::pair<uint32_t,uint32_t>(m[10],m[12]);
			}
		}
		std::cout << "Q7m10m12m13size: " << Q7m10m12m13.size() << " (out of " << Q7lookupstore << " good kept)" << std::endl;
		
		Q[0] = md5_iv[0]; Q[1] = md5_iv[3]; Q[2] = md5_iv[2]; Q[3] = md5_iv[1];
		Q2[0] = md5_iv[0]; Q2[1] = md5_iv[3]; Q2[2] = md5_iv[2]; Q2[3] = md5_iv[1];
		computeQtp1(0); // m0 is fixed
		if (!masked_value_QtQtm1(1).check(Qt(1)))
			throw std::runtime_error("Q1-inconsistency");
		computeQtp1(1); // m1 is fixed
		if (!masked_value_QtQtm1(2).check(Qt(2)))
			throw std::runtime_error("Q2-inconsistency");
		
		auto m2range = MA.word_range(2);
		auto m3range = MA.word_range(3);
		auto m4range = MA.word_range(4);
		uint64_t m4attempts = 0, Q7match = 0, Q7success = 0;
		uint64_t Q7ok = 0, Q8ok = 0, Q9ok = 0;
		uint64_t m4tunnelcnt = 0, m4tunnelok = 0;
		for (uint32_t m2 : m2range)
		{
			m[2] = m2;
			computeQtp1(2);
			if (!masked_value_QtQtm1(3).check(Qt(3)))
				continue;
			for (uint32_t m3 : m3range)
			{
				m[3] = m3;
				computeQtp1(3);
				if (!masked_value_QtQtm1(4).check(Qt(4)))
					continue;
//				for (uint32_t m4 : m4range)
				{
					++m4attempts;
//					m[4] = m4;
					computeQtp1(4); // m4 fixed
					if (!masked_value_QtQtm1(5).check(Qt(5)))
						continue;
					computeQtp1(5); // m5 fixed
					if (!masked_value_QtQtm1(6).check(Qt(6)))
						continue;
					computeQtp1(6); // m6 fixed
					if (!masked_value_QtQtm1(7).check(Qt(7)))
						continue;
					auto itend = Q7m10m12m13.equal_range(Qt(7));
				for (auto it = itend.first; it != itend.second; ++it)
				{
					m[10] = it->second[0];//.first;
					m[12] = it->second[1];//.second;
					m[13] = it->second[2];
					computeQtm3(13); // Q10
					computeQtm3(12); // Q9
					computeQtm3(11); // Q8
					uint32_t Q7org = Qt(7);
//					computeQtm3(10); // Q7
//					if (Qt(7) != Q7org)
//						throw std::runtime_error("Q7-inconsistency");
					if (hammingweight(++Q7match) == 1)
						std::cout << "Q7match=" << Q7match << " " << std::flush;

/*						
//					uint32_t Q4m4tunnel = (~Qvaluemask[offset+4] & ~Qprev[offset+5])
//								& Qt(6) & ~Qt(5);
//				uint32_t Q4org = Qt(4);
//				for (uint32_t Q4cur : masked_value(~Q4m4tunnel, Q4org).range())
//				{
					++m4tunnelcnt;
					Qt(4) = Q4cur;
					computeWt(3);
					if (!MA.checkword(3, m[3]))
						continue;
					computeWt(4);
					if (!MA.checkword(4, m[4]))
						continue;
					computeWt(5);
					if (!MA.checkword(5, m[5]))
						continue;
					computeWt(6);
					if (!MA.checkword(6, m[6]))
						continue;
					++m4tunnelok;
*/

					computeWt(7);
					if (!MA.checkword(7, m[7]))
						continue;
					if (hammingweight(++Q7ok) == 1)
					{
						std::cout << "Q7ok=" << Q7ok << " " << std::flush;
//						std::cout << "(" << log(double(m4tunnelok)/double(m4tunnelcnt))/log(2.0) << ")" << std::flush;
//						std::cout << "(" << log(double(m4tunnelcnt)/double(Q7match))/log(2.0) << ")" << std::flush;
					}
					computeWt(8);
					if (!MA.checkword(8, m[8]))
						continue;
					if (hammingweight(++Q8ok) == 1)
						std::cout << "Q8ok=" << Q8ok << " " << std::flush;
					computeWt(9);
					if (!MA.checkword(9, m[9]))
						continue;
					
//					if (hammingweight(++Q7success) == 1)
					std::cout << "\nQ7success=" << ++Q7success << " " << std::endl;
					for (unsigned t = 0; t < 16; ++t)
					{
						std::cout << "m" << t << "=";
						for (unsigned b = 0; b < 4; ++b)
							std::cout << char((m[t]>>(8*b))&0xFF);
						std::cout << " ";
					}
					std::cout << std::endl;
						
					computeQtp1(21); // m10 => Q22 already checked
					computeQtp1(22); // m15 => Q23 already checked
					
					computeQtp1(23); // m4 => Q24
					if (!masked_value_QtQtm1(24).check(Qt(24)))
						continue;
					if (!checkrotationQtQtp1(23))
						continue;
					std::cout << "Q24 satisfied!" << std::endl;
				}
				}
			}
		}
		std::cout << "m4:" << m4attempts << ", Q7success:" << Q7success << std::endl;
	}
	exit(0);
}

void collisionfinding_alphabet_thread::do_step16_statistics()
{
	static const uint64_t N = 1ULL<<20;
	vector< array<uint32_t, 32> > in, out;
	in.reserve(N);
	out.reserve(N);
	
	std::cout << "Generate Q13-Q16:";
	progress_display genpd(N);
	out.resize(N);
	for (size_t i = 0; i < N; ++i,++genpd)
	{
		Qt(13) = masked_value_Qt(13).sample();
		Qt(14) = masked_value_QtQtm1(14).sample();
		Qt(15) = masked_value_QtQtm1(15).sample();
		Qt(16) = masked_value_QtQtm1(16).sample();

		auto& cur = out[i];
		for (int t = 13; t <= 16; ++t)
			cur[offset+t] = Q[offset+t];
	}

	for (int t = 15; t >= 0; --t)
	{
		std::swap(in, out);
		out.clear(); out.reserve(N);
		std::cout << "\n=======\nStep " << t << ":";
		progress_display pd(N);
		uint64_t attempts = 0;
		while (out.size() < N)
		{
			out.emplace_back(); ++pd;
			auto& cur = out.back();
			while (true)
			{
				++attempts;
				// sample random input value
				auto& incur = in[ xrng64() % in.size() ];
				// copy just the Q values needed for step t
				for (unsigned i = offset+t-3; i <= offset+t+1; ++i)
					Q[i] = incur[i];
				auto Qtm3mv = masked_value_QtQtp1(t-3);
				auto mtrange = MA.word_range(t);
				if (Qtm3mv.count() < mtrange.count())
				{
					Qt(t-3) = Qtm3mv.sample();
					computeWt(t);
					if (!MA.checkword(t, m[t]))
						continue;
				} else {
					m[t] = MA.sampleword(t);
					computeQtm3(t);
					if (!Qtm3mv.check(Qt(t-3)))
						continue;
				}
				// found solution
				cur = incur;
				cur[offset+t-3] = Qt(t-3);
				break;
			}
		}
		std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << " success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
	}
	for (int t = 16; t < 26; ++t)
	{
		std::swap(in, out);
		out.clear(); out.reserve(N);
		std::cout << "\n=======\nStep " << t << ":";
		progress_display pd(N);
		uint64_t attempts = 0;
		while (out.size() < N)
		{
			out.emplace_back(); ++pd;
			auto& cur = out.back();
			while (true)
			{
				++attempts;
				// sample random input value
				auto& incur = in[ xrng64() % in.size() ];
				// copy just the Q values needed for step t
				for (unsigned i = offset+t-3; i <= offset+t+1; ++i)
					Q[i] = incur[i];
				auto Qtp1mv = masked_value_QtQtm1(t+1);
				int wt = md5_wt[t];
				auto mtrange = MA.word_range(wt);
				if (Qtp1mv.count() < mtrange.count())
				{
					Qt(t+1) = Qtp1mv.sample();
					computeWt(t);
					if (!MA.checkword(wt, m[wt]))
						continue;
				} else {
					m[wt] = MA.sampleword(wt);
					computeQtp1(t);
					if (!Qtp1mv.check(Qt(t+1)))
						continue;
				}
				// found solution
				cur = incur;
				cur[offset+t+1] = Qt(t+1);
				break;
			}
		}
		std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << " success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
	}
	differentialpath diffpath2 = diffpath;
	for (int t = 1; t <= 26; ++t)
	{
		for (unsigned b = 0; b < 32; ++b)
		{
			if (Qvaluemask[offset+t]&(1<<b))
				continue;
			uint64_t count1 = 0, countnprev = 0;
			for (size_t i = 0; i < N; ++i)
			{
				if (out[i][offset+t] & (uint32_t(1)<<b))
					++count1;
				if ((out[i][offset+t]^out[i][offset+t-1]) & (uint32_t(1)<<b))
					++countnprev;
			}
			double C1F = double(count1)/double(N);
			double CNPF = double(countnprev)/double(N);
			if (CNPF < 0.3 || CNPF > 0.7)
			{
				std::cout << "Q" << t << "[" << b << "] = ! with prob " << CNPF << std::endl;
				if (CNPF < 0.5)
					diffpath2.setbitcondition(t,b,bc_prev);
				else
					diffpath2.setbitcondition(t,b,bc_prevn);
			}
			if (C1F < 0.3 || C1F > 0.7)
			{
				std::cout << "Q" << t << "[" << b << "] = 1 with prob " << C1F << std::endl;
				if (C1F < 0.5)
					diffpath2.setbitcondition(t,b,bc_zero);
				else
					diffpath2.setbitcondition(t,b,bc_one);
			}
		}
	}
	show_path(diffpath2, m_diff);
	std::map<uint8_t,size_t> byte_set[4];
	for (int t = 0; t < 26; ++t)
	{
		int wt = md5_wt[t];
		for (unsigned b = 0; b < 4; ++b)
			byte_set[b].clear();
		for (size_t i = 0; i < N; ++i)
		{
			for (unsigned j = offset+t-3; j <= offset+t+1; ++j)
				Q[j] = out[i][j];
			computeWt(t);
			uint32_t Wt = m[wt];
			for (unsigned b = 0; b < 4; ++b)
				++ byte_set[b][ uint8_t((Wt>>(b*8))&0xFF) ];
		}
		for (unsigned b = 0; b < 4; ++b)
		{
			std::cout << "Byte " << (4*wt+b) << ": ";
			for (auto& cc : byte_set[b])
				if (cc.second >= N/256)
					std::cout << char(cc.first);
			std::cout << std::endl;
		}
	}
}

void collisionfinding_alphabet_thread::do_step16_old()
{
	std::cout << "0" << std::flush;

	uint64 loopcount0 = 0, loopcount1 = 0, loopcount2 = 0;
	while (true)
	{
		++loopcount0;
		// sample Q13-Q16,m1
		Qt(13) = masked_value_Qt(13).sample();
		Qt(14) = masked_value_QtQtm1(14).sample();
		Qt(15) = masked_value_QtQtm1(15).sample();
		Qt(16) = masked_value_QtQtm1(16).sample();

		m[1] = MA.sampleword(1);
		// check Q17
		computeQtp1(16);
		if (!masked_value_QtQtm1(17).check(Qt(17)))
			continue;

		if (!checkrotationQtQtp1(13) || !checkrotationQtQtp1(14) || !checkrotationQtQtp1(15) || !checkrotationQtQtp1(16)) continue;

		// iterate over valid Q12 until valid m15
		if (++loopcount1 == 1) std::cout << "Q17 " << std::flush;

		auto Q12mv = masked_value_QtQtp1(12);
		for (uint32_t Q12cur : Q12mv.range_rnd())
		{
			Qt(12) = Q12cur;
			if (!checkrotationQtQtp1(12))
				continue;
			computeWt(15);
			if (MA.checkword(15, m[15]))
				break;
		}
		if (!MA.checkword(15, m[15]))
			continue;
		if (++loopcount2 == 1) std::cout << "Q12 " << std::flush;
		
		static const unsigned minQ9m9tunnels = 12;
		// generate table of valid (m14,Q11) solutions with hw(Q11 & ~Q9valuemask) >= minQ9m9tunnel
		m14Q11.clear();
		uint32_t Q9m9tunnel = ~Qtvaluemask(9);
		auto Q11mv = masked_value_QtQtp1(11);
		// ensure at least minQ9m9tunnels number of 1-bits
		for (uint32_t Q11cur : bit_mask_minweight_range(~Q11mv.mask, minQ9m9tunnels))
		{
			Qt(11) = Q11cur ^ Q11mv.value;
			if (!checkrotationQtQtp1(11))
				continue;
			computeWt(14);
			if (!MA.checkword(14, m[14]))
				continue;
			if (hammingweight(Qt(11) & Q9m9tunnel) < minQ9m9tunnels)
				continue;
			m14Q11.emplace_back(m[14], Qt(11));
		};
		if (m14Q11.empty())
			continue;
		// sort (m14,Q11) solutions on decending number of 1-bits of Q11 to maximize Q9m9 tunnel
		std::sort(m14Q11.begin(), m14Q11.end(), [Q9m9tunnel](auto& l, auto& r) { return hw(l.second & Q9m9tunnel) > hw(r.second & Q9m9tunnel); });
		std::cout << " Q11hw:" << hw(m14Q11.front().second) << std::flush;

		m14m13.clear();
		for (auto&& mQ : m14Q11)
		{
//			if (m14m13.size() > (1<<16))
//				break;
			m[14] = mQ.first;
			Qt(11) = mQ.second;
			uint32_t Q9m9tunnel2 = Q9m9tunnel & Qt(11);
			if (hw(Q9m9tunnel2) < minQ9m9tunnels)
				continue;
			auto Q10mv = masked_value_QtQtp1(10);
			uint32_t Q10tunnelmask = ~Q10mv.mask & Q9m9tunnel2;
			// ensure at least minQ9m9tunnels number of 1-bits under Q10tunnelmask
			for (uint32_t Q10curneg : bit_mask_minweight_range(Q10tunnelmask, minQ9m9tunnels))
			{
				// flip bits as for Q9m9 tunnel the Q10 bits actually need to be 0
				uint32_t Q10cur = ~Q10curneg & Q10tunnelmask;
				for (uint32_t Q10cur2 : bit_mask_range(~Q10mv.mask & ~Q9m9tunnel2))
				{
					Qt(10) = Q10cur ^ Q10cur2 ^ Q10mv.value;
					computeWt(13);
					if (!MA.checkword(13, m[13]))
						continue;
					m14m13.emplace_back(m[14],m[13]);
				}
			}
		}
		if (m14m13.empty())
			continue;
		const unsigned Q789conds = hammingweight(masked_value_QtQtp1(9).mask) + hammingweight(masked_value_QtQtp1(8).mask) + hammingweight(masked_value_QtQtp1(7).mask);
		uint64_t expectedQ789sols = (uint64_t(m14m13.size()) * uint64_t(MA.word_size[12])) >> Q789conds;
		

		testcounts[0] += loopcount0;
		testcounts[100] = std::max<uint64>(loopcount0, testcounts[100]);
		testcounts[1] += loopcount1;
		testcounts[101] = std::max<uint64>(loopcount1, testcounts[101]);
		std::cout << "\nQ12-Q17,m1,m15 (l0:" << loopcount0 << ", l1:" << loopcount1 << ", m14m13:" << m14m13.size() << ", E|Q789|:" << expectedQ789sols << ")" << std::endl;
		if (expectedQ789sols >= Q7m12m13m14size)
			std::cout << "D" << std::endl;
		else
			continue;
		do_step17();
		break;
	}
}

void collisionfinding_alphabet_thread::filltables()
{
	isinfinite = true;
	for (int t = -3; t <= 64; ++t)
	{
		dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
		Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
		Qprev[offset+t] = Qprev2[offset+t] = 0;
	}
	Q5m5tunnel=0; Q4m5tunnel=0; Q4m4tunnel=0;
	Q10m10tunnel=0; Q9m10tunnel=0; Q9m9tunnel=0;
	Q8Q12m15tunnel=0; Q14m6Q3m5tunnel=0; Q14Q3m14tunnel=0;

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
	for (int t = 0; t < 64; ++t)
	{
		booleanfunction* F = 0;
		if (t < 16) F = & MD5_F_data;
		else if (t < 32) F = & MD5_G_data;
		else if (t < 48) F = & MD5_H_data;
		else F = & MD5_I_data;
		uint32 dF = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			bf_outcome outcome = F->outcome(diffpath(t,b), diffpath(t-1,b), diffpath(t-2,b));
			if (outcome.size()) {
				if (outcome[0] == bc_plus) 			dF += 1<<b;
				else if (outcome[0] == bc_minus)	dF -= 1<<b;
			} else throw std::runtime_error("ambiguous path!!");
		}
		dT[offset+t] = dQ[offset+t-3] + dF + m_diff[md5_wt[t]];
		dR[offset+t] = dQ[offset+t+1] - dQ[offset+t];
	}

	Q[0] = Qvalue[0]; Q2[0] = Q[0] + dQ[0];
	Q[1] = Qvalue[1]; Q2[1] = Q[1] + dQ[1];
	Q[2] = Qvalue[2]; Q2[2] = Q[2] + dQ[2];
	Q[3] = Qvalue[3]; Q2[3] = Q[3] + dQ[3];
}




































struct collfind_thread {
	collfind_thread(const vector<differentialpath>& paths, parameters_type& params)
		: diffpaths(paths), parameters(params)
	{}
	const vector<differentialpath>& diffpaths;
	const parameters_type& parameters;
	collisionfinding_thread worker;
	void operator()() {
		try {
			for (unsigned i = 0; i < 16; ++i)
				worker.m_diff[i] = parameters.m_diff[i];

#ifdef CPUPERFORMANCE
			for (unsigned t = 0; t < 64; ++t)
				worker.cpu_step_t[t] = 0;
			performance_counter_manager counter_man;
			counter_man.add_performance_counter(worker.cpu_step_t[16], "Step t=16");
			counter_man.add_performance_counter(worker.cpu_step_t[17], "Step t=17");
			counter_man.add_performance_counter(worker.cpu_step_t[18], "Step t=18");
			counter_man.add_performance_counter(worker.cpu_step_t[19], "Step t=19");
			counter_man.add_performance_counter(worker.cpu_step_t[20], "Step t=20");
			counter_man.add_performance_counter(worker.cpu_step_t[21], "Step t=21");
			counter_man.add_performance_counter(worker.cpu_step_t[22], "Step t=22");
			counter_man.add_performance_counter(worker.cpu_step_t[23], "Step t=23");
			counter_man.add_performance_counter(worker.cpu_step_t[24], "Step t=24");
			counter_man.add_performance_counter(worker.cpu_step_t[25], "Step t=25");
#endif
			timer sw(true);
			try {
				worker.findcollision(diffpaths, false);
				while (true) {
					worker.findcollision(diffpaths);
					if (sw.time() > 300) {
						sw.start();
#ifdef CPUPERFORMANCE
						for (unsigned t = 0; t+1 < 64; ++t)
							worker.cpu_step_t[t] -= worker.cpu_step_t[t+1];
						counter_man.show_results();
						for (unsigned t = 62; t > 0; --t)
							worker.cpu_step_t[t] += worker.cpu_step_t[t+1];
#endif
						mut.lock();
						for (unsigned i = 0; i < worker.testcounts.size(); ++i)
							if (worker.testcounts[i])
								cout << i << ":\t" << worker.testcounts[i] << endl;
						mut.unlock();
						worker.findcollision(diffpaths, false);
					}
				}
			} catch (...) {
#ifdef CPUPERFORMANCE
				for (unsigned t = 0; t+1 < 64; ++t)
					worker.cpu_step_t[t] -= worker.cpu_step_t[t+1];
				counter_man.show_results();
#endif
				throw;
			}
		} catch (std::exception & e) { cerr << "Worker thread: caught exception:" << endl << e.what() << endl; } catch (...) {}
	}
};
void collfind_threaded(const vector<differentialpath>& paths, parameters_type& parameters)
{
	boost::thread_group mythreads;
	for (unsigned i = 0; i < parameters.threads; ++i)
		mythreads.create_thread(collfind_thread(paths, parameters));
	mythreads.join_all();
}


struct collfind_alphabet_thread {
	collfind_alphabet_thread(const vector<differentialpath>& paths, parameters_type& params)
		: diffpaths(paths), parameters(params)
	{}
	const vector<differentialpath>& diffpaths;
	const parameters_type& parameters;
	collisionfinding_alphabet_thread worker;
	void operator()() {
		try {
			for (unsigned i = 0; i < 16; ++i)
				worker.m_diff[i] = parameters.m_diff[i];

#ifdef CPUPERFORMANCE
			for (unsigned t = 0; t < 64; ++t)
				worker.cpu_step_t[t] = 0;
			performance_counter_manager counter_man;
			counter_man.add_performance_counter(worker.cpu_step_t[16], "Step t=16");
			counter_man.add_performance_counter(worker.cpu_step_t[17], "Step t=17");
			counter_man.add_performance_counter(worker.cpu_step_t[18], "Step t=18");
			counter_man.add_performance_counter(worker.cpu_step_t[19], "Step t=19");
			counter_man.add_performance_counter(worker.cpu_step_t[20], "Step t=20");
			counter_man.add_performance_counter(worker.cpu_step_t[21], "Step t=21");
			counter_man.add_performance_counter(worker.cpu_step_t[22], "Step t=22");
			counter_man.add_performance_counter(worker.cpu_step_t[23], "Step t=23");
			counter_man.add_performance_counter(worker.cpu_step_t[24], "Step t=24");
			counter_man.add_performance_counter(worker.cpu_step_t[25], "Step t=25");
#endif
			timer sw(true);
			try {
				while (true) {
					worker.findcollision(diffpaths, parameters.alphabet);
					if (sw.time() > 300) {
						sw.start();
#ifdef CPUPERFORMANCE
						for (unsigned t = 0; t+1 < 64; ++t)
							worker.cpu_step_t[t] -= worker.cpu_step_t[t+1];
						counter_man.show_results();
						for (unsigned t = 62; t > 0; --t)
							worker.cpu_step_t[t] += worker.cpu_step_t[t+1];
#endif
						mut.lock();
						for (unsigned i = 0; i < worker.testcounts.size(); ++i)
							if (worker.testcounts[i])
								cout << i << ":\t" << worker.testcounts[i] << endl;
						mut.unlock();
					}
				}
			} catch (...) {
#ifdef CPUPERFORMANCE
				for (unsigned t = 0; t+1 < 64; ++t)
					worker.cpu_step_t[t] -= worker.cpu_step_t[t+1];
				counter_man.show_results();
#endif
				throw;
			}
		} catch (std::exception & e) { cerr << "Worker thread: caught exception:" << endl << e.what() << endl; } catch (...) {}
	}
};
void collfind_alphabet_threaded(const vector<differentialpath>& paths, parameters_type& parameters)
{
	collfind_alphabet_thread mythread(paths, parameters);
	mythread();
	return;
	boost::thread_group mythreads;
	for (unsigned i = 0; i < parameters.threads; ++i)
		mythreads.create_thread(collfind_alphabet_thread(paths, parameters));
	mythreads.join_all();
}

uint32 testtest(uint32 mask, uint32 value)
{
	uint32 y = 0;
	for (auto&& x : bit_mask_range(mask))
	{
		y += (x ^ value);
	}
	return y;
}

int collisionfinding(parameters_type& parameters)
{


	parameters.show_mdiffs();

	differentialpath diffpath;
	vector<differentialpath> vecpath;
	bool failed = true;
	try {
		load_gz(vecpath, binary_archive, parameters.infile1);
		failed = false;
	} catch (...) {}
	if (failed)
	{
		vecpath.clear();
		try {
			load_gz(diffpath, binary_archive, parameters.infile1);
			vecpath.push_back(diffpath);
			failed = false;
		} catch (...) {}
	}
	if (failed || vecpath.size() == 0) {
		cerr << "Error: could not load path(s) in '" << parameters.infile1 << "'!" << endl;
		return 1;
	}
	show_path(vecpath[0], parameters.m_diff);
	cout << "Starting..." << endl;

	if (parameters.alphabet.empty())
		collfind_threaded(vecpath, parameters);
	else
		collfind_alphabet_threaded(vecpath, parameters);
	return 0;
}

