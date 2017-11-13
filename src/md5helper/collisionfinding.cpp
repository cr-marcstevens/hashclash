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

#define MD5DETAIL_INLINE_IMPL
#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/timer.hpp>

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


int collisionfinding(parameters_type& parameters)
{
	parameters.show_mdiffs();
//	for (unsigned i = 0; i < 16; ++i)
//		m_diff[i] = parameters.m_diff[i];

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

	collfind_threaded(vecpath, parameters);
//	collisionfinding_thread worker;
//	for (unsigned i = 0; i < 16; ++i)
//		worker.m_diff[i] = parameters.m_diff[i];
//	worker.diffpath = diffpath;
//	worker.filltables();
//	while (true) {
//		worker.do_step16();
//	}
	return 0;
}

