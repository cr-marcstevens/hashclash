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

#include <hashclash/sdr.hpp>
#include <hashclash/saveload_gz.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/rng.hpp>

using namespace hashclash;
using namespace std;

extern boost::mutex mut;
extern std::string workdir;

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, uint32 block[]);

struct parameters_type {
	uint32 m_diff[16];
	string infile1, infile2, outfile1, outfile2;
	unsigned split;
	unsigned skipnc;
	unsigned pathtyperange;
	vector<string> files;
	int threads;

	void show_mdiffs()
	{
		for (unsigned k = 0; k < 16; ++k)
			if (m_diff[k] != 0)
				cout << "delta_m[" << k << "] = " << naf(m_diff[k]) << endl;
	}
};

int upperpaths(parameters_type& parameters);
int startnearcollision(parameters_type& parameters);
int convert(parameters_type& parameters);
int split(parameters_type& parameters);
int join(parameters_type& parameters);
int pathfromtext(parameters_type& parameters);
int pathfromcollision(parameters_type& parameters);
int collisionfinding(parameters_type& parameters);
int partialpathfromfile(parameters_type& parameters);

struct collisionfinding_thread {
	collisionfinding_thread(): testcounts(1<<20,0) {}

	int verifyconds();
	void do_step26();
	void do_step25();
	void do_step24();
	void do_step23();
	//void do_step22();
	void do_step21();
	void do_step20();
	//void do_step19();
	void do_step18();
	void do_step17();
	void do_step16();
	void filltables();
	void findcollision(const vector<differentialpath>& diffpaths, bool keeppath = true) 
	{
		if (!keeppath) {
			mut.lock();
			diffpath = diffpaths[xrng64() % diffpaths.size()];
			mut.unlock();
			filltables();
		}
		do_step16();
	}

	uint64 cpu_step_t[64];
	vector<uint64> testcounts/*(1<<20,0)*/;

	differentialpath diffpath;
	uint32 m_diff[16];

	vector< triple<uint32,uint32,uint32> > Q1Q2m0withm1ok;
	vector< triple<uint32,uint32,uint32> >::const_iterator Q1Q2m0it, Q1Q2m0itend;

	static const int offset = 3;
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

	uint32 Q5m5tunnel, Q4m5tunnel, Q4m4tunnel;
	uint32 Q10m10tunnel, Q9m10tunnel, Q9m9tunnel;
	uint32 Q8Q12m15tunnel, Q14m6Q3m5tunnel, Q14Q3m14tunnel;
	bool isinfinite;
};

#endif // MAIN_HPP
