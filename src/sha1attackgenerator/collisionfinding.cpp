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
#include <iomanip>
#include <algorithm>
#include <map>
#include <fstream>
#include <stdexcept>

#define SHA1DETAIL_INLINE_IMPL
#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1detail.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/bestof.hpp>
#include <hashclash/progress_display.hpp>
#include <hashclash/sha1messagespace.hpp>

#include <boost/lexical_cast.hpp>

#include "main.hpp"

sha1differentialpath maindiffpath;
sha1messagespace mainmespace;

uint32 m_diff[80];

bool throw_at_t16 = false;
uint32 firstmsg[16];

vector< uint32 > Qtp1valsvec[30];

vector< vector< vector<uint32> > > bitrelations;
vector< vector<uint32> > pathbitrelations;
vector< vector< vector<uint32> > > pathbitrelationsmatrix;

vector< vector<uint32> > firstmsg_vec;

void filter_bitconditions(vector< vector<uint32> >& bitrels, unsigned tbegin, unsigned tend)
{
	sweep_matrix(bitrels,80*32);
	vector<uint32> triv(81,0), impos(81,0);
	impos[80] = ~uint32(0);
	unsigned i = 0;
	while (i < bitrels.size()) {
		if (bitrels[i] == triv) {
			swap(bitrels[i], bitrels.back());
			bitrels.pop_back();
			continue;
		}
		if (bitrels[i] == impos)
			throw std::runtime_error("filter_bitconditions: contradicting bitrelations");
		
		int lastt = -1;
		for (int j = 0; j < bitrels[i].size() && j < 80; ++j)
			if (bitrels[i][j])
				lastt = j;
		if (lastt == -1)
			throw std::runtime_error("filter_bitconditions: huh?!?");
		if (lastt < tbegin || lastt >= tend) {
			swap(bitrels[i], bitrels.back());
			bitrels.pop_back();
		} else
			++i;
	}
}

int collisionfinding(parameters_type& parameters)
{
	bool usetunnelbitconditions = parameters.usetunnelbitconditions;

	sha1messagespace tmpspace;
	vector< vector<uint32> > bitrels, tmpbitrel, tmpbitrel2;
	for (unsigned i = 0; i < parameters.rnd234_m_bitrelationfiles.size(); ++i) {
		try {
			cout << "Loading '" << parameters.rnd234_m_bitrelationfiles[i] << "'..." << flush;
			load_bz2(tmpspace, text_archive, parameters.rnd234_m_bitrelationfiles[i]);
			cout << "done" << flush;
			tmpspace.tobitrelations_80(tmpbitrel);
			bitrels.insert(bitrels.end(), tmpbitrel.begin(), tmpbitrel.end());
			cout << " (" << tmpbitrel.size() << " bitrels, new total: " << bitrels.size() << ")" << endl;
		} catch (std::exception& e) {
			cerr << "Exception:" << endl << e.what() << endl;
			return 1;
		}
	}
	vector< sha1differentialpath > paths;

	cout << "Loading round 1 paths from '" << parameters.rnd1_pathsfile << "'..." << flush;
	try {
		load_bz2(paths, text_archive, parameters.rnd1_pathsfile);
		cout << "done (" << paths.size() << " paths)." << endl;
	} catch (std::exception& e) {
		cerr << "Exception:" << endl << e.what() << endl;
		return 1;
	}
	if (paths.size() == 0) exit(0);

	vector< vector<uint32> > bitrelsbu = bitrels;
for (unsigned pi = 0; pi < paths.size(); ++pi) { try {
	bitrels = bitrelsbu;
	sha1differentialpath upperpath;
	cout << "Loading round 2,3,4 path from '" << parameters.rnd234_pathfile << "'..." << flush;
	try {
		load_bz2(upperpath, text_archive, parameters.rnd234_pathfile);
		cout << "done." << endl;
	} catch (std::exception& e) {
		cerr << "Exception:" << endl << e.what() << endl;
		return 1;
	}
	for (int t = paths[pi].tbegin(); t < paths[pi].tend(); ++t)
		upperpath[t] = paths[pi][t];
	for (unsigned t = 0; t < paths[pi].tend()-1; ++t)
		upperpath.getme(t) = paths[pi].getme(t);
	cleanup_path(upperpath);
	maindiffpath = upperpath;

	bool bad = false;
	for (int t = maindiffpath.tbegin()+4; t < maindiffpath.tend()-1; ++t) {
		if (maindiffpath.getme(t).mask != parameters.m_mask[t]) {
			uint32 dm = maindiffpath.getme(t).adddiff();
			uint32 mcur = 0;
			uint32 mmask = parameters.m_mask[t];
			uint32 madd = (~mmask)+1;
			sdr msdr;
			msdr.mask = mmask;
			do {
				mcur += madd; mcur &= mmask;
				msdr.sign = mcur;
				if (msdr.adddiff() == dm) {
					cout << "corrected: \t" << t << ":\t" << msdr << " == " << sdr(parameters.m_mask[t]) << endl;
					break;
				}
			} while (mcur != 0);
			if (msdr.adddiff() == dm) {
				maindiffpath.getme(t) = msdr;
			} else {
				bad = true;
				cout << "failed: \t" << t << ":\t" << maindiffpath.getme(t) << " != " << sdr(parameters.m_mask[t]) << endl;
			}
		}
	}
	if (bad) exit(0);
	show_path(maindiffpath);

	
	// remove bitrelations possibly limiting me[0],...,me[19]
	const unsigned tend_fix_me = parameters.tend_rnd1_me;
	cout << "Fixed me diffs for t=[0," << tend_fix_me << "): (" << bitrels.size() << "=>" << flush;
	filter_bitconditions(bitrels, tend_fix_me, 80);
	cout << bitrels.size() << ")" << endl;
	for (unsigned i = 0; i < bitrels.size(); ++i) {
		cout << " - ";
		bool firstone = true;
		for (unsigned t = 0; t < 80; ++t)
			for (unsigned b = 0; b < 32; ++b)
				if (bitrels[i][t] & (1<<b)) {
					if (firstone)
						firstone = false;
					else
						cout << " + ";
					cout << "M[" << t << "," << b << "]";
				}
		cout << " = " << (bitrels[i][80]&1) << endl;
	}
	
	tmpspace.clear();
	for (unsigned t = 0; t < 80; ++t) {
		if (t < tend_fix_me) {
			for (unsigned b = 0; b < 31; ++b) {
				int bit = maindiffpath.getme(t).get(b);
				if (bit == 0)
					tmpspace.buildbasis_addfreebit(t,b);
				else
					tmpspace.buildbasis_setbit(t,b,bit==-1);
			}
			tmpspace.buildbasis_addfreebit(t,31);
		} else {
			for (unsigned b = 0; b < 32; ++b)
				tmpspace.buildbasis_addfreebit(t,b);
		}
	}
	tmpspace.tobitrelations_80(tmpbitrel);
	bitrels.insert(bitrels.end(), tmpbitrel.begin(), tmpbitrel.end());
	cout << "Extra me [0," << tend_fix_me << ") bitrelations: " << tmpbitrel.size() << endl;
	{
		cout << "Extra tunnel me bitrelations: " << endl;
		ifstream ifs("tunnel_bitconditions.txt");
		read_message_bitconditions(ifs, tmpbitrel);
		for (unsigned i = 0; i < tmpbitrel.size(); ++i) {
			bool firstone = true;
			for (unsigned t = 0; t < 80; ++t)
				if (tmpbitrel[i][t]) {
					if (firstone) firstone = false;
					else cout << "+ ";
					cout << "m" << t << sdr(tmpbitrel[i][t]) << " ";
				}
			cout << "= " << (tmpbitrel[i][80]&1) << endl;
		}
		bitrels.insert(bitrels.end(), tmpbitrel.begin(), tmpbitrel.end());
	}
	tmpspace.frombitrelations_80(bitrels);
	mainmespace = tmpspace;
	tmpspace.tobitrelations_16(bitrels);
	cout << "Total bitrelations: " << bitrels.size() << endl;
#if 0
	cout << "Independent message bits: " << endl;
	uint32 meindep[16];
	for (unsigned i = 0; i < 16; ++i)
		meindep[16] = ~uint32(0);
	for (unsigned i = 0; i < bitrels.size(); ++i)
		for (unsigned t = 0; t < 16; ++t)
			meindep[t] &= ~bitrels[i][t];
	for (unsigned i = 0; i < 16; ++i) {
		cout << "m[" << i << "] indep:\t";
		for (int b = 31; b >= 0; --b)
			cout << ((meindep[i]&(1<<b))?"1":"0");
		cout << endl;
	}
#endif
	pathbitrelations = bitrels;
	pathbitrelationsmatrix.clear();
	pathbitrelationsmatrix.resize(16);
	for (unsigned i = 0; i < 16; ++i)
		pathbitrelationsmatrix[i].resize(32);
	for (unsigned i = 0; i < pathbitrelations.size(); ++i) {
		int lastcol = -1;
		for (int col = 16*32-1; col >= 0; --col)
			if (pathbitrelations[i][col>>5]&(1<<(col&31))) {
				lastcol = col;
				unsigned t = lastcol>>5;
				unsigned b = lastcol&31;
				pathbitrelationsmatrix[t][b] = pathbitrelations[i];
				break;
			}
		if (lastcol == -1) throw;
	}
#if 0
	//needs diffpath & pathbitrelationsmatrix
	analyze_indepsection_prob();
#endif
	if (parameters.tunnelfile.size()) {
		vector<sha1differentialpath> tunnels;
		cout << "Loading tunnels from '" << parameters.tunnelfile << "'..." << flush;
		try {
			load_bz2(tunnels, text_archive, parameters.tunnelfile);
		} 
		catch (std::exception &e) { tunnels.clear(); cout << "failed:" << endl << e.what() << endl; }
		catch (...) { tunnels.clear(); cout << "failed." << endl; }
		if (tunnels.size())
			analyze_tunnels_diffpath(maindiffpath, bitrels, tunnels);
		exit(0);
	}
	if (!usetunnelbitconditions) {
		// if not using tunnel bit conditions then we'll assume we want to analyze tunnels
		analyze_tunnels_diffpath(maindiffpath, bitrels);
		exit(0);
	}
	// if we're using tunnel bit conditions we'll assume we want to generate the collision finding program

	generate_program();
	break;
} catch (std::exception& e) { cerr << "c: " << e.what() << endl; } catch (...) {} 
} // for (unsigned pi = 0; pi < paths.size(); ++pi)
}


