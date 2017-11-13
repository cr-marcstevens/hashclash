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

#include <hashclash/sdr.hpp>
#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/sha1messagespace.hpp>

using namespace hashclash;
using namespace std;


struct parameters_type {
	parameters_type()
	{
		usetunnelbitconditions = false;
		tunneldeepanalysis = false;
		changebitstat = false;
		mecarry = 0;
		threebittunnels = true;
		simpletunnels = false;
		mintunnelpov = -1;
	}

	vector<string> rnd234_m_bitrelationfiles;
	string rnd234_pathfile; // single path
	string rnd1_pathsfile; // possibly multiple paths allowing choice without m bitrelation contradictions
	unsigned tend_rnd1_me;
			
	unsigned mod,index;
	uint32 m_mask[80];
	bool usetunnelbitconditions;
	std::string tunnelfile;
	unsigned mecarry;
	bool tunneldeepanalysis;
	bool changebitstat;
	bool tunnelfilterbitcondition;
	bool threebittunnels;
	bool simpletunnels;
	int mintunnelpov;
	void show_mdiffs()
	{
	}
};
extern std::string workdir;
extern sha1differentialpath maindiffpath;
extern sha1messagespace mainmespace;
extern parameters_type mainparameters;
const int offset = 4;

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, uint32 block[]);

int collisionfinding(parameters_type& parameters);


// filters.cpp
void mespace_to_pathbitrelationsmatrix(const sha1messagespace& mespace, vector< vector< vector<uint32> > >& pathbitrelationsmatrix);
void random_me(uint32 me[80], const vector< vector< vector<uint32> > >& pathbitrelationsmatrix);

void filter_tunnels_bitconditions(vector< sha1differentialpath >& tunnels, const sha1differentialpath& path);
bool filter_tunnel_bitrelations(const sha1differentialpath& tunnel, const vector< vector<uint32> >& bitrels);
void filter_tunnels_bitrelations(vector< sha1differentialpath >& tunnels, const vector< vector<uint32> >& bitrels);

// program_generator.cpp
string mask2bit_to_string(const string& variable, const uint32 mask, const unsigned bit);
void generate_program();

// tunnel_analysis.cpp
unsigned analyze_bc_tunnel(const sha1differentialpath& tunnel);
unsigned analyze_tunnel(const sha1differentialpath& tunnel, bool verbose = false, bool deepanalyze = false);
void analyze_tunnels_diffpath(const sha1differentialpath& mypath, const vector< vector<uint32> >& bitrels, vector<sha1differentialpath> tunnels = vector<sha1differentialpath>());

// checkoppaths.cpp
void checkokpaths();

#endif // MAIN_HPP
