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

using namespace hashclash;
using namespace std;

extern std::string workdir;

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, uint32 block[]);

struct parameters_type {
	unsigned mod,index;
	uint32 m_mask[80];
	string infile1, infile2, outfile1, outfile2;
	int cpuaffinity;
	unsigned split;
	vector<string> files;
	bool unique;
	int filtert;
	int seli;
	bool invert;

	void show_mdiffs()
	{
	}
};

int convert(parameters_type& parameters);
int split(parameters_type& parameters);
int join(parameters_type& parameters);
int pathfromtext(parameters_type& parameters);
int showmespaceconditions(parameters_type& parameters);

int filterfeasible(parameters_type& parameters);

int selectpath(parameters_type& parameters);

#endif // MAIN_HPP
