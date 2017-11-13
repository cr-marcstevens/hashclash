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

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "main.hpp"

#include <hashclash/sdr.hpp>
#include <hashclash/timer.hpp>

using namespace hashclash;
using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

std::string workdir;

int main(int argc, char** argv) 
{
	hashclash::timer runtime(true);

	try {
		parameters_type parameters;
		vector< vector<unsigned> > msgmask(16);
		unsigned msgoffset = 0;
		bool msglcext = false;

				
		// Define program options
		po::options_description 
			cmds("Allowed commands"),
			desc("Allowed options"), 
			msg("Define message differences (as +bitnr and -bitnr, bitnr=1-32)"), 
			all("Allowed options");

		cmds.add_options()
			("help,h",				"Show options\n")
			("convert",				"Convert files between binary and text\n")

			("split", po::value<unsigned>(&parameters.split),
									"Split inputfile1 in given # files\n")
			("join,j", po::value<vector<string> >(&parameters.files),
									"Join files and save to outputfile1\n"
									"Each filename has to be proceeded by -j\n")
			("pathfromtext",		"Load path in text form from inputfile1\n"
									"   and save as paths set to outputfile1\n")
			("showmespaceconditions", "Load me space from inputfile1\n   and show bitconditions\n")
			("filterfeasible", "Filter non-feasible differential paths")
			("select,s", "Select i-th path from set.")
			;
		desc.add_options()
			("workdir,w"
				, po::value<string>(&workdir)->default_value("./data")
				, "Set working directory.")

			("inputfile1"
				, po::value<string>(&parameters.infile1)->default_value("")
				, "Set inputfile 1.")
			("inputfile2"
				, po::value<string>(&parameters.infile2)->default_value("")
				, "Set inputfile 1.")
			("outputfile1"
				, po::value<string>(&parameters.outfile1)->default_value("")
				, "Set outputfile 1.")
			("outputfile2"
				, po::value<string>(&parameters.outfile2)->default_value("")
				, "Set outputfile 2.")
			("cpuaffinity"
				, po::value<int>(&parameters.cpuaffinity)->default_value(-1)
				, "Limit findcoll. to use only given cpu")
			("unique"
				, po::bool_switch(&parameters.unique)
				, "Only keep unique paths when joining or splitting.")
			("filtert"
				, po::value<int>(&parameters.filtert)->default_value(-4)
				, "Filter paths by having solutions up to Q_filtert.")
			("seli"
				, po::value<int>(&parameters.seli)->default_value(0)
				, "Select i-th path.")
			("invert"
				, po::bool_switch(&parameters.invert)
				, "Inverted selected path")
			("mod,m"
				, po::value<unsigned>(&parameters.mod)->default_value(1)
				, "Number of processes to split over.")
			("index,i"
				, po::value<unsigned>(&parameters.index)->default_value(0)
				, "Index of this process (icm --mod).")
			;

		msg.add_options()	
			("diffm0", po::value< vector<unsigned> >(&msgmask[0]), "mask m0")
			("diffm1", po::value< vector<unsigned> >(&msgmask[1]), "mask m1")
			("diffm2", po::value< vector<unsigned> >(&msgmask[2]), "mask m2")
			("diffm3", po::value< vector<unsigned> >(&msgmask[3]), "mask m3")
			("diffm4", po::value< vector<unsigned> >(&msgmask[4]), "mask m4")
			("diffm5", po::value< vector<unsigned> >(&msgmask[5]), "mask m5")
			("diffm6", po::value< vector<unsigned> >(&msgmask[6]), "mask m6")
			("diffm7", po::value< vector<unsigned> >(&msgmask[7]), "mask m7")
			("diffm8", po::value< vector<unsigned> >(&msgmask[8]), "mask m8")
			("diffm9", po::value< vector<unsigned> >(&msgmask[9]), "mask m9")
			("diffm10", po::value< vector<unsigned> >(&msgmask[10]), "mask m10")
			("diffm11", po::value< vector<unsigned> >(&msgmask[11]), "mask m11")
			("diffm12", po::value< vector<unsigned> >(&msgmask[12]), "mask m12")
			("diffm13", po::value< vector<unsigned> >(&msgmask[13]), "mask m13")
			("diffm14", po::value< vector<unsigned> >(&msgmask[14]), "mask m14")
			("diffm15", po::value< vector<unsigned> >(&msgmask[15]), "mask m15")
			("diffoffset", po::value< unsigned >(&msgoffset), "mask offset")
			("difflcext", po::bool_switch(&msglcext), "extend to local coll.")
			;
		all.add(cmds).add(desc).add(msg);
	
		// Parse program options
		po::positional_options_description p;
		p.add("inputfile1", 1);
		p.add("inputfile2", 1);
		p.add("outputfile1", 1);
		p.add("outputfile2", 1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv)
			.options(all).positional(p).run(), vm);
		{
			std::ifstream ifs("sha1diffpathhelper.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, all), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help")
			|| 0 == 0
					+vm.count("convert")
					+vm.count("split")
					+vm.count("join")
					+vm.count("pathfromtext")
					+vm.count("showmespaceconditions")
					+vm.count("filterfeasible")
					+vm.count("select")
					) {
			cout << cmds << desc << endl;
			return 0;
		}

		uint32 me[16];
		for (unsigned k = 0; k < 16; ++k) {			
			me[k] = 0;
			for (unsigned j = 0; j < msgmask[k].size(); ++j)
				me[k] |= 1<<msgmask[k][j];
		}
		sha1_me_generalised(parameters.m_mask, me, msgoffset);
		if (msglcext) {
			uint32 met[16];
			uint32 mex[80];
			for (unsigned i = 0; i < 16; ++i)
				met[i] = rotate_left(me[i],5);
			sha1_me_generalised(mex, met, msgoffset+1);
			for (unsigned i = 0; i < 80; ++i)
				parameters.m_mask[i] ^= mex[i];

			sha1_me_generalised(mex, me, msgoffset+2);
			for (unsigned i = 0; i < 80; ++i)
				parameters.m_mask[i] ^= mex[i];

			for (unsigned i = 0; i < 16; ++i)
				met[i] = rotate_left(me[i],30);
			for (unsigned j = 3; j < 6; ++j) {
				sha1_me_generalised(mex, met, msgoffset+j);
				for (unsigned i = 0; i < 80; ++i)
					parameters.m_mask[i] ^= mex[i];
			}			
		}
		if (vm.count("convert"))
			return convert(parameters);
		if (vm.count("split"))
			return split(parameters);
		if (vm.count("join"))
			return join(parameters);
		if (vm.count("pathfromtext"))
			return pathfromtext(parameters);
		if (vm.count("showmespaceconditions"))
			return showmespaceconditions(parameters);
		if (vm.count("filterfeasible"))
			return filterfeasible(parameters);
		if (vm.count("select"))
			return selectpath(parameters);
		cout << "huh?!?" << endl;
	} catch (exception& e) {
		cout << "Runtime: " << runtime.time() << endl;
		cerr << "Caught exception!!:" << endl << e.what() << endl;
		throw;
	} catch (...) {
		cout << "Runtime: " << runtime.time() << endl;
		cerr << "Unknown exception caught!!" << endl;
		throw;
	}
	cout << "Runtime: " << runtime.time() << endl;
	return 0;
}



unsigned load_block(istream& i, uint32 block[])
{	
	for (unsigned k = 0; k < 16; ++k)
		block[k] = 0;

	unsigned len = 0;
	char uc;
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			i.get(uc);
			if (i) {
				++len;
				block[k] += uint32((unsigned char)(uc)) << ((3-c)*8);
			} else {
				i.putback(uc);
				i.setstate(ios::failbit);
				return len;
			}
		}
	return len;
}

void save_block(ostream& o, uint32 block[])
{
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
			o << (unsigned char)((block[k]>>((3-c)*8))&0xFF);
}
