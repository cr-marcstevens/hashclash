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

boost::mutex mut;
std::string workdir;

int main(int argc, char** argv) 
{
	hashclash::timer runtime(true);

        cout <<
                "MD5 differential path toolbox\n"
                "Copyright (C) 2009 Marc Stevens\n"
                "http://homepages.cwi.nl/~stevens/\n"
                << endl;

	try {
		parameters_type parameters;
		vector< vector<int> > msgdiff(16);
				
		// Define program options
		po::options_description 
			cmds("Allowed commands"),
			desc("Allowed options"), 
			msg("Define message differences (as +bitnr and -bitnr, bitnr=1-32)"), 
			all("Allowed options");

		cmds.add_options()
			("help,h",				"Show options\n")
			("startnearcollision",	"Use inputfile{1,2} to:\n"
									" - determine next near-collision template\n"
									" - construct partial lower diff. path\n"
									" - construct partial upper diff. path\n"
									" - write md5diffpath_forward.cfg\n"
									" - write md5diffpath_backward.cfg\n"
									" - write md5diffpath_connect.cfg\n")
			("upperpaths",			"Write all partial upper diff. paths\n")
			("findcollision",		"Find nearcollision using path\n"
									"   given by inputfile1\n")
			("convert",				"Convert files between binary and text\n")
			("split", po::value<unsigned>(&parameters.split),
									"Split inputfile1 in given # files\n")
			("join,j", po::value<vector<string> >(&parameters.files),
									"Join files and save to outputfile1\n"
									"Each filename has to be proceeded by -j\n")
			("pathfromtext",		"Load path in text form from inputfile1\n"
									"   and save as paths set to outputfile1\n")
			("pathfromcollision",	"Reconstruct paths from colliding inputfiles")
			("startpartialpathfromfile", "Create partial path from binary file")
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
			("pathtyperange"
				, po::value<unsigned>(&parameters.pathtyperange)->default_value(0)
				, "Increases potential # diffs eliminated per n.c.")

			("skipnc"
				, po::value<unsigned>(&parameters.skipnc)->default_value(0)
				, "Skip a number of near-collision templates")

			("threads"
				, po::value<int>(&parameters.threads)->default_value(-1)
				, "Number of worker threads")
			;

		msg.add_options()	
			("diffm0", po::value< vector<int> >(&msgdiff[0]), "delta m0")
			("diffm1", po::value< vector<int> >(&msgdiff[1]), "delta m1")
			("diffm2", po::value< vector<int> >(&msgdiff[2]), "delta m2")
			("diffm3", po::value< vector<int> >(&msgdiff[3]), "delta m3")
			("diffm4", po::value< vector<int> >(&msgdiff[4]), "delta m4")
			("diffm5", po::value< vector<int> >(&msgdiff[5]), "delta m5")
			("diffm6", po::value< vector<int> >(&msgdiff[6]), "delta m6")
			("diffm7", po::value< vector<int> >(&msgdiff[7]), "delta m7")
			("diffm8", po::value< vector<int> >(&msgdiff[8]), "delta m8")
			("diffm9", po::value< vector<int> >(&msgdiff[9]), "delta m9")
			("diffm10", po::value< vector<int> >(&msgdiff[10]), "delta m10")
			("diffm11", po::value< vector<int> >(&msgdiff[11]), "delta m11")
			("diffm12", po::value< vector<int> >(&msgdiff[12]), "delta m12")
			("diffm13", po::value< vector<int> >(&msgdiff[13]), "delta m13")
			("diffm14", po::value< vector<int> >(&msgdiff[14]), "delta m14")
			("diffm15", po::value< vector<int> >(&msgdiff[15]), "delta m15")
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
			std::ifstream ifs("md5diffpathhelper.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, all), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help")
			|| 0 == vm.count("startnearcollision")
					+vm.count("upperpaths")
					+vm.count("enhancepath")
					+vm.count("findcollision")
					+vm.count("convert")
					+vm.count("split")
					+vm.count("join")
					+vm.count("pathfromtext")
					+vm.count("pathfromcollision")
					+vm.count("startpartialpathfromfile")
					) {
			cout << cmds << desc << endl;
			return 0;
		}

		for (unsigned k = 0; k < 16; ++k)
		{
			parameters.m_diff[k] = 0;
			for (unsigned j = 0; j < msgdiff[k].size(); ++j)
				if (msgdiff[k][j] > 0)
					parameters.m_diff[k] += 1<<(msgdiff[k][j]-1);
				else
					parameters.m_diff[k] -= 1<<(-msgdiff[k][j]-1);
		}

		if (parameters.threads <= 0 || parameters.threads > boost::thread::hardware_concurrency())
			parameters.threads = boost::thread::hardware_concurrency();

		if (vm.count("startnearcollision"))
			return startnearcollision(parameters);
		if (vm.count("upperpaths"))
			return upperpaths(parameters);
		if (vm.count("findcollision"))
			return collisionfinding(parameters);
		if (vm.count("convert"))
			return convert(parameters);
		if (vm.count("split"))
			return split(parameters);
		if (vm.count("join"))
			return join(parameters);
		if (vm.count("pathfromtext"))
			return pathfromtext(parameters);
		if (vm.count("pathfromcollision"))
			return pathfromcollision(parameters);
		if (vm.count("startpartialpathfromfile"))
			return partialpathfromfile(parameters);

		// Start job with given parameters
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
			if (!!i) {
				++len;
				block[k] += uint32((unsigned char)(uc)) << (c*8);
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
			o << (unsigned char)((block[k]>>(c*8))&0xFF);
}
