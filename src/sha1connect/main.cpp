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
#include <hashclash/booleanfunction.hpp>

using namespace hashclash;
using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

boost::mutex mut;
std::string workdir;
uint64 pbcount = 0;
int lastlowpathindex = -1;
vector< vector<unsigned> > badcounts(85, vector<unsigned>(32,0));

int main(int argc, char** argv) 
{
	hashclash::timer runtime(true);

	try {
		path_container container;
		vector< vector<unsigned> > msgmask(16);
		unsigned msgoffset = 0;
		bool msglcext = false;
		string tunnelsfilename;
				
		// Define program options
		po::options_description 
			desc("Allowed options"), 
			msg("Define message differences (as +bitnr and -bitnr, bitnr=1-32)"), 
			all("Allowed options");

		desc.add_options()
			("help,h", "Show options.")
			("mod,m"
				, po::value<unsigned>(&container.modn)->default_value(1)
				, "Do only 1/m of all work.")

			("index,i"
				, po::value<unsigned>(&container.modi)->default_value(0)
				, "Do i-th part of all work (i=0,...,m-1).")

			("workdir,w"
				, po::value<string>(&workdir)->default_value("./data")
				, "Set working directory.")

			("inputfilelow"
				, po::value<string>(&container.inputfilelow)
				, "Use specified inputfile for lower paths.")
                        ("loworder"
                                , po::value<unsigned>(&container.loworder)->default_value(0)
                                , "Order of processing low paths (0=random, 1=sorted, 2=w-sorted)")

			("inputfilehigh"
				, po::value<string>(&container.inputfilehigh)
				, "Use specified inputfile for upper paths.")
				
			("inputfileredo"
				, po::value<string>(&container.inputfileredo)
				, "Use specified inputfile to retest\nspecific lower and upper paths.")

			("showinputpaths,s"
				, po::bool_switch(&container.showinputpaths)
				, "Show all input paths.\n")

			("tstep,t"
				, po::value<unsigned>(&container.t)
				, "Do step t (=0,...,15).")

			("Qcondstart"
				, po::value<int>(&container.Qcondstart)->default_value(64)
				, "Count conditions from Q_t up.")
			("maxcond"
				, po::value<unsigned>(&container.bestpathcond)->default_value(262144)
				, "Maximum # conds")
			("keepall"
				, po::bool_switch(&container.keepall)
				, "Store all paths with cond <= maxcond")
			("determinelowestcond"
				, po::bool_switch(&container.determinelowestcond)
				, "Do not generate paths, only determine cond.")
			("splitmode"
				, po::value<unsigned>(&container.splitmode)
				, "Set splitmode: 0=upperpath, 2=lowerpath, +1=random")
                                
                        ("tunnelconditionsfile"
                                , po::value<string>(&tunnelsfilename)->default_value("")
                                , "Force compatibility with tunnels")
			("threads"
				, po::value<int>(&container.threads)->default_value(-1)
				, "Number of worker threads")
                        ("showstats"
                                , po::bool_switch(&container.showstats)
                                , "Show internal statistics for every thread every hour")
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
		all.add(desc).add(msg);
	
		// Parse program options
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, all), vm);
		{
			std::ifstream ifs("sha1diffpathconnect.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, all), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help") || vm.count("tstep")==0) {
			cout << desc << endl;
			return 0;
		}
		if (container.inputfilelow.size() == 0)
		{
			cerr << "Error: inputfilelow must be given!" << endl;
			return 2;
		}
		if (container.inputfilehigh.size() == 0)
		{
			cerr << "Error: inputfilehigh must be given!" << endl;
			return 2;
		}
		if (container.modi >= container.modn) {
			cerr << "Error: i must be strictly less than m!" << endl;
			return 1;
		}
		cout << "Total bitconditions bound: " << container.bestpathcond << endl;
		if (container.modn != 1)
			cout << "Using set " << container.modi << " of " << container.modn << endl;

		if (container.threads <= 0 || container.threads > boost::thread::hardware_concurrency())
			container.threads = boost::thread::hardware_concurrency();
                cout << "Using " << container.threads << " threads." << endl;
                
		uint32 me[16];
		for (unsigned k = 0; k < 16; ++k) {			
			me[k] = 0;
			for (unsigned j = 0; j < msgmask[k].size(); ++j)
				me[k] |= 1<<msgmask[k][j];
		}
		sha1_me_generalised(container.m_mask, me, msgoffset);
		if (msglcext) {
			uint32 met[16];
			uint32 mex[80];
			for (unsigned i = 0; i < 16; ++i)
				met[i] = rotate_left(me[i],5);
			sha1_me_generalised(mex, met, msgoffset+1);
			for (unsigned i = 0; i < 80; ++i)
				container.m_mask[i] ^= mex[i];

			sha1_me_generalised(mex, me, msgoffset+2);
			for (unsigned i = 0; i < 80; ++i)
				container.m_mask[i] ^= mex[i];

			for (unsigned i = 0; i < 16; ++i)
				met[i] = rotate_left(me[i],30);
			for (unsigned j = 3; j < 6; ++j) {
				sha1_me_generalised(mex, met, msgoffset+j);
				for (unsigned i = 0; i < 80; ++i)
					container.m_mask[i] ^= mex[i];
			}			
		}		

                if (tunnelsfilename.size()) {
                        cout << "Loading tunnel conditions: " << tunnelsfilename << "..." << flush;
                        bool loadedtunnels = false;
                        try {
                                vector<sha1differentialpath> tmptunnel;
                                load_bz2(tmptunnel, text_archive, tunnelsfilename);
                                container.tunnelconditions = tmptunnel[0];
                                loadedtunnels = true;
                        } catch (exception &e) { cout << "error: " << endl << e.what() << endl; } catch (...) { cout << "unknown error!" << endl;}
                        if (loadedtunnels) {
                                cout << "done." << endl;
                                show_path(container.tunnelconditions);
                        }
                }


		// Start job with given parameters
		dostep(container);

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
