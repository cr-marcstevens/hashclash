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

#include <hashclash/rng.hpp>
#include <hashclash/sdr.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/booleanfunction.hpp>

#include "main.hpp"

using namespace hashclash;
using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

std::string workdir;

int main(int argc, char** argv) 
{
	hashclash::timer runtime(true);

	try {
		path_container_autobalance container;
		vector< vector<unsigned> > msgmask(16);
		unsigned msgoffset = 0;
		bool msglcext = false;
		string tunnelsfilename;
				
		// Define program options
		po::options_description 
			desc("Allowed options"), 
			msg("Define message differences (as bitnr, bitnr=0-31)"), 
			ihv("Define start IV (necessary for t=0 and t=1)"),
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

			("inputfile,f"
				, po::value<string>(&container.inputfile)->default_value("")
				, "Use specified inputfile.")

			("newinputpath,n"
				, po::bool_switch(&container.newinputpath)
				, "Start with empty path.")

			("showinputpaths,s"
				, po::bool_switch(&container.showinputpaths)
				, "Show all input paths.\n")

			("tstep,t"
				, po::value<unsigned>(&container.t)
				, "Do step t (=0,...,15).")
				
			("trange"
				, po::value<unsigned>(&container.trange)->default_value(0)
				, "Number of additional steps to perform.")

			("maxconditions,c"
				, po::value<unsigned>(&container.maxcond)->default_value(2176)
				, "Limit total amount of bitconditions.")

			("condtend,b"
				, po::value<int>(&container.tend)->default_value(80)
				, "Set starting Q_t for maxconditions.")

			("maxQ26upcond"
				, po::value<unsigned>(&container.maxQ26upcond)->default_value(0)
				, "Limit conditions on Q_26 and up")

			("includenaf"
				, po::bool_switch(&container.includenaf)
				, "Include naf(Q_t+1) in # conditions.")

			("halfnafweight"
				, po::bool_switch(&container.halfnafweight)
				, "naf(Q_t+1) only weighs half in # cond.")

			("maxweight,p"
				, po::value<unsigned>(&container.maxweight)->default_value(32)
				, "Limit carry propagation.")
			("minweight"
				, po::value<unsigned>(&container.minweight)->default_value(0)
				, "Limit carry propagation.")

			("maxsdrs,q"
				, po::value<unsigned>(&container.maxsdrs)->default_value(2000)
				, "Limit total number of sdrs.")

			("autobalance,a"
				, po::value<unsigned>(&container.ubound)->default_value(0)
				, "Do autobalancing using target amount.")

			("fillfraction"
				, po::value<double>(&container.fillfraction)->default_value(0)
				, "Fill up to fraction of target amount\n\twith heavier paths than maxcond.")

			("estimate,e"
				, po::value<double>(&container.estimatefactor)->default_value(0)
				, "Do an estimate for maxcond.")

			("nafestimate"
				, po::value<unsigned>(&container.nafestweight)->default_value(8)
				, "Estimate of naf(Q_t+1) in # conditions.")

			("onemessagediff"
				, po::bool_switch(&container.onemessagediff)
				, "Use only one messagediff dmt")
			("expandprevmessagediff"
				, po::bool_switch(&container.expandprevmessagediff)
				, "Expand previous messagediff")
			("tunnelconditionsfile"
				, po::value<string>(&tunnelsfilename)->default_value("")
				, "Force compatibility with tunnels")
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
			std::ifstream ifs("sha1diffpathbackward.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, all), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help") || vm.count("tstep")==0) {
			cout << desc << endl;
			return 0;
		}
		if (container.modi >= container.modn) {
			cerr << "Error: i must be strictly less than m!" << endl;
			return 1;
		}
		if (container.modn != 1) {
			cout << "Using set " << container.modi << " of " << container.modn << endl;
			container.trange = 0;
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

		// do not estimate if not autobalancing
		if (container.ubound == 0)
			container.estimatefactor = 0;

		// only estimate for larger ubound, not lower
		if (container.estimatefactor < 1)
			container.estimatefactor = 0;

		if (container.fillfraction < 0)
			container.fillfraction = 0;
		if (container.fillfraction > 1)
			container.fillfraction = 1;

		// Start job with given parameters
		container.set_parameters();
		for (unsigned tt = container.t; tt > container.t-container.trange; --tt) {
			cout << "delta_m[" << tt << "] = " << sdr(0,container.m_mask[tt]) << endl;
			path_container_autobalance containertmp = container;
			containertmp.t = tt;
			dostep(containertmp, true);
		}
		container.t -= container.trange;
		cout << "delta_m[" << container.t << "] = " << sdr(0,container.m_mask[container.t]) << endl;
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
