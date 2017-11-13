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

#include <hashclash/sdr.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/booleanfunction.hpp>

#include "main.hpp"

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
		"Extend MD5 differential paths backward\n"
		"Copyright (C) 2009 Marc Stevens\n"
		"http://homepages.cwi.nl/~stevens/\n"
		<< endl;

	try {
		path_container_autobalance container;
		vector< vector<int> > msgdiff(16);
				
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

			("newinputpath,n"
				, po::bool_switch(&container.newinputpath)
				, "Start with empty path.")

			("inputfile,f"
				, po::value<string>(&container.inputfile)->default_value("")
				, "Use specified inputfile.")

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
				, po::value<int>(&container.tend)->default_value(64)
				, "Set ending Q_t for maxconditions.")

			("maxQ26upcond"
				, po::value<unsigned>(&container.maxQ26upcond)->default_value(0)
				, "Limit conditions on Q_26 and up")

			("includenaf"
				, po::bool_switch(&container.includenaf)
				, "Include naf(Q_t-3) in # conditions.")

			("halfnafweight"
				, po::bool_switch(&container.halfnafweight)
				, "naf(Q_t-3) only weighs half in # cond.")

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
				, po::value<double>(&container.fillfraction)->default_value(1)
				, "Fill up to fraction of target amount\n\twith heavier paths than maxcond.")

			("estimate,e"
				, po::value<double>(&container.estimatefactor)->default_value(2)
				, "Do an estimate for maxcond.")

			("nafestimate"
				, po::value<unsigned>(&container.nafestweight)->default_value(8)
				, "Estimate of naf(Q_t-3) in # conditions.")

			("noverify"
				, po::bool_switch(&container.noverify)
				, "Do not verify rotations etc.")

			("threads"
				, po::value<int>(&container.threads)->default_value(-1)
				, "Number of worker threads")

			("uct"
				, po::value<int>(&container.uct)->default_value(-4)
				, "Disallow condition: Q_t[b] = c")
			("ucb"
				, po::value<int>(&container.ucb)->default_value(-1)
				, "Disallow condition: Q_t[b] = c")
			("ucc"
				, po::value<char>(&container.ucc)->default_value('.')
				, "Disallow condition: Q_t[b] = c")
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
		all.add(desc).add(msg);
	
		// Parse program options
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, all), vm);
		{
			std::ifstream ifs("md5diffpathbackward.cfg");
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

		for (unsigned k = 0; k < 16; ++k)
		{
			container.m_diff[k] = 0;
			for (unsigned j = 0; j < msgdiff[k].size(); ++j)
				if (msgdiff[k][j] > 0)
					container.m_diff[k] += 1<<(msgdiff[k][j]-1);
				else
					container.m_diff[k] -= 1<<(-msgdiff[k][j]-1);
			if (container.m_diff[k] != 0)
				cout << "delta_m[" << k << "] = " << naf(container.m_diff[k]) << endl;
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

		if (container.threads <= 0 || container.threads > boost::thread::hardware_concurrency())
			container.threads = boost::thread::hardware_concurrency();

		// Start job with given parameters
                container.set_parameters();
		for (unsigned tt = container.t; tt > container.t-container.trange; --tt) {
			path_container_autobalance containertmp = container;
			containertmp.t = tt;
			dostep(containertmp, true);
		}
		container.t -= container.trange;
		dostep(container);
                                                                                                                                                                        
	} catch (exception& e) {
		cout << "Runtime: " << runtime.time() << endl;
		cerr << "Caught exception!!:" << endl << e.what() << endl;
		return 1;
	} catch (...) {
		cout << "Runtime: " << runtime.time() << endl;
		cerr << "Unknown exception caught!!" << endl;
		return 1;
	}
	cout << "Runtime: " << runtime.time() << endl;
	return 0;
}
