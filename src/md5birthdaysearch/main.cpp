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

#include <boost/program_options.hpp>

#include "main.hpp"

#include <hashclash/md5detail.hpp>
#include <hashclash/timer.hpp>

using namespace hashclash;
using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

std::string workdir;

int main(int argc, char** argv) 
{
	int result = 0;
	timer runtime(true);

	cout <<
		"Birthday search for MD5 chosen-prefix collisions\n"
		"Copyright (C) 2009 Marc Stevens\n"
		"http://homepages.cwi.nl/~stevens/\n"
		<< endl;

	try {
		birthday_parameters parameters;

		// Define program options
		po::options_description 
			desc("Allowed options"), 			
			all("Allowed options");

		desc.add_options()
			("help,h", "Show options.")
			("mod,m"
				, po::value<unsigned>(&parameters.modn)->default_value(1)
				, "Do only 1/m of all work.")

			("index,i"
				, po::value<unsigned>(&parameters.modi)->default_value(0)
				, "Do i-th part of all work (i=0,...,m-1).")

			("workdir,w"
				, po::value<string>(&workdir)->default_value("./data")
				, "Set working directory.")

			("inputfile1"
				, po::value<string>(&parameters.inputfile1)
				, "Use specified inputfile as first message.")

			("inputfile2"
				, po::value<string>(&parameters.inputfile2)
				, "Use specified inputfile as second message.")

			("outputfile1"
				, po::value<string>(&parameters.outputfile1)->default_value("file1.bin")
				, "Use specified outputfile for first message.")

			("outputfile2"
				, po::value<string>(&parameters.outputfile2)->default_value("file2.bin")
				, "Use specified outputfile for second message.")

			("distribution"
				, po::bool_switch(&parameters.distribution)
				, "Show workdistribution.\nDepends on hybridbits and pathtyperange.")

			("hybridbits"
				, po::value<unsigned>(&parameters.hybridbits)->default_value(0)
				, "Set 0 for 64-bit, 32 for 96-bit birthdaying.")

			("pathtyperange"
				, po::value<unsigned>(&parameters.pathtyperange)->default_value(2)
				, "Increases potential # diffs eliminated per n.c.")

			("maxblocks"
				, po::value<unsigned>(&parameters.maxblocks)->default_value(9)
				, "Upperbound on the amount of near-collisions.")

			("logtraillength"
				, po::value<int>(&parameters.logpathlength)->default_value(-1)
				, "Specify avg. trail length in bits.")

			("maxmemory"
				, po::value<unsigned>(&parameters.maxmemory)->default_value(100)
				, "Max. memory in MB used for storing trails.")

			("memhardlimit"
				, po::bool_switch(&parameters.memhardlimit)
				, "Hard limit max. memory instead of average.")
			("threads"
				, po::value<unsigned>(&parameters.threads)->default_value(0)
				, "Number of computing threads to start.")
			;

#ifdef HAVE_CUDA
		desc.add_options()
			("cuda_dev_query", "Query CUDA devices")
			("cuda_enable", po::bool_switch(&parameters.cuda_enabled), "Enable CUDA")
			;
#else
		all.add_options()			
			("cuda_dev_query", "Query CUDA devices")
			("cuda_enable", po::bool_switch(&parameters.cuda_enabled), "Enable CUDA")
			;
#endif

		all.add(desc);
	
		// Parse program options
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, all), vm);
		{
			std::ifstream ifs("md5birthdaysearch.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, all), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help")) {
			cout << desc << endl;
			return 0;
		}
		if (parameters.distribution) {
			determine_nrblocks_distribution(parameters);
			return 0;
		}

		if (parameters.modi >= parameters.modn) {
			cerr << "Error: i must be strictly less than m!" << endl;
			return 1;
		}
		if (parameters.modn != 1)
			cout << "System " << parameters.modi << " of " << parameters.modn << endl;
		if (parameters.inputfile1.size() == 0)
		{
			cerr << "Error: inputfile1 must be given!" << endl;
			return 2;
		}
		if (parameters.inputfile2.size() == 0)
		{
			cerr << "Error: inputfile2 must be given!" << endl;
			return 2;
		}
		if (parameters.logpathlength < 0 && parameters.maxmemory == 0)
		{
			cerr << "Cannot determine logtraillength when maxmemory is unspecified!" << endl;
			return 3;
		}
		if (parameters.modi == 0)
		{
			if (parameters.outputfile1.size() == 0)
			{
				cerr << "Error: outputfile1 must be given!" << endl;
				return 2;
			}
			if (parameters.outputfile2.size() == 0)
			{
				cerr << "Error: outputfile2 must be given!" << endl;
				return 2;
			}			
		}
#ifdef HAVE_CUDA
		if (vm.count("cuda_dev_query")) {
			cuda_device_query();
			return 0;
		}
#endif // CUDA
		result = dostep(parameters);
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
	return result;
}
