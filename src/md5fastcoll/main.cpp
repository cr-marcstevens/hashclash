/*

MD5 collision generator
=======================
Source code files:
  block0.cpp
  block1.cpp
  main.cpp
  main.hpp
  md5.cpp
  block1wang.cpp
  block1stevens00.cpp
  block1stevens01.cpp
  block1stevens10.cpp
  block1stevens11.cpp
Win32 executable:
  fastcoll_v1.0.0.5.exe

Version
=======
version 1.0.0.5-1, April 2006.

Copyright
=========
© M. Stevens, 2006. All rights reserved.

Disclaimer
==========
This software is provided as is. Use is at the user's risk.
No guarantee whatsoever is given on how it may function or malfunction.
Support cannot be expected.
This software is meant for scientific and educational purposes only.
It is forbidden to use it for other than scientific or educational purposes.
In particular, commercial and malicious use is not allowed.
Further distribution of this software, by whatever means, is not allowed
without our consent.
This includes publication of source code or executables in printed form,
on websites, newsgroups, CD-ROM's, etc.
Changing the (source) code without our consent is not allowed.
In all versions of the source code this disclaimer, the copyright
notice and the version number should be present.

*/

#include <iostream>
#include <fstream>
#include <time.h>

#include "main.hpp"

using namespace std;

const uint32 MD5IV[] = { 0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476 };

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, const uint32 block[]);
void find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[], uint32 msg2block0[], uint32 msg2block1[], bool verbose = false);

#if 0

// example trivial version with md5 initial value
int main()
{
	seed32_1 = uint32(time(NULL));
	seed32_2 = 0x12345678;

	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];
	find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, true);

	ofstream ofs1("msg1", ios::binary);
	save_block(ofs1, msg1block0);
	save_block(ofs1, msg1block1);
	ofstream ofs2("msg2", ios::binary);
	save_block(ofs2, msg2block0);
	save_block(ofs2, msg2block1);
}

#else

#include <sstream>
#include <string>
#include <utility>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <boost/cstdint.hpp>

typedef boost::uint64_t uint64;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
using boost::lexical_cast;

void test_md5iv(bool single = false);
void test_rndiv(bool single = false);
void test_reciv(bool single = false);
void test_all();

int main(int argc, char** argv)
{
	seed32_1 = uint32(time(NULL));
	seed32_2 = 0x12345678;

	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };

	string outfn1 = "msg1.bin";
	string outfn2 = "msg2.bin";
	string ihv;
	string prefixfn;
	bool verbose = true;

	cout <<
		"MD5 collision generator v1.5\n"
		"by Marc Stevens (http://www.win.tue.nl/hashclash/)\n"
		<< endl;

	try
	{
		boost::timer runtime;

		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "Show options.")
			("quiet,q", "Be less verbose.")
			("ihv,i", po::value<string>(&ihv), "Use specified initial value. Default is MD5 initial value.")
			("prefixfile,p", po::value<string>(&prefixfn), "Calculate initial value using given prefixfile. Also copies data to output files.")			
			("out,o", po::value<vector<string> >()->multitoken(), "Set output filenames. This must be the last option and exactly 2 filenames must be specified. \nDefault: -o msg1.bin msg2.bin")
			;

		po::options_description hidden;
		hidden.add_options()
			("testmd5iv", "Endlessly time collision generation using MD5 initial value.")
			("testrndiv", "Endlessly time collision generation using arbitrary random initial values.")
			("testreciv", "Endlessly time collision generation using recommended random initial values.")
			("testall", "Endlessly time collision generation for each case.")
			;

		po::options_description cmdline;
		cmdline.add(desc).add(hidden);
		po::positional_options_description p;
		p.add("prefixfile", 1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(cmdline).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("quiet"))
			verbose = false;

		if (vm.count("help") || argc == 1) {
			cout << desc << endl;
			return 1;
		}

		if (vm.count("testmd5iv"))
			test_md5iv();
		if (vm.count("testrndiv"))
			test_rndiv();
		if (vm.count("testreciv"))
			test_reciv();
		if (vm.count("testall"))
			test_all();

		if (vm.count("prefixfile"))
		{
			unsigned l = prefixfn.size();
			if (l >= 4 && prefixfn[l-4]=='.' && prefixfn[l-3]!='.' && prefixfn[l-2]!='.' && prefixfn[l-1]!='.')
			{
				outfn1 = prefixfn.substr(0, l-4) + "_msg1" + prefixfn.substr(l-4);
				outfn2 = prefixfn.substr(0, l-4) + "_msg2" + prefixfn.substr(l-4);
				unsigned i = 1;
				while ( fs::exists(fs::path(outfn1)) 
					 || fs::exists(fs::path(outfn2)))
				{
					outfn1 = prefixfn.substr(0, l-4) + "_msg1_" + lexical_cast<string>(i) + prefixfn.substr(l-4);
					outfn2 = prefixfn.substr(0, l-4) + "_msg2_" + lexical_cast<string>(i) + prefixfn.substr(l-4);
					++i;
				}
			}
		}

		if (vm.count("out"))
		{
			vector<string> outfns = vm["out"].as< vector<string> >();
			if (outfns.size() != 2)
			{
				cerr << "Error: exactly two output filenames should be specified." << endl;
				return 1;
			}
			outfn1 = outfns[0];
			outfn2 = outfns[1];
		}

		if (verbose)
			cout << "Using output filenames: '" << outfn1 << "' and '" << outfn2 << "'" << endl;
		ofstream ofs1(outfn1.c_str(), ios::binary);
		if (!ofs1)
		{
			cerr << "Error: cannot open outfile: '" << outfn1 << "'" << endl;
			return 1;
		}
		ofstream ofs2(outfn2.c_str(), ios::binary);
		if (!ofs2)
		{
			cerr << "Error: cannot open outfile: '" << outfn2 << "'" << endl;
			return 1;
		}

		if (vm.count("prefixfile"))
		{
			if (verbose)
				cout << "Using prefixfile: '" << prefixfn << "'" << endl;
			ifstream ifs(prefixfn.c_str(), ios::binary);
			if (!ifs)
			{
				cerr << "Error: cannot open inputfile: '" << prefixfn << "'" << endl;
				return 1;
			}
			uint32 block[16];
			while (true)
			{
				unsigned len = load_block(ifs, block);
				if (len)
				{
					save_block(ofs1, block);
					save_block(ofs2, block);
					md5_compress(IV, block);
				} else
					break;
			}
		}
		else
		{
			if (vm.count("ihv") == 0)
				ihv = "0123456789abcdeffedcba9876543210";
			if (ihv.size() != 32)
			{
				cerr << "Error: an initial value must be specified as a hash value of 32 hexadecimal characters." << endl;
				return 1;
			} else
			{
				uint32 c;
				for (unsigned i = 0; i < 4; ++i)
				{
					IV[i] = 0;
					for (unsigned b = 0; b < 4; ++b)
					{
						stringstream ss;
						ss << ihv.substr(i*8+b*2,2);
						ss >> hex >> c;					
						IV[i] += c << (b*8);
					}
				}

			}
		}
		if (verbose)
		{
			cout << "Using initial value: " << hex;
			unsigned oldwidth = cout.width(2);
			char oldfill = cout.fill('0');
			
			for (unsigned i = 0; i < 4; ++i)
			{
				for (unsigned b = 0; b < 4; ++b)
				{
					cout.width(2);
					cout.fill('0');
					cout << ((IV[i]>>(b*8))&0xFF);
				}
			}
			cout.width(oldwidth);
			cout.fill(oldfill);
			cout << dec << endl;
		}

		if (verbose)
			cout << endl;

		uint32 msg1block0[16];
		uint32 msg1block1[16];
		uint32 msg2block0[16];
		uint32 msg2block1[16];
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, true);

		save_block(ofs1, msg1block0);
		save_block(ofs1, msg1block1);
		save_block(ofs2, msg2block0);
		save_block(ofs2, msg2block1);
		if (verbose)
			cout << "Running time: " << runtime.elapsed() << " s" << endl;
		return 0;
	} catch (exception& e)
	{
		cerr << "\nException caught:\n" << e.what() << endl;
		return 1;
	} catch (...)
	{
		cerr << "\nUnknown exception caught!" << endl;
		return 1;
	}
}

void test_md5iv(bool single)
{
	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	boost::timer runtime;
	while (true)
	{
		runtime.restart();
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.elapsed();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_md5iv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_rndiv(bool single)
{
	uint32 IV[4];
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	boost::timer runtime;
	while (true)
	{
		runtime.restart();
		IV[0] = xrng64(); IV[1] = xrng64(); IV[2] = xrng64(); IV[3] = xrng64();
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.elapsed();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_rndiv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_reciv(bool single)
{
	uint32 IV[4];
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	boost::timer runtime;
	while (true)
	{
		runtime.restart();
		IV[0] = xrng64(); IV[1] = xrng64(); IV[2] = xrng64(); IV[3] = xrng64();
		IV[2] |= 1<<25; IV[2] ^= ((IV[2] & (1<<24))<<1);
		IV[3] &= ~(1<<25); IV[3] ^= ((IV[3] & (1<<24))<<1);

		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.elapsed();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_reciv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_all()
{
	while (true)
	{
		test_md5iv(true);
		test_reciv(true);
		test_rndiv(true);
	}
}

#endif


unsigned load_block(istream& i, uint32 block[])
{
	unsigned len = 0;
	char uc;
	for (unsigned k = 0; k < 16; ++k)
	{
		block[k] = 0;
		for (unsigned c = 0; c < 4; ++c)
		{
			i.get(uc);
			if (i) 
				++len;
			else
				uc = 0;
			block[k] += uint32((unsigned char)(uc))<<(c*8);
		}
	}
	return len;
}

void save_block(ostream& o, const uint32 block[])
{
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
			o << (unsigned char)((block[k] >> (c*8))&0xFF);
}

void find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[], uint32 msg2block0[], uint32 msg2block1[], bool verbose)
{
	if (verbose)
		cout << "Generating first block: " << flush;
	find_block0(msg1block0, IV);

	uint32 IHV[4] = { IV[0], IV[1], IV[2], IV[3] };
	md5_compress(IHV, msg1block0);

	if (verbose)
		cout << endl << "Generating second block: " << flush;
	find_block1(msg1block1, IHV);

	for (int t = 0; t < 16; ++t)
	{
		msg2block0[t] = msg1block0[t];
		msg2block1[t] = msg1block1[t];
	}
	msg2block0[4] += 1 << 31; msg2block0[11] += 1 << 15; msg2block0[14] += 1 << 31;
	msg2block1[4] += 1 << 31; msg2block1[11] -= 1 << 15; msg2block1[14] += 1 << 31;
	if (verbose)
		cout << endl;
}
