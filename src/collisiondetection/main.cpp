#include <stdexcept>
#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <hc/differentialpath.hpp>
#include <hc/booleanfunction.hpp>
#include <hc/md5detail.hpp>
#include <hc/timer.hpp>
#include <hc/sdr.hpp>
#include <hc/rng.hpp>

namespace tmp {
#define UINT32_DEFINED
	typedef boost::uint32_t uint32;
	extern "C" {
		#include "md5detectcoll.c"
	}
}

using namespace hc;
using namespace std;
namespace po = boost::program_options;
int param_tb = 60;

void show_diffpath(const vector<uint32>& ihv, const vector<uint32>& msg, const vector<uint32>& ihv2, const vector<uint32>& msg2)
{
	uint32 Q[68], Q2[68];
	Q[0] = ihv[0]; Q[1] = ihv[3]; Q[2] = ihv[2]; Q[3] = ihv[1];
	Q2[0] = ihv2[0]; Q2[1] = ihv2[3]; Q2[2] = ihv2[2]; Q2[3] = ihv2[1];
	for (unsigned t = 0; t < 64; ++t) {
		Q[3+t+1] = md5_step(t, Q[3+t], Q[3+t-1], Q[3+t-2], Q[3+t-3], msg[md5_wt[t]]);
		Q2[3+t+1] = md5_step(t, Q2[3+t], Q2[3+t-1], Q2[3+t-2], Q2[3+t-3], msg2[md5_wt[t]]);
	}
	uint32 mdiff[16];
	for (unsigned i = 0; i < 16; ++i) {
		mdiff[i] = msg2[i] - msg[i];
		if (mdiff[i]) {
			cout << "dm" << i << "=" << naf(mdiff[i]) << " ";
		}
	}
	cout << endl;
	differentialpath path;
	for (int t = -3; t <= 64; ++t)
		path[t] = sdr(Q[3+t], Q2[3+t]);

	bf_conditions cond;
	for (int t = 63; t >= 0; --t) {
		sdr deltaF;
		booleanfunction* F;
		if (t < 16) {
			F = &MD5_F_data;
			deltaF = sdr(md5_ff(Q[3+t], Q[3+t-1], Q[3+t-2]), md5_ff(Q2[3+t], Q2[3+t-1], Q2[3+t-2]));
		} else if (t < 32) {
			F = &MD5_G_data;
			deltaF = sdr(md5_gg(Q[3+t], Q[3+t-1], Q[3+t-2]), md5_gg(Q2[3+t], Q2[3+t-1], Q2[3+t-2]));
		} else if (t < 32) {
			F = &MD5_H_data;
			deltaF = sdr(md5_hh(Q[3+t], Q[3+t-1], Q[3+t-2]), md5_hh(Q2[3+t], Q2[3+t-1], Q2[3+t-2]));
		} else {
			F = &MD5_I_data;
			deltaF = sdr(md5_ii(Q[3+t], Q[3+t-1], Q[3+t-2]), md5_ii(Q2[3+t], Q2[3+t-1], Q2[3+t-2]));
		}
		for (unsigned b = 0; b < 32; ++b) {
			cond = F->backwardconditions(path(t,b), path(t-1,b), path(t-2,b), fromdiff(deltaF.get(b)));
			path.setbitcondition(t, b, cond.get<0>());
			path.setbitcondition(t-1, b, cond.get<1>());
			path.setbitcondition(t-2, b, cond.get<2>());
		}
	}
	show_path(path, mdiff);
	cout << " 9: ";
	for (int b = 31; b >= 0; --b)
		cout << path[9][b];
	cout << endl;
	cout << "10: ";
	for (int b = 31; b >= 0; --b)
		cout << ((Q[3+10]>>b)&1);
	cout << endl;
	cout << "11: ";
	for (int b = 31; b >= 0; --b)
		cout << ((Q[3+11]>>b)&1);
	cout << endl;
	return;
	
	uint32 QQ[68], QQ2[68], mm[16], mm2[16];
	uint64 cnt = 0, okcnt = 0;
	const int tb = param_tb; //60;

	// measure default wang prob
	while (okcnt < 256) {
		for (int i = 0; i < 4; ++i) {
			QQ[3+tb-i] = xrng64();
			QQ2[3+tb-i] = QQ[3+tb-i]  + Q2[3+tb-i]-Q[3+tb-i]; //(uint32(1)<<31);
		}
		for (unsigned i = 0; i < 16; ++i) {
			mm[i] = xrng64();
			mm2[i] = mm[i] + mdiff[i];
		}
		for (int t = tb; t < 64; ++t) {
			QQ[3+t+1] = md5_step(t, QQ[3+t], QQ[3+t-1], QQ[3+t-2], QQ[3+t-3], mm[md5_wt[t]]);
			QQ2[3+t+1] = md5_step(t, QQ2[3+t], QQ2[3+t-1], QQ2[3+t-2], QQ2[3+t-3], mm2[md5_wt[t]]);
		}
		bool ok = true;
		for (unsigned i = 68-4; i < 68; ++i) {
			if (naf(QQ2[i]-QQ[i]).mask & ~uint32((1<<31)|(1<<25))) 
				ok = false;
		}
		if (ok) ++okcnt;
		++cnt;
//		if (((++cnt) & 0xFFFFFF) == 0) cout << cnt << " " << okcnt << " 2^" << log(double(okcnt)/double(cnt))/log(2.0) << endl;
	}
	double wangprob = double(okcnt)/double(cnt);
	cout << "Prob: [" << tb << ",64):" << endl;
	cout << "- Wang: 2^" << log(double(okcnt)/double(cnt))/log(2.0) << endl;

	cnt = okcnt = 0;
	// measure this block prob
	while (okcnt < 16) {
		for (int i = 0; i < 4; ++i) {
			QQ[3+tb-i] = xrng64();
			QQ2[3+tb-i] = QQ[3+tb-i]  + Q2[3+tb-i]-Q[3+tb-i]; // (uint32(1)<<31);
		}
		for (unsigned i = 0; i < 16; ++i) {
			mm[i] = xrng64();
			mm2[i] = mm[i] + mdiff[i];
		}
		for (int t = tb; t < 64; ++t) {
			QQ[3+t+1] = md5_step(t, QQ[3+t], QQ[3+t-1], QQ[3+t-2], QQ[3+t-3], mm[md5_wt[t]]);
			QQ2[3+t+1] = md5_step(t, QQ2[3+t], QQ2[3+t-1], QQ2[3+t-2], QQ2[3+t-3], mm2[md5_wt[t]]);
		}
		bool ok = true;
		for (unsigned i = 68-4; i < 68; ++i) {
			if (Q2[i] - Q[i] != QQ2[i] - QQ[i])		
				ok = false;
		}
		if (ok) {
			++okcnt;
		}
//		++cnt;
		if (( (++cnt) & 0xFFFFF ) == 0) {
			cout << cnt << " " << okcnt << " 2^" << log(double(okcnt)/double(cnt))/log(2.0) << "\r";
		}
	}
	cout << endl;
	cout << "- This: " << cnt << " " << okcnt << " 2^" << log(double(okcnt)/double(cnt))/log(2.0) << endl;
	double thisprob = double(okcnt)/double(cnt);
	cout << "- Wang/prob: 2^" << log(wangprob/thisprob)/log(2.0) << endl;
	for (unsigned i = 0; i < 4; ++i)
		cout << hex << ihv2[i] << " ";
	cout << endl;
	for (unsigned i = 0; i < 4; ++i)
		cout << hex << md5_iv[i] << dec << " ";
	cout << endl; 
	for (unsigned i = 0; i < 16; ++i)
		cout << "m" << i << "=" << hex << msg2[i] << dec << " ";
	cout << endl;
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
				block[k] += uint32((unsigned char)(uc)) << (c*8);
			} else {
				i.putback(uc);
				i.setstate(ios::failbit);
				return len;
			}
		}
	return len;
}

inline unsigned load_block(ifstream& ifs, uint32 block[], uint32 pos)
{
	ifs.seekg(pos*64, ios::beg);
	return load_block(ifs, block);
}
inline unsigned load_block(istream& i, vector<uint32>& block)
{
	block.resize(16);
	return load_block(i, &block[0]);
}
inline unsigned load_block(ifstream& ifs, vector<uint32>& block, uint32 pos)
{
	ifs.seekg(pos, ios::beg);
	return load_block(ifs, block);
}
inline void md5compress(vector<uint32>& ihv, const vector<uint32>& msg) {
	md5compress(&ihv[0], &msg[0]);
}

uint32 getfilelen(ifstream& ifs) {
	uint32 pos = uint32(ifs.tellg());
	ifs.seekg(0, ios::end);
	uint32 len = uint32(ifs.tellg());
	ifs.seekg(pos, ios::beg);
	return len;
}

bool detectmd5coll(ifstream& ifs, bool verbose, bool diffpath, unsigned Boffset) {
	bool hascoll = false;
	vector<uint32> ihv2(4), block2(16);

	uint32 filelen = getfilelen(ifs);
	if (Boffset <= filelen)
		filelen -= Boffset;
	else
		filelen = 0;
	cout << "Filelen=" << filelen << " blocks=" << (filelen>>6) << flush;

	filelen >>= 6;
	vector< vector<uint32> > blocks(filelen, vector<uint32>(16));
	vector< vector<uint32> > ihvs(filelen+1, vector<uint32>(4));
	for (unsigned i = 0; i < 4; ++i)
		ihvs[0][i] = md5_iv[i];
	for (uint32 i = 0; i < filelen; ++i) {
		load_block(ifs, blocks[i], (i<<6) + Boffset);
		ihvs[i+1] = ihvs[i];
		md5compress(ihvs[i+1], blocks[i]);
		if (0) { //i == 0) {
			unsigned char* p = (unsigned char*)(&blocks[0][0]);
			for (unsigned int k = 0; k < 64; ++k)
				cout << (unsigned int)(p[k]) << " ";
			cout << endl;
		}
	}

	vector< vector<uint32> > ihvs2 = ihvs;
	vector< vector<uint32> > blocks2 = blocks;
	for (uint32 i = 0; i < filelen; ++i) {
		tmp::collision_type coll = tmp::detect_nearcollision(&ihvs[i][0], &blocks[i][0], &ihvs2[i][0], &blocks2[i][0], &ihvs[i+1][0]);
		if (coll == tmp::ct_dbb) {
			// if dbb is the only near-collision block then it is a false positive
			if (i == 0) 
				continue;
			tmp::collision_type coll2 = tmp::detect_nearcollision(&ihvs[i-1][0], &blocks[i-1][0], &ihvs2[i-1][0], &blocks2[i-1][0], &ihvs2[i][0]);
			if (coll2 == tmp::ct_none) 
				continue;
		}
		if (coll != tmp::ct_none) {
			hascoll = true;
			if (verbose) {
				cout << " [(block " << i << ": ";
				switch (coll) {
					case tmp::ct_dbb: cout << "dbb)" << flush; break;
					case tmp::ct_wang: cout << "wang)" << flush; break;
					case tmp::ct_cpc: cout << "cpc)" << flush; break;
					case tmp::ct_general: cout << "general)" << flush; break;
					case tmp::ct_special: cout << "special)" << flush; break;
				}
				if (diffpath && coll != tmp::ct_none) {
					cout << endl;
					show_diffpath(ihvs[i], blocks[i], ihvs2[i], blocks2[i]);
				}
				for (uint32 j = i; j > 0; --j) {
					coll = tmp::detect_nearcollision(&ihvs[j-1][0], &blocks[j-1][0], &ihvs2[j-1][0], &blocks2[j-1][0], &ihvs2[j][0]);
					switch (coll) {
						case tmp::ct_dbb: cout << "(block " << j-1 << ": dbb)" << flush; break;
						case tmp::ct_wang: cout << "(block " << j-1 << ": wang)" << flush; break;
						case tmp::ct_cpc: cout << "(block " << j-1 << ": cpc)" << flush; break;
						case tmp::ct_general: cout << "(block " << j-1 << ": general)" << flush; break;
						case tmp::ct_special: cout << "(block " << j-1 << ": special)" << flush; break;
					}
					if (diffpath && coll != tmp::ct_none) {
						cout << endl;
						show_diffpath(ihvs[j-1], blocks[j-1], ihvs2[j-1], blocks2[j-1]);
					}
					if (coll == tmp::ct_none || j-1==0) {
						cout << "]" << flush;
						break;
					}
				}
			}
		}
	}
	if (!verbose && hascoll)
		cout << "possible collision detected." << endl;
	else if (hascoll) cout << endl;
	return !hascoll;
}

int main(int argc, char** argv) 
{
	hc::timer runtime(true);
	vector<string> files;
	bool verbose = false, diffpath = false;
	unsigned Boffset = 0;
	try {
		// Define program options
		po::options_description 
			desc("Allowed options");

		desc.add_options()
			("help,h",		"Show options")
			("verbose,v", po::bool_switch(&verbose), "Be verbose")
			("diffpath,p", po::bool_switch(&diffpath), "Reconstruct diff.paths")
			("file,f", po::value<vector<string> >(&files), "File(s) to be checked")
			("offset,o", po::value<unsigned>(&Boffset)->default_value(0), "Byte offset")
			("tb,t", po::value<int>(&param_tb)->default_value(32), "t_begin")
			;
	
		// Parse program options
		po::positional_options_description p;
		p.add("file", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv)
			.options(desc).positional(p).run(), vm);
		{
			std::ifstream ifs("detectmd5coll.cfg");
			if (ifs) po::store(po::parse_config_file(ifs, desc), vm);
		}
		po::notify(vm);

		// Process program options
		if (vm.count("help") || vm.count("file") == 0) {
			cout << desc << endl;
			return 0;
		}

		for (unsigned i = 0; i < files.size(); ++i) {
			cout << "Checking file '" << files[i] << "': " << flush;
			try {
				ifstream ifs(files[i].c_str(), ios::binary);
				if (!ifs) {
					cout << "could not open file!" << endl;
					continue;
				}
				if (detectmd5coll(ifs, verbose, diffpath, Boffset))
					cout << " no collisions detected." << endl;
			} catch (exception& e) {
				cout << e.what() << endl;
			}
		}

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
