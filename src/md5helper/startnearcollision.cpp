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

#include <string>
#include <stdexcept>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>

#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>

#include "main.hpp"

#define HASHCLASH_MD5COMPRESS_STEP(f, a, b, c, d, m, ac, rc) \
        a += f(b, c, d) + m + ac; a = rotate_left(a,rc); a += b;

double testnearcollprob(uint32 dm11, uint32 diffcd, uint32 diffb) {
        uint32 count = 0;
        for (unsigned k = 0; k < (1<<23); ++k)
        {
                uint32 a = xrng64(), b = xrng64()+xrng64()*11, c = xrng64(), d = xrng64()+xrng64()*11;
                uint32 a2 = a, b2 = b, c2 = c, d2 = d;
                uint32 m11 = xrng64(), m2 = xrng64(), m9 = xrng64();
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, m11, 0xbd3af235, 10);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, m2, 0x2ad7d2bb, 15);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, m9, 0xeb86d391, 21);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d2, a2, b2, c2, m11+dm11, 0xbd3af235, 10);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c2, d2, a2, b2, m2, 0x2ad7d2bb, 15);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b2, c2, d2, a2, m9, 0xeb86d391, 21);
                if (d2-d==diffcd && c2-c==diffcd && b2-b==diffb)
                        ++count;
        }
        return double(count)/double(1<<23);
}

void constructupperpath_sbcpc(differentialpath& path, uint32 dm[16], uint32 diffihv[4])
{
	uint32 a,b,c,d,a2,b2,c2,d2,oa,ob,oc,od,m11,m2,m9,m4;
	a = xrng64(); b = xrng64()+xrng64()*11; c = xrng64(); d = xrng64()+xrng64()*11;
	while (true) {
                oa = a2 = a; ob = b2 = b; oc = c2 = c; od = d2 = d;
                m4 = xrng64();
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, m4, 0xf7537e82,  6);
		HASHCLASH_MD5COMPRESS_STEP(md5_ii, a2, b2, c2, d2, m4+dm[4], 0xf7537e82,  6);
		if (a2-a != diffihv[0]) {
			b = xrng64(); c = xrng64(); d = xrng64();
			continue;
		}
		m11 = xrng64(); 
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, m11, 0xbd3af235, 10);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d2, a2, b2, c2, m11+dm[11], 0xbd3af235, 10);
		if (d2-d != diffihv[3]) {
			b = xrng64(); c = xrng64();
			continue;
		}
		m2 = xrng64(); 
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, m2, 0x2ad7d2bb, 15);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c2, d2, a2, b2, m2+dm[2], 0x2ad7d2bb, 15);
		if (c2-c != diffihv[2]) {
			b = xrng64();
			continue;
		}
		m9 = xrng64(); 
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, m9, 0xeb86d391, 21);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b2, c2, d2, a2, m9+dm[9], 0xeb86d391, 21);
		if (b2-b == diffihv[1]) break;
        }
	for (unsigned i = 55; i <= 64; ++i)
		path[i].clear();
	path[61] = sdr(a,a2);
	path[62] = sdr(d,d2);
	path[63] = sdr(c,c2);
	path[64] = sdr(b,b2);

	a2 = oa; b2 = ob; c2 = oc; d2 = od;
	a = oa; b = ob; c = oc; d = od;
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, a, b, c, d, m4, 0xf7537e82,  6);
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, a2, b2, c2, d2, m4+dm[4], 0xf7537e82,  6);
	sdr dF61 = sdr(md5_ii(a,b,c),md5_ii(a2,b2,c2));
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, m11, 0xbd3af235, 10);
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, d2, a2, b2, c2, m11+dm[11], 0xbd3af235, 10);
	sdr dF62 = sdr(md5_ii(d,a,b),md5_ii(d2,a2,b2));
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, m2, 0x2ad7d2bb, 15);
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, c2, d2, a2, b2, m2+dm[2], 0x2ad7d2bb, 15);
	sdr dF63 = sdr(md5_ii(c,d,a),md5_ii(c2,d2,a2));

	for (unsigned b = 0; b < 32; ++b) {
		bitcondition f61 = bc_constant, f62 = bc_constant, f63 = bc_constant;
		if (dF61.get(b) == +1) f61 = bc_plus;
		if (dF61.get(b) == -1) f61 = bc_minus;
		if (dF62.get(b) == +1) f62 = bc_plus;
		if (dF62.get(b) == -1) f62 = bc_minus;
		if (dF63.get(b) == +1) f63 = bc_plus;
		if (dF63.get(b) == -1) f63 = bc_minus;
		bf_conditions bc0 = MD5_I_data.forwardconditions(path[61][b], bc_constant, bc_constant, f61);
		bf_conditions bc1 = MD5_I_data.forwardconditions(path[62][b], bc0.first, bc0.second, f62);
		bf_conditions bc2 = MD5_I_data.forwardconditions(path[63][b], bc1.first, bc1.second, f63);
		path.setbitcondition(59, b, bc0.third);
		path.setbitcondition(60, b, bc1.third);
		path.setbitcondition(61, b, bc2.third);
		path.setbitcondition(62, b, bc2.second);
		path.setbitcondition(63, b, bc2.first);
	}
}

void constructupperpath(differentialpath& path, uint32 dm11, uint32 diffcd, uint32 diffb)
{
	uint32 a,b,c,d,a2,b2,c2,d2,oa,ob,oc,od,m11,m2,m9;
	while (true) {
                a = xrng64(); b = xrng64()+xrng64()*11; c = xrng64(); d = xrng64()+xrng64()*11;
                a2 = a; b2 = b; c2 = c; d2 = d;
                oa = a; ob = b; oc = c; od = d;
                m11 = xrng64(); m2 = xrng64(); m9 = xrng64();
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, m11, 0xbd3af235, 10);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, m2, 0x2ad7d2bb, 15);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b, c, d, a, m9, 0xeb86d391, 21);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, d2, a2, b2, c2, m11+dm11, 0xbd3af235, 10);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, c2, d2, a2, b2, m2, 0x2ad7d2bb, 15);
                HASHCLASH_MD5COMPRESS_STEP(md5_ii, b2, c2, d2, a2, m9, 0xeb86d391, 21);
                if (d2-d==diffcd && c2-c==diffcd && b2-b==diffb)
                	break;
        }
	for (unsigned i = 58; i <= 64; ++i)
		path[i].clear();
	path[62] = sdr(d,d2);
	path[63] = sdr(c,c2);
	path[64] = sdr(b,b2);

	a2 = oa; b2 = ob; c2 = oc; d2 = od;
	a = oa; b = ob; c = oc; d = od;
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, d, a, b, c, m11, 0xbd3af235, 10);
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, d2, a2, b2, c2, m11+dm11, 0xbd3af235, 10);
	sdr dF62 = sdr(md5_ii(d,a,b),md5_ii(d2,a2,b2));
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, c, d, a, b, m2, 0x2ad7d2bb, 15);
	HASHCLASH_MD5COMPRESS_STEP(md5_ii, c2, d2, a2, b2, m2, 0x2ad7d2bb, 15);
	sdr dF63 = sdr(md5_ii(c,d,a),md5_ii(c2,d2,a2));

	for (unsigned b = 0; b < 32; ++b) {
		bitcondition f62 = bc_constant, f63 = bc_constant;
		if (dF62.get(b) == +1) f62 = bc_plus;
		if (dF62.get(b) == -1) f62 = bc_minus;
		if (dF63.get(b) == +1) f63 = bc_plus;
		if (dF63.get(b) == -1) f63 = bc_minus;
		bf_conditions bc1 = MD5_I_data.forwardconditions(path[62][b], bc_constant, bc_constant, f62);
		bf_conditions bc2 = MD5_I_data.forwardconditions(path[63][b], bc1.first, bc1.second, f63);
		path.setbitcondition(60, b, bc1.third);
		path.setbitcondition(61, b, bc2.third);
		path.setbitcondition(62, b, bc2.second);
		path.setbitcondition(63, b, bc2.first);
	}
}

int startnearcollision(parameters_type& parameters)
{
	uint32 ihv1[4];
	uint32 ihv2[4];
	uint32 msg1[16];
	uint32 msg2[16];
	{
		ifstream if1(parameters.infile1.c_str(), ios::binary);
		if (!if1) {
			cerr << "Error: cannot open inputfile 1 '" << parameters.infile1 << "'!" << endl;
			return 1;
		}
		ifstream if2(parameters.infile2.c_str(), ios::binary);
		if (!if2) {
			cerr << "Error: cannot open inputfile 2 '" << parameters.infile2 << "'!" << endl;
			return 1;
		}
		
		for (unsigned k = 0; k < 4; ++k)
			ihv1[k] = ihv2[k] = md5_iv[k];

		// load, md5 and save inputfile1
		unsigned file1blocks = 0;
		while (load_block(if1, msg1) > 0) {
			md5compress(ihv1, msg1);
			++file1blocks;
		}
		// load, md5 and save inputfile2
		unsigned file2blocks = 0;
		while (load_block(if2, msg2) > 0) {
			md5compress(ihv2, msg2);			
			++file2blocks;
		}
		if (file1blocks != file2blocks) {
			cerr << "Error: inputfile 1 and 2 are not of equal size" << endl;
			return 2;
		}
	}

	cout << "IHV1   = {" << ihv1[0] << "," << ihv1[1] << "," << ihv1[2] << "," << ihv1[3] << "}" << endl;
	cout << "IHV1   = " << hex;
	for (unsigned k = 0; k < 4; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout.width(2); cout.fill('0');
			cout << ((ihv1[k]>>(c*8))&0xFF);
		}
	cout << dec << endl << endl;

	cout << "IHV2   = {" << ihv2[0] << "," << ihv2[1] << "," << ihv2[2] << "," << ihv2[3] << "}" << endl;
	cout << "IHV2   = " << hex;
	for (unsigned k = 0; k < 4; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout.width(2); cout.fill('0');
			cout << ((ihv2[k]>>(c*8))&0xFF);
		}
	cout << dec << endl << endl;

	uint32 dihv[4] = { ihv2[0] - ihv1[0], ihv2[1] - ihv1[1], ihv2[2] - ihv1[2], ihv2[3] - ihv1[3] };
	cout << "dIHV   = {" << dihv[0] << "," << dihv[1] << "," << dihv[2] << "," << dihv[3] << "}" << endl;

	bool sbcpc = false;
	if (dihv[0]==0 && dihv[1]==0 && dihv[2]==0 && dihv[3]==0)
	{
		cerr << "dIHV is zero!" << endl;
		return 4;
	}
	if (dihv[0] == +(1<<5) && dihv[3] == +(1<<5) -(1<<25))
		sbcpc = true;
	else if (dihv[0] != 0 || dihv[2] != dihv[3])
	{
		cerr << "Error: dIHV is not of the required form dIHV[0]=0, dIHV[2]=dIHV[3]" << endl;
		return 3;
	}
	uint32 m_diff[16] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
	vector<uint32> dm11;
	vector< pair<uint32,uint32> > ncdiff;
	unsigned j = 0;
	if (sbcpc) {
		m_diff[11] = 1<<15;
		m_diff[4] = 1<<31;
		m_diff[14] = 1<<31;
		m_diff[2] = 1<<8;
	} else {
		sdr diffcd = naf(0-dihv[3]);
		sdr diffb  = naf(dihv[3] - dihv[1]);		
		unsigned bb = 0;
		while (bb < 32) {
			if (diffcd.get(bb) != 0) {
				unsigned bend=bb+1;
				while (bend < 32 && bend<=bb+parameters.pathtyperange && diffcd.get(bend)==0)
					++bend;
				uint32 diff1 = 1<<bb;
				if (diffcd.get(bb) == -1)
					diff1 = -(1<<bb);
				uint32 diff2 = diff1;
				for (unsigned b2 = bb; b2 < bend; ++b2)
					if (diffb.get((b2+21)%32) == +1)
						diff2 += 1<<((b2+21)%32);
					else if (diffb.get((b2+21)%32) == -1)
						diff2 -= 1<<((b2+21)%32);
				ncdiff.push_back( pair<uint32,uint32>(diff1,diff2) );
				bb = bend;
			} else if (diffb.get((bb+21)%32) != 0) {
				uint32 diff1 = 0, diff2 = 0;
				unsigned bend=bb+1;
				unsigned tempinc=0;
				if (bb+1 < 32 && diffcd.get(bb+1) != 0) {
					if (diffcd.get(bb+1) == +1)
						diff1 += 1<<bb;
					else
						diff1 -= 1<<bb;
					ncdiff.push_back( pair<uint32,uint32>(diff1,diff1) );
					bend=bb+2;
					tempinc=1;				
				} else {
					diff1 = 1<<bb;
					if (diffb.get((bb+21)%32) == +1)
						diff1 = 0-(1<<bb);
					ncdiff.push_back( pair<uint32,uint32>(0-diff1,0-diff1) );
				}
				diff2 = diff1;
				while (bend < 32 && bend<=bb+parameters.pathtyperange+tempinc && diffcd.get(bend)==0)
					++bend;
				for (unsigned b2 = bb; b2 < bend; ++b2)
					if (diffb.get((b2+21)%32) == +1)
						diff2 += 1<<((b2+21)%32);
					else if (diffb.get((b2+21)%32) == -1)
						diff2 -= 1<<((b2+21)%32);
				ncdiff.push_back( pair<uint32,uint32>(diff1,diff2) );
				bb = bend;
			} else ++bb;
		}
		if (ncdiff.size()> 1 && 0!=naf(ncdiff[ncdiff.size()-1].first).get(31) && 0!=(ncdiff[0].first&1)) {
			ncdiff[0].first += ncdiff[ncdiff.size()-1].first;
			ncdiff[0].second += ncdiff[ncdiff.size()-1].second;
			ncdiff.pop_back();
		}
		uint32 temp1 = 0, temp2 = 0;
		for (unsigned i = 0; i < ncdiff.size(); ++i)
		{
			cout << "NC" << i << ": delta c = delta d = " << naf(ncdiff[i].first) << endl;
			cout << "NC" << i << ": delta b           = " << naf(ncdiff[i].second) << endl;
			temp1 += ncdiff[i].first;
			temp2 += ncdiff[i].second;
		}
		if (dihv[3] != uint32(0-temp1) || dihv[1] != uint32(0-temp2)) {
			cerr << "Internal error in determining near-collision blocks:" << endl;
			cerr << dihv[3] << "=?" << uint32(0-temp1) << endl;
			cerr << dihv[1] << "=?" << uint32(0-temp2) << endl;
		}
		dm11.resize(ncdiff.size(),0);
		vector<double> prob(ncdiff.size(),0);
		for (unsigned i = 0; i < ncdiff.size(); ++i)
		{
			sdr fnaf = naf(ncdiff[i].first);
			for (unsigned b = 0; b < 32; ++b)
				if (fnaf.get(b) == +1) {
					dm11[i] = uint32(1)<<((b-10)%32);
					break;
				} else if (fnaf.get(b) == -1) {
					dm11[i] = 0-(uint32(1)<<((b-10)%32));
					break;
				}
			if (fnaf.get(0) != 0 && fnaf.get(31) != 0) {
				if (fnaf.get(0) == +1)
					dm11[i] = uint32(1)<<21;
				else
					dm11[i] = 0-(uint32(1)<<21);
			}
			prob[i] = testnearcollprob(dm11[i], ncdiff[i].first, ncdiff[i].second);
			cout << "NC" << i << ": prob=" << log(prob[i])/log(double(2)) << endl;
		}
		j = 0;
		for (unsigned i = 1; i < prob.size(); ++i)
			if (prob[i] > prob[j]) j = i;
		unsigned skipnum = parameters.skipnc % ncdiff.size();
		for (unsigned skip = 0; skip < skipnum; ++skip) {
			swap(dm11[j],dm11[dm11.size()-1]); dm11.pop_back();
			swap(ncdiff[j],ncdiff[ncdiff.size()-1]); ncdiff.pop_back();
			swap(prob[j],prob[prob.size()-1]); prob.pop_back();
			j = 0;
			for (unsigned i = 1; i < prob.size(); ++i)
				if (prob[i] > prob[j]) j = i;
		}
		unsigned b = 0;
		int m11sign = +1;
		for (bb = 0; bb < 32; ++bb)
			if (naf(dm11[j]).get(bb) != 0) {
				b = bb;
				if (naf(dm11[j]).get(bb) == -1) 
					m11sign = -1;
				break;
			}

		cout << "delta m_11 = " << naf(dm11[j]) << endl;
		m_diff[11] = dm11[j];
	}

	/*** Construct lower diff. path ***/
	differentialpath lowerpath;
	uint32 Q1[4] = { ihv1[0], ihv1[3], ihv1[2], ihv1[1] };
	uint32 Q2[4] = { ihv2[0], ihv2[3], ihv2[2], ihv2[1] };
	for (int i = 0; i < 4; ++i)
		for (unsigned k = 0; k < 32; ++k)
		{
			if ((Q1[i]>>k)&1) {
				if ((Q2[i]>>k)&1)
					lowerpath.setbitcondition(i-3,k,bc_one);
				else
					lowerpath.setbitcondition(i-3,k,bc_minus);
			} else {
				if ((Q2[i]>>k)&1)
					lowerpath.setbitcondition(i-3,k,bc_plus);
				else
					lowerpath.setbitcondition(i-3,k,bc_zero);
			}
		}
	uint32 dF = 0;
	for (unsigned k = 0; k < 32; ++k)
	{
		bf_outcome outcome = MD5_F_data.outcome(lowerpath(0,k), lowerpath(-1,k), lowerpath(-2,k));
		if (outcome.size()) {
			if (outcome[0] == bc_plus) 			dF += 1<<k;
			else if (outcome[0] == bc_minus)	dF -= 1<<k;
		}
	}
	uint32 dQtm3 = lowerpath[-3].diff();
	uint32 dT = dQtm3 + dF + m_diff[0];
	std::vector<std::pair<uint32,double> > rotateddiff;
	rotate_difference(dT, md5_rc[0], rotateddiff);
	double bestrot = 0;
	for (unsigned i = 0; i < rotateddiff.size(); ++i)
	{
		uint32 dR = rotateddiff[i].first;
		wordconditions Q1 = naf(lowerpath[0].diff() + dR);
		rotateddiff[i].second = check_rotation(dR, dT, md5_rc[0], lowerpath[0], Q1, 1<<12);
		if (rotateddiff[i].second > bestrot) {
			lowerpath[1] = Q1;
			bestrot = rotateddiff[i].second;
		}
	}

	/*** Show and store lower diff. path ***/
	show_path(lowerpath, m_diff);
	try {
		vector<differentialpath> temp;
		temp.push_back(lowerpath);
		save_gz(temp, workdir +  "/lowerpath", binary_archive);
		cout << "Saved lower diff. path to '" << workdir + "/lowerpath.bin.gz'." << endl;
	} catch (...) {
		cerr << "Error: could not write '" << workdir + "/lowerpath.bin.gz'!" << endl;
	}
	cout << endl;
	
	/*** Construct upper diff. path ***/
	differentialpath upperpath, bestpath;
	unsigned bestcond = 1<<31;
	for (unsigned i = 0; i < 1<<10; ++i) {
		if (sbcpc) {
			uint32 dihvmin[4] = { -dihv[0], -dihv[1], -dihv[2], -dihv[3] };
			constructupperpath_sbcpc(upperpath, m_diff, dihvmin);
		} else {
			upperpath[32].clear();
			constructupperpath(upperpath, dm11[j], ncdiff[j].first, ncdiff[j].second);
		}
		if (upperpath.nrcond() < bestcond) {
			bestcond = upperpath.nrcond();
			bestpath = upperpath;
		}
	}
	upperpath = bestpath;

	/*** Show and store upper diff. path ***/
	show_path(upperpath, m_diff);
	try {
		vector<differentialpath> temp;
		temp.push_back(upperpath);
		save_gz(temp, workdir +  "/upperpath", binary_archive);
		cout << "Saved upper diff. path to '" << workdir + "/upperpath.bin.gz'." << endl;
	} catch (...) {
		cerr << "Error: could not write '" << workdir + "/upperpath.bin.gz'!" << endl;
	}

	/*** Write md5diffpath_forward.cfg ***/
	{
		ofstream off("md5diffpathforward.cfg");
		if (!off)
			cerr << "Error: could not write md5diffpathforward.cfg!" << endl;
		else {
			/* write IHV info */
			off << "# IHV" << endl;
			off << "ihv0a = " << ihv1[0] << endl;
			off << "ihv1a = " << ihv1[1] << endl;
			off << "ihv2a = " << ihv1[2] << endl;
			off << "ihv3a = " << ihv1[3] << endl << endl;
			off << "# IHV'" << endl;
			off << "ihv0b = " << ihv2[0] << endl;
			off << "ihv1b = " << ihv2[1] << endl;
			off << "ihv2b = " << ihv2[2] << endl;
			off << "ihv3b = " << ihv2[3] << endl << endl;
			/* write delta m11 */
			off << "# message block difference" << endl;
			for (unsigned i = 0; i < 16; ++i)
				for (unsigned b = 0; b < 32; ++b)
					if (naf(m_diff[i]).get(b) != 0)
						off << "diffm" << i << " = " << int(b+1)*naf(m_diff[i]).get(b) << endl;
			/* other settings */
			off << "# parameters" << endl;

			/* copy template file if any */
			ifstream iff("md5diffpathforward.cfg.template");
			if (!iff)
				cerr << "Warning: could not read md5diffpathforward.cfg.template" << endl;
			else
				off << iff.rdbuf() << endl;
			cout << "Saved 'md5diffpathforward.cfg'." << endl;
		}
	}

	/*** Write md5diffpath_backward.cfg ***/
	{
		ofstream ofb("md5diffpathbackward.cfg");
		if (!ofb)
			cerr << "Error: could not write md5diffpathbackward.cfg!" << endl;
		else {
			/* write delta m11 */
			ofb << "# message block difference" << endl;
			for (unsigned i = 0; i < 16; ++i)
				for (unsigned b = 0; b < 32; ++b)
					if (naf(m_diff[i]).get(b) != 0)
						ofb << "diffm" << i << " = " << int(b+1)*naf(m_diff[i]).get(b) << endl;
			/* other settings */
			ofb << "# parameters" << endl;

			/* copy template file if any */
			ifstream ifb("md5diffpathbackward.cfg.template");
			if (!ifb)
				cerr << "Warning: could not read md5diffpathbackward.cfg.template" << endl;
			else
				ofb << ifb.rdbuf() << endl;
			cout << "Saved 'md5diffpathbackward.cfg'." << endl;
		}
	}

	/*** Write md5diffpath_connect.cfg ***/
	{
		ofstream ofc("md5diffpathconnect.cfg");
		if (!ofc)
			cerr << "Error: could not write md5diffpathconnect.cfg!" << endl;
		else {
			/* write delta m11 */
			ofc << "# message block difference" << endl;
			for (unsigned i = 0; i < 16; ++i)
				for (unsigned b = 0; b < 32; ++b)
					if (naf(m_diff[i]).get(b) != 0)
						ofc << "diffm" << i << " = " << int(b+1)*naf(m_diff[i]).get(b) << endl;
			/* other settings */
			ofc << "# parameters" << endl;

			/* copy template file if any */
			ifstream ifc("md5diffpathconnect.cfg.template");
			if (!ifc)
				cerr << "Warning: could not read md5diffpathconnect.cfg.template" << endl;
			else
				ofc << ifc.rdbuf() << endl;
			cout << "Saved 'md5diffpathconnect.cfg'." << endl;
		}
	}

	/*** Write md5diffpath_helper.cfg ***/
	{
		ofstream ofc("md5diffpathhelper.cfg");
		if (!ofc)
			cerr << "Error: could not write md5diffpathhelper.cfg!" << endl;
		else {
			/* write delta m11 */
			ofc << "# message block difference" << endl;
			for (unsigned i = 0; i < 16; ++i)
				for (unsigned b = 0; b < 32; ++b)
					if (naf(m_diff[i]).get(b) != 0)
						ofc << "diffm" << i << " = " << int(b+1)*naf(m_diff[i]).get(b) << endl;

			cout << "Saved 'md5diffpathhelper.cfg'." << endl;
		}
	}

	return 0;
}











void writeupperpath(std::string dir, parameters_type& parameters, unsigned b, int m11sign)
{
	try {
		boost::filesystem::create_directory(dir);
	} catch( std::exception& e) {
		std::cerr << e.what() << endl;
		std::cerr << "failed to create directory: " << dir << endl;
		throw;
	} catch(...) {
		std::cerr << "failed to create directory: " << dir << endl;
		throw;
	}
	bitcondition plus = bc_plus, minus = bc_minus;
	if (m11sign == -1)
		swap(plus, minus);
	uint32 m_diff[16] = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
	m_diff[11] = m11sign * (1<<b);

	/*** Construct upper diff. path ***/
	differentialpath upperpath;
	upperpath[32].clear();
	upperpath.setbitcondition(60, (b+10)&31, bc_zero);
	upperpath.setbitcondition(61, (b+10)&31, bc_one);
	upperpath.setbitcondition(62, (b+10)&31, plus);
	upperpath.setbitcondition(63, (b+10)&31, plus);
	upperpath.setbitcondition(64, (b+10)&31, plus);

	/*** Show and store upper diff. path ***/
	try {
		vector<differentialpath> temp;
		temp.push_back(upperpath);
		save_gz(temp, dir +  "/upperpath", binary_archive);
		//cout << "Saved upper diff. path to '" << dir + "/upperpath.bin'." << endl;
	} catch (...) {
		cerr << "Error: could not write '" << dir << "/upperpath.bin'!" << endl;
	}

	/*** Write md5diffpath_backward.cfg ***/
	{
		ofstream ofb((dir+"/md5diffpathbackward.cfg").c_str());
		if (!ofb)
			cerr << "Error: could not write " << dir << "md5diffpathbackward.cfg!" << endl;
		else {
			/* write delta m11 */
			ofb << "# message block difference" << endl;
			ofb << "diffm11 = " << m11sign*int(b+1) << endl << endl;
			/* other settings */
			ofb << "# parameters" << endl;

			/* copy template file if any */
			ifstream ifb("md5diffpathbackward.cfg.template");
			if (!ifb)
				cerr << "Error: could not read md5diffpathbackward.cfg.template" << endl;
			else {
				ofb << ifb.rdbuf() << endl;
			}
			cout << "Saved '" << dir << "md5diffpathbackward.cfg'." << endl;
		}
	}
}

int upperpaths(parameters_type& parameters)
{
	for (unsigned b = 0; b < 32; ++b) 
	{
		string bitstr = boost::lexical_cast<string>(b);
		writeupperpath(workdir + "plus2power" + bitstr, parameters, b, 1);
		writeupperpath(workdir + "minus2power" + bitstr, parameters, b, -1);
	}
	return 0;
}

