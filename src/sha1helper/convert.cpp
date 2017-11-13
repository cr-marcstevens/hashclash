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

#include <cmath>
#include <algorithm>
#include <map>
#include <set>

#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/sha1detail.hpp>
#include <hashclash/sha1messagespace.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/progress_display.hpp>
#include <hashclash/rng.hpp>

#include <boost/lexical_cast.hpp>

#include "main.hpp"

void getline(std::istream& ifs)
{
	char c = 0;
	while (ifs && c != '\n' && c != '\r')
		ifs.get(c);
	while (ifs && (c == '\n' || c == '\r'))
		ifs.get(c);
	if (ifs)
		ifs.putback(c);
}

int pathfromtext(parameters_type& parameters)
{
	if (parameters.outfile1 == "") {
		cout << "No outputfile1 given!" << endl;
		return 2;
	}
	if (parameters.infile1 == "") {
		cout << "No inputfile1 given!" << endl;
		return 2;
	}
	parameters.show_mdiffs();
	sha1differentialpath path;
	vector<sha1differentialpath> vecpath;
	ifstream ifs(parameters.infile1.c_str());
	if (!ifs) {
		cerr << "Error: could not open " << parameters.infile1 << "!" << endl;
		return 1;
	}
	cout << "Parsing inputfile:";
	while (ifs) {
		cout << endl;
		char c;
		ifs >> c;
		cout << c << " ";
		if (!ifs) break;
		if (c != 'Q') {
			cout << "expected 'Q' here, going to next line";
			getline(ifs);
			continue;
		}
		int t;
		if (!(ifs >> t)) break;
		cout << t << " ";
		if (t < -4 || t > 80) {
			cout << "expected integer t with -4 <= t <= 80 here, going to next line";
			getline(ifs);
			continue;
		}
		while (ifs && c != '|') {
			ifs >> c;
			cout << c << " ";
		}
		if (c != '|') {
			cout << "expected '|' here, going to next line";
			getline(ifs);
			continue;
		}
		wordconditions wc;
		ifs >> wc.bytes[3] >> wc.bytes[2] >> wc.bytes[1] >> wc.bytes[0] >> c;
		cout << wc << " " << c << " ";
		if (!ifs) break;
		if (c != '|') {
			cout << "expected '|' here, going to next line";
			getline(ifs);
			continue;
		}
		path[t] = wc;
		path.getme(t).clear();
		if (t >= 0 && t < 80) {
			sdr met;
			ifs >> met;
			if (ifs)
				path.getme(t) = met;
		}
		getline(ifs);
	}
	cout << endl << "Parsed path:" << endl;
	show_path(path);
	vecpath.push_back(path);
	cout << "Saving " << parameters.outfile1 << "..." << flush;
	try {
		save_bz2(vecpath, binary_archive, parameters.outfile1);
		cout << "done." << endl;
	} catch (...) {
		cout << "failed." << endl;
		return 1;
	}
	return 0;
}

int split(parameters_type& parameters)
{
	vector<sha1differentialpath> vecpath, splitvec;
	cout << "Loading " << parameters.infile1 << "..." << flush;
	try {
		load_bz2(vecpath, binary_archive, parameters.infile1);
		cout << "done (loaded " << vecpath.size() << " paths)." << endl;
	} catch (...) {
		cout << "failed." << endl;
		return 1;
	}
	if (parameters.unique) {
		sort(vecpath.begin(), vecpath.end());
		vecpath.erase( unique(vecpath.begin(), vecpath.end()), vecpath.end() );
		cout << "Reduced to " << vecpath.size() << " unique paths." << endl;
	}
	double step = double(vecpath.size()) / double(parameters.split);
	uint32 start = 0;
	uint32 end = vecpath.size();
	for (unsigned i = 0; i < parameters.split; ++i)
	{
		uint32 iend = uint32(double(step*double(i+1)));
		if (i+1 == parameters.split || iend > end)
			iend = end;
		splitvec.clear();
		splitvec.reserve(iend-start);
		for (unsigned j = start; j < iend; ++j)
			splitvec.push_back(vecpath[j]);
		start = iend;
		string filename = parameters.infile1 + "_" 
					+ boost::lexical_cast<string>(i) + "of"
					+ boost::lexical_cast<string>(parameters.split) + ".bin";
		cout << "Saving " << splitvec.size() << " paths to " << filename << "..." << flush;
		try {
			save_bz2(splitvec, binary_archive, filename);
			cout << "done." << endl;
		} catch (...) {
			cout << "failed." << endl;
			return 1;
		}
	}
	return 0;
}

int join(parameters_type& parameters)
{
	if (parameters.outfile1 == "") {
		cout << "No outputfile1 given!" << endl;
		return 2;
	}
	vector<sha1differentialpath> vecpath, joinvec;
	for (unsigned i = 0; i < parameters.files.size(); ++i)
	{
		joinvec.clear();
		cout << "Loading " << parameters.files[i] << "..." << flush;
		try {
			load_bz2(joinvec, binary_archive, parameters.files[i]);
			cout << "done (loaded " << joinvec.size() << " paths)." << endl;
			for (unsigned j = 0; j < joinvec.size(); ++j)
				vecpath.push_back(joinvec[j]);
		} catch (...) {
			cout << "failed." << endl;
		}
	}
	if (parameters.unique) {
		for (unsigned i = 0; i < vecpath.size(); ++i) {
			for (int t = vecpath[i].tbegin(); t < vecpath[i].tend(); ++t) {
				if (!(t-4 >= vecpath[i].tbegin() && t+1 < vecpath[i].tend()))
					vecpath[i].getme(t) = sdr();
			}
		}
		sort(vecpath.begin(), vecpath.end());
		vecpath.erase( unique(vecpath.begin(), vecpath.end()), vecpath.end() );
		cout << "Reduced to " << vecpath.size() << " unique paths." << endl;
	}
	cout << "Saving " << vecpath.size() << " paths to " << parameters.outfile1 << "..." << flush;
	try {
		save_bz2(vecpath, binary_archive, parameters.outfile1);
		cout << "done." << endl;
	} catch (...) {
		cout << "failed." << endl;
		return 1;
	}
	return 0;
}

int convert(parameters_type& parameters)
{	
	sha1differentialpath path;	
	vector<sha1differentialpath> vecpath;
	set<sha1differentialpath> setpath;
	bool failed = true;
	if (failed) {
		cout << "Trying to read vector of differential paths in binary..." << flush;
		try {
			load_bz2(vecpath, binary_archive, parameters.infile1);
			if (vecpath.size() > 0) {
				failed = false;
				cout << "success." << endl;
				cout << "Trying to save vector of differential paths in text..." << flush;
				save_bz2(vecpath, text_archive, parameters.infile1 + ".txt.bz2");
				cout << "success." << endl;
			}
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed) {
		cout << "Trying to read vector of differential paths in text..." << flush;
		try {
			load_bz2(vecpath, text_archive, parameters.infile1);
			if (vecpath.size() > 0) {
				failed = false;
				cout << "success." << endl;
				cout << "Trying to save vector of differential paths in binary..." << flush;
				save_bz2(vecpath, binary_archive, parameters.infile1 + ".bin");
				cout << "success." << endl;
			}
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed) {
		cout << "Trying to read set of differential paths in binary..." << flush;
		try {
			load_bz2(setpath, binary_archive, parameters.infile1);
			if (setpath.size() > 0) {
				vecpath.assign(setpath.begin(),setpath.end());
				failed = false;
				cout << "success." << endl;
				cout << "Trying to save vector of differential paths in text and binary..." << flush;
				save_bz2(vecpath, text_archive, parameters.infile1 + ".txt.bz2");
				save_bz2(vecpath, binary_archive, parameters.infile1 + ".bin");
				cout << "success." << endl;
			}
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed) {
		cout << "Trying to read set of differential paths in text..." << flush;
		try {
			load_bz2(setpath, text_archive, parameters.infile1);
			if (setpath.size() > 0) {
				vecpath.assign(setpath.begin(),setpath.end());
				failed = false;
				cout << "success." << endl;
				cout << "Trying to save vector of differential paths in text and binary..." << flush;
				save_bz2(vecpath, text_archive, parameters.infile1 + ".txt.bz2");
				save_bz2(vecpath, binary_archive, parameters.infile1 + ".bin");
				cout << "success." << endl;
			}
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed) {
		cout << "Trying to read differential path in binary..." << flush;
		try {
			load_bz2(path, binary_archive, parameters.infile1);
			failed = false;
			cout << "success." << endl;
			cout << "Trying to save differential path in text..." << flush;
			save_bz2(path, text_archive, parameters.infile1 + ".txt.bz2");
			cout << "success." << endl;
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed) {
		cout << "Trying to read differential path in text..." << flush;
		try {
			load_bz2(path, text_archive, parameters.infile1);
			failed = false;
			cout << "success." << endl;
			cout << "Trying to save differential path in binary..." << flush;
			save_bz2(path, binary_archive, parameters.infile1 + ".bin");
			cout << "success." << endl;
		} catch (std::exception& e) {
			cout << "failed." << e.what() << endl;
		}
	}
	if (failed)
		return 1;
	return 0;
}


struct dqt_sort {
	uint32 data[5];

	bool operator< (const dqt_sort& r) const
	{
		if (data[4] < r.data[4]) return true;
		if (data[4] > r.data[4]) return false;
		if (data[3] < r.data[3]) return true;
		if (data[3] > r.data[3]) return false;
		if (data[2] < r.data[2]) return true;
		if (data[2] > r.data[2]) return false;
		if (data[1] < r.data[1]) return true;
		if (data[1] > r.data[1]) return false;
		if (data[0] < r.data[0]) return true;
		return false;
	}
	bool operator== (const dqt_sort& r) const
	{
		for (int i = 4; i >= 0; --i)
			if (data[i] != r.data[i]) return false;
		return true;
	}

	template<class Archive>
	void serialize(Archive& ar, const unsigned int file_version)
	{			
		ar & boost::serialization::make_nvp("data", data);
	}
};
void analyze_indepsection_prob(sha1messagespace& mespace, int tbegin, const uint32 dmmask[80], vector<dqt_sort>& target_dihvs);
int showmespaceconditions(parameters_type& parameters)
{
	sha1messagespace space;
	vector< vector<uint32> > bitrels;
	bool failed = true;
	cout << "Trying to read space in text..." << flush;
	try {
		load_bz2(space, text_archive, parameters.infile1);
		space.tobitrelations_80(bitrels);
		save(bitrels, text_archive, parameters.infile1 + ".bitrel.txt");
		failed = false;
		cout << "success." << endl;
	} catch (...) {
		cout << "failed." << endl;
	}
	if (failed) {
		cout << "Trying to read space in binary..." << flush;
		try {
			load_bz2(space, binary_archive, parameters.infile1);
			space.tobitrelations_80(bitrels);
			failed = false;
			cout << "success." << endl;
		} catch (...) {
			cout << "failed." << endl;
		}
	}
	if (failed) {
		cout << "Trying to read bitconditions in text..." << flush;
		try {
			load(bitrels, text_archive, parameters.infile1);
			space.frombitrelations_80(bitrels);
			failed = false;
			cout << "success." << endl;
			save_bz2(space, text_archive, parameters.infile1 + ".mespace.txt.bz2");
		} catch (std::exception& e) {
			cout << "failed:" << endl << e.what() << endl;
		} catch (...) {
			cout << "failed." << endl;
		}
	}
	if (failed) return 1;
        for (unsigned i = 0; i < bitrels.size(); ++i) {
                cout << " - ";
                bool firstone = true;
                for (unsigned t = 0; t < 80; ++t)
                        for (unsigned b = 0; b < 32; ++b)  
                                if (bitrels[i][t] & (1<<b)) {
                                        if (firstone)
                                                firstone = false;
                                        else
                                                cout << " + ";
                                        cout << "M[" << t << "," << b << "]";
                                }
                cout << " = " << (bitrels[i][80]&1) << endl;
        }
//        return 0;
	int tbegin = 67;
	vector<dqt_sort> tmp;
#if 1
	vector<dqt_sort> base(6);
	for (unsigned i = 0; i < 6; ++i) {
		base[i].data[0] = 1<<31;
		base[i].data[1] = 1<<1;
		base[i].data[2] = (i<2) ? (1<<31) : 0;
		base[i].data[3] = base[i].data[4] = 0;
		if (i==1 || i==4 || i==5) base[i].data[3] += 1<<4;
		if (i==0 || i==1) base[i].data[3] += 1<<6; else base[i].data[3] += 1<<7;
		if (i==0 || i==1) base[i].data[4] += (1<<11)+(1<<4)-(1<<2); else base[i].data[4] += (1<<12);
		if (i==1 || i==4 || i==5) base[i].data[4] += (1<<9);
		if (i==2 || i==4) base[i].data[4] += (1<<1)+(1<<3);
		if (i==3 || i==5) base[i].data[4] += (1<<4)-(1<<1);
	}
	for (unsigned t794 = 0; t794 < 2; ++t794)
	 for (unsigned t792 = 0; t792 < 2; ++t792)
  	  for (unsigned t787 = 0; t787 < 2; ++t787)
	   for (unsigned t783 = 0; t783 < 2; ++t783)
	    for (unsigned t763778 = 0; t763778 < 2; ++t763778)
	     for (unsigned t771772 = 0; t771772 < 2; ++t771772)
	     {
	   	for (unsigned i = 0; i < base.size(); ++i) {
	   		dqt_sort tmptmp = base[i];
	   		if (t794) tmptmp.data[4] -= (1<<5);
	   		if (t792) tmptmp.data[4] += (1<<3);
	   		if (t787) { tmptmp.data[3] -= (1<<8); tmptmp.data[4] -= (1<<13); }
	   		if (t783) { tmptmp.data[3] -= (1<<4); tmptmp.data[4] -= (1<<9); }
	   		if (t763778) tmptmp.data[1] -= (1<<2);
	   		if (t771772 && tmptmp.data[2]) { tmptmp.data[3] += (1<<7); tmptmp.data[4] += (1<<12); }
	   		tmp.push_back(tmptmp);
	   	}
	     }
	cout << "# dIHVs: " << tmp.size() << endl;
#endif
	analyze_indepsection_prob(space, tbegin, parameters.m_mask, tmp);
}

void analyze_indepsection_prob(sha1messagespace& mespace, int tbegin, const uint32 dmmask[80], vector<dqt_sort>& target_dihvs)
{
	const unsigned int tend = 80;
	cout << "Analyzing probability for t=[" << tbegin << "-" << tend << ")." << endl;

	vector< vector<uint32> > pathbitrelations16;
	vector< vector< vector<uint32> > > pathbitrelationsmatrix;
	mespace.tobitrelations_16(pathbitrelations16);
	cout << "Bitrelations: " << pathbitrelations16.size() << endl;
	pathbitrelationsmatrix.resize(16);
	for (unsigned i = 0; i < 16; ++i)
		pathbitrelationsmatrix[i].resize(32);
	for (unsigned i = 0; i < pathbitrelations16.size(); ++i) {
		int lastcol = -1;
		for (int col = 16*32-1; col >= 0; --col)
			if (pathbitrelations16[i][col>>5]&(1<<(col&31))) {
				lastcol = col;
				unsigned t = lastcol>>5;
				unsigned b = lastcol&31;
				pathbitrelationsmatrix[t][b] = pathbitrelations16[i];
				break;
			}
		if (lastcol == -1) throw;
	}

	const int offset = 4;
	uint32 Q[85];
	uint32 Q2[85];
	uint32 m[80];
	uint32 m2[80];

	bool forceddihvs = false;
	vector<dqt_sort> target_dihvs2;
	if (target_dihvs.size()) {
		forceddihvs = true;
		sort(target_dihvs.begin(), target_dihvs.end());
		cout << "Forcing target dIHVs (#=" << target_dihvs.size() << ")." << endl;
	}

	if (0) //target_dihvs.size() == 0)
	{
		cout << "Determining target dIHVs..." << endl;
		uint64 oksize = 1;
		while (true) {
			if (oksize >= (1<<28)) break;
			bool oksizeok = true;
			try {
				vector<dqt_sort> tmpdihvs(oksize*2);
			} catch (std::exception&) {oksizeok = false;} catch(...) {oksizeok = false;}
			if (!oksizeok)
				break;
			oksize *= 2;
		}
		uint64 bufsize = oksize;
		cout << "Maximum size temporary dIHV buffer: " << bufsize << endl;
		vector<dqt_sort> dihvs(bufsize);
		progress_display pd(dihvs.size());
		for (unsigned i = 0; i < dihvs.size(); ++i,++pd) {
			for (int t = 0; t < 16; ++t) {
				uint32 metmask = 0;
				uint32 metset1 = 0;
				for (unsigned b = 0; b < 32; ++b) {
					if (pathbitrelationsmatrix[t][b].size()) {
						metmask |= pathbitrelationsmatrix[t][b][t];
						uint32 v = pathbitrelationsmatrix[t][b][16]&1;
						for (int t1 = 0; t1 < t; ++t1)
							v ^= m[t1]&pathbitrelationsmatrix[t][b][t1];
						if (hw(v)&1)
							metset1 |= 1<<b;
					}
				}
				m[t] = (xrng128() & ~metmask) | metset1;
				m2[t] = m[t] ^ dmmask[t];
			}
			for (int t = 16; t < 80; ++t) {
				m[t]=rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);
				m2[t]=rotate_left(m2[t-3] ^ m2[t-8] ^ m2[t-14] ^ m2[t-16], 1);
			}
			for (int t = tbegin-4; t <= tbegin; ++t)
				Q2[offset+t] = Q[offset+t] = xrng128();
			for (int t = tbegin; t < 60; ++t) {
				sha1_step_round3(t, Q, m);
				sha1_step_round3(t, Q2, m2);
			}
			for (int t = (tbegin<60?60:tbegin); t < 80; ++t) {
				sha1_step_round4(t, Q, m);
				sha1_step_round4(t, Q2, m2);
			}
			dihvs[i].data[0] = rotate_left(Q2[offset+80-4],30)-rotate_left(Q[offset+80-4],30);
			dihvs[i].data[1] = rotate_left(Q2[offset+80-3],30)-rotate_left(Q[offset+80-3],30);
			dihvs[i].data[2] = rotate_left(Q2[offset+80-2],30)-rotate_left(Q[offset+80-2],30);
			dihvs[i].data[3] = Q2[offset+80-1]-Q[offset+80-1];
			dihvs[i].data[4] = Q2[offset+80]-Q[offset+80];
		}
		sort(dihvs.begin(), dihvs.end());
		unsigned i = 0;
		vector<unsigned> bestcnts;
		while (i < dihvs.size()) {
			unsigned j = i+1;
			while (j < dihvs.size() && dihvs[j] == dihvs[i]) ++j;
			bestcnts.push_back(j-i);
			i = j;
		}
		cout << "Prob t=[" << tbegin << "-" << 80 << "): " << flush;
		sort(bestcnts.begin(),bestcnts.end());
		unsigned bestcnt = bestcnts[bestcnts.size()-1];
		if (bestcnt == 1) {
			cerr << "Warning: bestcount = 1: no best dIHVs discernable!" << endl;
			exit(0);
		}

		i = 0;
		while (i < dihvs.size()) {
			unsigned j = i+1;
			while (j < dihvs.size() && dihvs[j] == dihvs[i]) ++j;
			if (forceddihvs) {
				if (binary_search(target_dihvs.begin(), target_dihvs.end(), dihvs[i])) {
					cout << "\tprob=" << log(double(j-i)/double(dihvs.size()))/log(2.0) << ":\t";
					for (unsigned k = 0; k < 5; ++k)
						cout << naf(dihvs[i].data[k]) << " ";
					cout << endl;
				} else if ((j-i)*2 >= bestcnt) {
					cout << "N.I.\tprob=" << log(double(j-i)/double(dihvs.size()))/log(2.0) << ":\t";
					for (unsigned k = 0; k < 5; ++k)
						cout << naf(dihvs[i].data[k]) << " ";
					cout << endl;
					target_dihvs2.push_back(dihvs[i]);
					sort(target_dihvs2.begin(), target_dihvs2.end());
				}
			} else {
				if ((j-i)*2 >= bestcnt) {
 					cout << "\tprob=" << log(double(j-i)/double(dihvs.size()))/log(2.0) << ":\t";
					for (unsigned k = 0; k < 5; ++k)
						cout << naf(dihvs[i].data[k]) << " ";
					cout << endl;
					target_dihvs.push_back(dihvs[i]);
				}
			}
			i = j;
		}
		{ vector<dqt_sort> tmptmp(1); dihvs.swap(tmptmp); } // free memory
		sort(target_dihvs.begin(), target_dihvs.end());
	}
	
	cout << "Checking cumulative probability target dIHVs: " << flush;
	uint64 cnt = 0, okcnt = 0, okcnt2 = 0;
	dqt_sort dihv;
	while (true) {
		for (int t = 0; t < 16; ++t) {
			uint32 metmask = 0;
			uint32 metset1 = 0;
			for (unsigned b = 0; b < 32; ++b) {
				if (pathbitrelationsmatrix[t][b].size()) {
					metmask |= pathbitrelationsmatrix[t][b][t];
					uint32 v = pathbitrelationsmatrix[t][b][16]&1;
					for (int t1 = 0; t1 < t; ++t1)
						v ^= m[t1]&pathbitrelationsmatrix[t][b][t1];
					if (hw(v)&1)
						metset1 |= 1<<b;
				}
			}
			m[t] = (xrng128() & ~metmask) | metset1;
			m2[t] = m[t] ^ dmmask[t];
		}
		for (int t = 16; t < 80; ++t) {
			m[t]=rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);
			m2[t]=rotate_left(m2[t-3] ^ m2[t-8] ^ m2[t-14] ^ m2[t-16], 1);
		}
		for (int t = tbegin-4; t <= tbegin; ++t)
			Q2[offset+t] = Q[offset+t] = xrng64();
		for (int t = tbegin; t < 60; ++t) {
			sha1_step_round3(t, Q, m);
			sha1_step_round3(t, Q2, m2);
		}
		for (int t = (tbegin<60?60:tbegin); t < 80; ++t) {
			sha1_step_round4(t, Q, m);
			sha1_step_round4(t, Q2, m2);
		}
		dihv.data[0] = rotate_left(Q2[offset+80-4],30)-rotate_left(Q[offset+80-4],30);
		dihv.data[1] = rotate_left(Q2[offset+80-3],30)-rotate_left(Q[offset+80-3],30);
		dihv.data[2] = rotate_left(Q2[offset+80-2],30)-rotate_left(Q[offset+80-2],30);
		dihv.data[3] = Q2[offset+80-1]-Q[offset+80-1];
		dihv.data[4] = Q2[offset+80]-Q[offset+80];
		++cnt;
		if (binary_search(target_dihvs.begin(), target_dihvs.end(), dihv)) {
			++okcnt;
			if (forceddihvs)
				++okcnt2;
		} else {
			if (forceddihvs && binary_search(target_dihvs2.begin(), target_dihvs2.end(), dihv))
				++okcnt2;
		}
		if (hw(cnt)+hw(cnt>>32)==1 && okcnt > 1) {
			cout << "[" << cnt << ":p=" << log(double(okcnt)/double(cnt))/log(2.0);
			if (forceddihvs)
				cout << "|p2=" << log(double(okcnt2)/double(cnt))/log(2.0);
			cout << "]" << flush;
		}
	}
	cout << endl;
}
