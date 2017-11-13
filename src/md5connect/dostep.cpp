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

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <utility>
#include <algorithm>
#include <string>
#include <iostream>
#include <time.h>

#include <boost/lexical_cast.hpp>

#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/progress_display.hpp>
#include <hashclash/timer.hpp>

#include "main.hpp"

vector<uint32> lowdQt, lowdQtm1, lowdQtm2, lowdQtm3;

void random_permutation(vector<differentialpath>& paths)
{
	// use a pseudo-random permutation fixed by the number of paths
	seed(paths.size());
	for (unsigned i = 0; i < paths.size(); ++i)
	{
		unsigned k = xrng64() % paths.size();
		paths[i].swap(paths[k]);
	}
	addseed(time(NULL));
}

inline std::string pathsstring(const std::string& basepath, unsigned modi, unsigned modn)
{
	return workdir +  "/"  + basepath 
		+ "_" + boost::lexical_cast<std::string>(modi) 
		+ "of" + boost::lexical_cast<std::string>(modn);
}


progress_display* dostep_progress = 0;
unsigned dostep_index = 0;
struct dostep_thread {
	dostep_thread(vector<differentialpath>& inlow, vector<differentialpath>& inhigh, path_container& out)
		: pathsinlow(inlow), pathsinhigh(inhigh), container(out)   
	{}
	vector<differentialpath>& pathsinlow;
	vector<differentialpath>& pathsinhigh;
	path_container& container;
	void operator()() {
		md5_connect_thread* worker = new md5_connect_thread;
		try {			
			while (true) {
				mut.lock();
				unsigned i = dostep_index;
				unsigned iend = i + (unsigned(pathsinhigh.size() - dostep_index)>>4);
				if (iend > i+4) iend = i+4;
				if (iend == i && i < pathsinhigh.size())
					iend = i+1;
				if (iend != i)
					(*dostep_progress) += iend-i;
				dostep_index = iend;
				mut.unlock();
				if (i >= pathsinhigh.size())
					break;
				for (; i < iend; ++i)
					worker->md5_connect(pathsinlow, pathsinhigh[i], container);
			}
		} catch (std::exception & e) { cerr << "Worker thread: caught exception:" << endl << e.what() << endl; } catch (...) {}
		delete worker;
	}       
};
void dostep_threaded(vector<differentialpath>& inlow, vector<differentialpath>& inhigh, path_container& out)
{
	dostep_index = 0;
	std::string tstring = "t=" + boost::lexical_cast<std::string>(out.t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	dostep_progress = new progress_display(inhigh.size(), true, cout, tstring, "      ", "      ");
	boost::thread_group mythreads;
	for (unsigned i = 0; i < out.threads; ++i)
		mythreads.create_thread(dostep_thread(inlow,inhigh,out));
	mythreads.join_all();
	if (dostep_progress->expected_count() != dostep_progress->count())
		*dostep_progress += dostep_progress->expected_count() - dostep_progress->count();
	delete dostep_progress;
}


void dostep(path_container& container)
{
	const unsigned t = container.t;
	const unsigned modn = container.modn;
	const unsigned modi = container.modi;

	vector< differentialpath > pathsinhigh, pathsout;
	if (modn > 1) {
		try {
			vector< differentialpath > pathstmp2;
			cout << "Loading " << container.inputfilehigh << "..." << flush;
			load_gz(pathstmp2, binary_archive, container.inputfilehigh);
			random_permutation(pathstmp2);
			for (unsigned j = modi; j < pathstmp2.size(); j += modn)
				pathsinhigh.push_back(pathstmp2[j]);
			sort(pathsinhigh.begin(), pathsinhigh.end(), diffpathupper_less());
			cout << "done: " << pathsinhigh.size() << "." << endl;
		} catch(...) {
			cout << "failed." << endl;
		}
	} else {
		try {
			cout << "Loading " << container.inputfilehigh << "..." << flush;
			load_gz(pathsinhigh, binary_archive, container.inputfilehigh);
			sort(pathsinhigh.begin(), pathsinhigh.end(), diffpathupper_less());
			cout << "done: " << pathsinhigh.size() << "." << endl;
		} catch(...) {
			cout << "failed." << endl;
		}
	}
	if (pathsinhigh.size() == 0)
		throw std::runtime_error("No upper differential paths loaded!");
	vector< differentialpath> pathsinlow;
	try {
		cout << "Loading " << container.inputfilelow << "..." << flush;
		load_gz(pathsinlow, binary_archive, container.inputfilelow);
		sort(pathsinlow.begin(), pathsinlow.end(), diffpathlower_less());
		cout << "done: " << pathsinlow.size() << "." << endl;
	} catch(...) {
		cout << "failed." << endl;
	}
	if (pathsinlow.size() == 0)
		throw std::runtime_error("No lower differential paths loaded!");	

	lowdQt.resize(pathsinlow.size());
	lowdQtm1.resize(pathsinlow.size());
	lowdQtm2.resize(pathsinlow.size());
	lowdQtm3.resize(pathsinlow.size());
	for (unsigned j = 0; j < pathsinlow.size(); ++j)
	{
		lowdQt[j] = pathsinlow[j][t].diff();
		lowdQtm1[j] = pathsinlow[j][t-1].diff();
		lowdQtm2[j] = pathsinlow[j][t-2].diff();
		lowdQtm3[j] = pathsinlow[j][t-3].diff();
	}

	if (container.showinputpaths)
		for (unsigned r = 0; r < pathsinlow.size(); ++r)
			show_path(pathsinlow[r], container.m_diff);

	std::string tstring = "t=" + boost::lexical_cast<std::string>(t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	
//	progress_display show_progress(pathsinhigh.size(), true, cout, tstring, "      ", "      ");
//	for (unsigned k = 0; k < pathsinhigh.size(); ++k,++show_progress)
//		md5_connect(pathsinlow, pathsinhigh[k], container);
	dostep_threaded(pathsinlow, pathsinhigh, container);
}

class check_path_collfind_helper {
public:
	differentialpath diffpath;
	uint32 m_diff[16];

	static const int offset = 3;
	uint32 m[16];
	uint32 m2[16];
	uint32 Q[68];
	uint32 Q2[68];

	uint32 dQ[68];
	uint32 dT[68];
	uint32 dR[68];

	uint32 Qvaluemask[68];
	uint32 Qvalue[68];
	uint32 Qprev[68];
	uint32 Qprev2[68];

	uint32 Q5m5tunnel, Q4m5tunnel, Q4m4tunnel;
	uint32 Q10m10tunnel, Q9m10tunnel, Q9m9tunnel;
	uint32 Q8Q12m15tunnel, Q14m6Q3m5tunnel, Q14Q3m14tunnel;

	inline bool checkrotation(uint32 Qt, uint32 Qtp1, uint32 dT, uint32 dR, unsigned rc)
	{
		uint32 R1 = Qtp1-Qt;
		uint32 R2 = R1 + dR;
		uint32 T1 = rotate_right(R1, rc);
		uint32 T2 = rotate_right(R2, rc);
		return (T2-T1 == dT);
	}

	bool check_step16()
	{
		hashclash::timer sw(true);
		for (unsigned i = 0; i < (1<<13); ++i) {
			Q[offset+1] = ((xrng128() & ~Qvaluemask[offset+1]) | Qvalue[offset+1]) ^ (Qprev[offset+1]&Q[offset+0]);
			if (!checkrotation(Q[offset+0], Q[offset+1], dT[offset+0], dR[offset+0], md5_rc[0]))
				continue;
			Q[offset+2] = ((xrng128() & ~Qvaluemask[offset+2]) | Qvalue[offset+2]) ^ (Qprev[offset+2]&Q[offset+1]);
			if (!checkrotation(Q[offset+1], Q[offset+2], dT[offset+1], dR[offset+1], md5_rc[1]))
				continue;
			// determine m0
			uint32 R0 = Q[offset+1] - Q[offset+0];
			uint32 T0 = rotate_right(R0, md5_rc[0]);
			uint32 F0 = md5_ff(Q[offset+0],Q[offset-1],Q[offset-2]);
			m[0] = T0 - F0 - Q[offset-3] - md5_ac[0];
			// determine m1
			uint32 R1 = Q[offset+2] - Q[offset+1];
			uint32 T1 = rotate_right(R1, md5_rc[1]);
			uint32 F1 = md5_ff(Q[offset+1],Q[offset+0],Q[offset-1]);
			m[1] = T1 - F1 - Q[offset-2] - md5_ac[1];

			Q[offset+13] = (xrng128() & (~Qvaluemask[offset+13] | Qprev[offset+13])) | Qvalue[offset+13];
			Q[offset+14] = ((xrng128() & ~Qvaluemask[offset+14]) | Qvalue[offset+14]) ^ (Qprev[offset+14]&Q[offset+13]);
			Q[offset+15] = ((xrng128() & ~Qvaluemask[offset+15]) | Qvalue[offset+15]) ^ (Qprev[offset+15]&Q[offset+14]);
			for (unsigned k = 0; k < 128; ++k)
			{
				Q[offset+16] = ((xrng128() & ~Qvaluemask[offset+16]) | Qvalue[offset+16]) ^ (Qprev[offset+16]&Q[offset+15]);
				if (!checkrotation(Q[offset+15], Q[offset+16], dT[offset+15], dR[offset+15], md5_rc[15]))
					continue;
				uint32 F16 = md5_gg(Q[offset+16],Q[offset+15],Q[offset+14]);
				uint32 T16 = F16 + m[1] + Q[offset+13] + md5_ac[16];
				uint32 R16 = rotate_left(T16, md5_rc[16]);
				Q[offset+17] = Q[offset+16] + R16;
				if (Qvalue[offset+17] != ((Q[offset+17]&Qvaluemask[offset+17])^(Qprev[offset+17]&Q[offset+16])))
					continue;
				uint32 T16b = T16 + dT[offset+16];
				uint32 R16b = rotate_left(T16b, md5_rc[16]);
				if (R16b - R16 != dR[offset+16]) 
					continue;

				return true;
			}
			if (sw.time() > 1.0) return false;
		}	
		return false;
	}

	void filltables()
	{
		for (int t = -3; t <= 64; ++t)
		{
			dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
			Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
			Qprev[offset+t] = Qprev2[offset+t] = 0;
		}
		Q5m5tunnel=0; Q4m5tunnel=0; Q4m4tunnel=0;
		Q10m10tunnel=0; Q9m10tunnel=0; Q9m9tunnel=0;
		Q8Q12m15tunnel=0; Q14m6Q3m5tunnel = 0; Q14Q3m14tunnel = 0;


		// determine tunnels and update differential path
		// Q14Q3m14 tunnel
		if (0)
		for (unsigned b = 0; b < 32; ++b)
		{
			bitcondition Q3b = diffpath(3,b), Q4b = diffpath(4,b), Q5b = diffpath(5,b);
			bitcondition Q14b = diffpath(14,b), Q15b = diffpath(15,b), Q16b = diffpath(16,b);
			if (Q3b == bc_constant && Q14b == bc_constant
				&& (Q4b == bc_constant || Q4b == bc_zero || Q4b == bc_plus)
				&& (Q5b == bc_constant || Q5b == bc_one || Q5b == bc_minus || Q5b == bc_prevn)
				&& (Q15b == bc_constant || Q15b == bc_zero || Q15b == bc_plus)
				&& (Q16b == bc_constant || Q16b == bc_constant || Q16b == bc_plus || Q16b == bc_prev))
			{
				Q14Q3m14tunnel |= 1<<b;
				diffpath.setbitcondition(3,b,bc_zero);
				diffpath.setbitcondition(14,b,bc_zero);
				if (Q4b == bc_constant)
					diffpath.setbitcondition(4,b,bc_zero);
				if (Q5b == bc_constant || Q5b == bc_prevn)
					diffpath.setbitcondition(5,b,bc_one);
				if (Q15b == bc_constant)
					diffpath.setbitcondition(15,b,bc_zero);
				if (Q16b == bc_constant || Q16b == bc_prev)
					diffpath.setbitcondition(16,b,bc_zero);
			}
		}

		// Q8Q12m15 tunnel
		for (unsigned b = 0; b < 32; ++b)
		{
			bitcondition Q8b = diffpath(8,b), Q9b = diffpath(9,b), Q10b = diffpath(10,b);
			bitcondition Q12b = diffpath(12,(b+22)&31), Q13b = diffpath(13,(b+22)&31);
			if (Q8b == bc_constant && Q12b == bc_constant
					&& Q9b != bc_prev && Q9b != bc_prevn
					&& Q13b != bc_prev && Q13b != bc_prevn
					&& (Q10b == bc_constant || Q10b == bc_minus || Q10b == bc_one)
				)
			{
				Q8Q12m15tunnel |= 1<<b;
				diffpath.setbitcondition(8,b,bc_zero);
				diffpath.setbitcondition(12,(b+22)&31,bc_zero);
				if (Q10b == bc_constant)
					diffpath.setbitcondition(10,b,bc_one);
			}
		}
		// Q14m6Q3m5 tunnel
		for (unsigned b = 0; b < 32; ++b)
		{
			bitcondition Q14b = diffpath(14,b), Q3b = diffpath(3,b)
						, Q15b = diffpath(15,b), Q16b = diffpath(16,b)
						, Q4b = diffpath(4,b);
			if (Q14b == bc_constant && Q3b == bc_constant
					&& Q15b != bc_prev && Q15b != bc_prevn
					&& Q4b != bc_prev && Q4b != bc_prevn)
			{
				bool istunnel = false;
				if (Q16b == bc_prev)
					istunnel = true;
				else if (Q16b == bc_constant) {
					istunnel = true;
					if (Q15b == bc_constant)
						diffpath.setbitcondition(16,b,bc_prev);
					else if (Q15b == bc_zero || Q15b == bc_plus)
						diffpath.setbitcondition(16,b,bc_zero);
					else if (Q15b == bc_one || Q15b == bc_minus)
						diffpath.setbitcondition(16,b,bc_one);
					else
						istunnel = false;
				} else if (Q16b == bc_zero || Q16b == bc_plus) {
					if (Q15b == bc_constant)
						diffpath.setbitcondition(15,b,bc_zero);
					if (Q15b == bc_zero || Q15b == bc_plus)
						istunnel = true;
				} else if (Q16b == bc_one || Q16b == bc_minus) {
					if (Q15b == bc_constant)
						diffpath.setbitcondition(15,b,bc_one);
					if (Q15b == bc_one || Q15b == bc_minus)
						istunnel = true;
				}
				if (istunnel) {
					diffpath.setbitcondition(14,b,bc_zero);
					diffpath.setbitcondition(3,b,bc_zero);
					Q14m6Q3m5tunnel |= 1<<b;
				}
			}
		}
		// Q4m4, Q4m5, Q5m5 tunnels
		for (unsigned b = 0; b < 32; ++b)
		{
			bitcondition Q4b = diffpath(4,b), Q5b = diffpath(5,b), Q6b = diffpath(6,b);
			if ((Q4b == bc_constant) 
				&& (Q5b == bc_constant || Q5b == bc_zero || Q5b == bc_plus)
				&& (Q6b == bc_constant || Q6b == bc_one || Q6b == bc_minus || Q6b == bc_prevn))
			{
				if (Q5b == bc_constant) 
					diffpath.setbitcondition(5,b,bc_zero);
				if (Q6b == bc_constant || Q6b == bc_prevn)
					diffpath.setbitcondition(6,b,bc_one);
				diffpath.setbitcondition(4,b,bc_zero);
				Q4m4tunnel |= 1<<b;
			} 
			else if ((Q4b == bc_constant)
				&& (Q5b == bc_constant || Q5b == bc_one || Q5b == bc_minus)
				&& (Q6b == bc_constant || Q6b == bc_one || Q6b == bc_minus || Q6b == bc_prev))
			{
				if (Q5b == bc_constant) 
					diffpath.setbitcondition(5,b,bc_one);
				if (Q6b == bc_constant || Q6b == bc_prev)
					diffpath.setbitcondition(6,b,bc_one);
				diffpath.setbitcondition(4,b,bc_zero);
				Q4m5tunnel |= 1<<b;
			} 
			else if ((Q5b == bc_constant)
				&& (Q6b == bc_constant || Q6b == bc_zero || Q6b == bc_plus))
			{
				if (Q6b == bc_constant)
					diffpath.setbitcondition(6,b,bc_zero);
				diffpath.setbitcondition(5,b,bc_zero);
				Q5m5tunnel |= 1<<b;			
			}
		}
		// Q9m9, Q9m10, Q10m10 tunnels
		for (unsigned b = 0; b < 32; ++b)
		{
			bitcondition Q9b = diffpath(9,b), Q10b = diffpath(10,b), Q11b = diffpath(11,b);
			if ((Q9b == bc_constant) 
				&& (Q10b == bc_constant || Q10b == bc_zero || Q10b == bc_plus)
				&& (Q11b == bc_constant || Q11b == bc_one || Q11b == bc_minus || Q11b == bc_prevn))
			{
				// if possible we maintain a dynamic tunnel on Q9m9 Q9m10
				if (Q10b == bc_constant && Q11b == bc_prevn) 
					diffpath.setbitcondition(10,b,bc_zero);

				if (Q11b == bc_constant || Q11b == bc_prevn)
					diffpath.setbitcondition(11,b,bc_one);
				diffpath.setbitcondition(9,b,bc_zero);
				Q9m9tunnel |= 1<<b;			
			} 
			else if ((Q9b == bc_constant)
				&& (Q10b == bc_constant || Q10b == bc_one || Q10b == bc_minus)
				&& (Q11b == bc_constant || Q11b == bc_one || Q11b == bc_minus || Q11b == bc_prev))
			{
				// if possible we maintain a dynamic tunnel on Q9m9 Q9m10
				if (Q10b == bc_constant && Q11b == bc_prev) 
					diffpath.setbitcondition(10,b,bc_one);

				if (Q11b == bc_constant || Q11b == bc_prev)
					diffpath.setbitcondition(11,b,bc_one);
				diffpath.setbitcondition(9,b,bc_zero);
				Q9m10tunnel |= 1<<b;			
			} 
			else if ((Q10b == bc_constant)
				&& (Q11b == bc_constant || Q11b == bc_zero || Q11b == bc_plus))
			{
				if (Q11b == bc_constant)
					diffpath.setbitcondition(11,b,bc_zero);
				diffpath.setbitcondition(10,b,bc_zero);
				Q10m10tunnel |= 1<<b;			
			}
		}

		// build tables
		for (int t = diffpath.tbegin(); t < diffpath.tend(); ++t)
		{
			dQ[offset+t] = diffpath[t].diff();
			Qprev[offset+t] = diffpath[t].prev() | diffpath[t].prevn();
			Qprev2[offset+t] = diffpath[t].prev2() | diffpath[t].prev2n();
			Qvaluemask[offset+t] = (~diffpath[t].set0()) | diffpath[t].set1()
								| diffpath[t].prev() | diffpath[t].prevn()
								| diffpath[t].prev2() | diffpath[t].prev2n();
			Qvalue[offset+t] = diffpath[t].set1() | diffpath[t].prevn() | diffpath[t].prev2n();
		}
		for (int t = 0; t < 64; ++t)
		{
			booleanfunction* F = 0;
			if (t < 16) F = & MD5_F_data;
			else if (t < 32) F = & MD5_G_data;
			else if (t < 48) F = & MD5_H_data;
			else F = & MD5_I_data;
			uint32 dF = 0;
			for (unsigned b = 0; b < 32; ++b)
			{
				bf_outcome outcome = F->outcome(diffpath(t,b), diffpath(t-1,b), diffpath(t-2,b));
				if (outcome.size()) {
					if (outcome[0] == bc_plus) 			dF += 1<<b;
					else if (outcome[0] == bc_minus)	dF -= 1<<b;
				} else throw std::runtime_error("ambiguous path!!");
			}
			dT[offset+t] = dQ[offset+t-3] + dF + m_diff[md5_wt[t]];
			dR[offset+t] = dQ[offset+t+1] - dQ[offset+t];
		}

		Q[0] = Qvalue[0]; Q2[0] = Q[0] + dQ[0];
		Q[1] = Qvalue[1]; Q2[1] = Q[1] + dQ[1];
		Q[2] = Qvalue[2]; Q2[2] = Q[2] + dQ[2];
		Q[3] = Qvalue[3]; Q2[3] = Q[3] + dQ[3];
	}
};

bool check_path_collfind(const differentialpath& diffpath, const uint32 mdiff[16])
{
	check_path_collfind_helper helper;
	helper.diffpath = diffpath;
	for (unsigned k = 0; k < 16; ++k)
		helper.m_diff[k] = mdiff[k];
	helper.filltables();
	return helper.check_step16();
}