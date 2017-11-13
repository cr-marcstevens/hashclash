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

#include "main.hpp"

using namespace hashclash;
using namespace std;

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







// map: (dQ1, dF1) => (count, Q1_set0, Q1_set1)
typedef map< pair<sdr,sdr>, pair<unsigned, pair<uint32, uint32> > > step1diffs_type;
void generate_step1diffs(step1diffs_type& step1diffs, const path_container_autobalance& container)
{
	const unsigned m = container.modn;
	const unsigned i = container.modi;
	const uint32* m_diff = container.m_diff;

	uint32 Q1[4], Q2[4];
	Q1[0] = container.IHV1[0]; Q1[1] = container.IHV1[3]; Q1[2] = container.IHV1[2]; Q1[3] = container.IHV1[1];
	Q2[0] = container.IHV2[0]; Q2[1] = container.IHV2[3]; Q2[2] = container.IHV2[2]; Q2[3] = container.IHV2[1];

	step1diffs_type::iterator step1diffs_it, step1diffs_it2;
	uint32 Ta = md5_ff(Q1[3], Q1[2], Q1[1]) + Q1[0] + md5_ac[0];
	uint32 Tb = md5_ff(Q2[3], Q2[2], Q2[1]) + Q2[0] + md5_ac[0] + m_diff[0];

	uint64 endcount = uint64(1)<<32;
	uint64 count = 0;
	unsigned k = 0;
	progress_display show_progress(endcount, true, cout, "t=0:  ", "      ", "      ");
	uint32 Q1a, Fa, Q1b, Fb;
	sdr dQ1, dF;
	pair<sdr,sdr> index;	
	for (uint64 count = 0; count < endcount; ++count,++show_progress,++Ta,++Tb)
	{
		++k;
		if (k == m) k = 0;
		if (k != i) continue;
		
		Q1a = rotate_left(Ta, 7) + Q1[Qoff+0];
		Fa = md5_ff(Q1a, Q1[Qoff+0], Q1[Qoff-1]);
		Q1b = rotate_left(Tb, 7) + Q2[Qoff+0];
		Fb = md5_ff(Q1b, Q2[Qoff+0], Q2[Qoff-1]);
		dQ1.set(Q1a, Q1b);
		dF.set(Fa, Fb);

		index.first = dQ1; index.second = dF;
		step1diffs_it = step1diffs.find(index);
		if (step1diffs_it == step1diffs.end()) {
			pair<unsigned, pair<uint32, uint32> > tmp;
			tmp.first = 1;
			tmp.second.first = Q1a;
			tmp.second.second = Q1a;
			step1diffs[index] = tmp;
		} else {
			++step1diffs_it->second.first;
			step1diffs_it->second.second.first |= Q1a;
			step1diffs_it->second.second.second &= Q1a;
		}
	}		
}

void combine_step1diffs(step1diffs_type& global_step1diffs, const step1diffs_type& step1diffs)
{
	step1diffs_type::const_iterator step1diffs_it;
	step1diffs_type::iterator step1diffs_it2;
	for (step1diffs_it = step1diffs.begin(); step1diffs_it != step1diffs.end(); ++step1diffs_it)
	{
		step1diffs_it2 = global_step1diffs.find(step1diffs_it->first);
		if (step1diffs_it2 == global_step1diffs.end())
			global_step1diffs.insert(*step1diffs_it);
		else {
			step1diffs_it2->second.first += step1diffs_it->second.first;
			step1diffs_it2->second.second.first |= step1diffs_it->second.second.first;
			step1diffs_it2->second.second.second &= step1diffs_it->second.second.second;
		}
	}
}

void dostep0(const path_container_autobalance& container)
{
	// map: (dQ1, dF1) => (count, Q1_set0, Q1_set1)
	map< pair<sdr,sdr>, pair<unsigned, pair<uint32, uint32> > >			step1diffs;
	map< pair<sdr,sdr>, pair<unsigned, pair<uint32, uint32> > >::iterator	step1diffs_it, step1diffs_it2;

	cout << "Searching first step differentials: " << endl;
	generate_step1diffs(step1diffs, container);
	hashclash::save_gz(step1diffs, pathsstring("paths0", container.modi, container.modn), binary_archive);
	cout << "Saved " << step1diffs.size() << " first step differentials." << endl;
}

void dostep1(const path_container_autobalance& container)
{
	const unsigned modn = container.modn;
	const unsigned modi = container.modi;

	map< pair<sdr,sdr>, pair<unsigned, pair<uint32, uint32> > >	step1diffs, step1diffs2;
	map< pair<sdr,sdr>, pair<unsigned, pair<uint32, uint32> > >::iterator	step1diffs_it, step1diffs_it2;
	
	if (modi != 0) return;
	for (unsigned j = 0; j < modn; ++j)
	{
		cout << "Loading " << pathsstring("paths0", j, modn) << "..." << flush;
		hashclash::load_gz(step1diffs2, pathsstring("paths0", j, modn), binary_archive);
		cout << "done." << endl;
		combine_step1diffs(step1diffs, step1diffs2);
	}

	unsigned count_max = 0;
	unsigned count_good = 0;

	step1diffs_it2 = step1diffs.begin();
	vector<unsigned> countvec;
	for (step1diffs_it = step1diffs.begin(); step1diffs_it != step1diffs.end(); ++step1diffs_it)
	{
		countvec.push_back(step1diffs_it->second.first);
		if (step1diffs_it->second.first > count_max)
			count_max = step1diffs_it->second.first;
	}
	std::sort(countvec.begin(), countvec.end());
	int countindex = int(countvec.size()) - int(65536);
	if (countindex < 0) countindex = 0;
	unsigned count_min = countvec[countindex];
		
	uint32 Q1[4], Q2[4];
	Q1[0] = container.IHV1[0]; Q1[1] = container.IHV1[3]; Q1[2] = container.IHV1[2]; Q1[3] = container.IHV1[1];
	Q2[0] = container.IHV2[0]; Q2[1] = container.IHV2[3]; Q2[2] = container.IHV2[2]; Q2[3] = container.IHV2[1];
	differentialpath tmp;
	for (int i = -3; i <= 0; ++i)
		tmp[i].set(Q2[3+i]-Q1[3+i], Q1[3+i], Q1[3+i]);

	vector< differentialpath > goodpaths;	
	for (step1diffs_it = step1diffs.begin(); step1diffs_it != step1diffs.end(); ++step1diffs_it)
	{
		if (step1diffs_it->second.first >= count_min)
		{
			++count_good;
			uint32 dQ1 = step1diffs_it->first.first.adddiff();
			uint32 dF1 = step1diffs_it->first.second.adddiff();
			uint32 dT1 = dF1 + (Q2[Qoff-2]-Q1[Qoff-2]) + container.m_diff[1];
			uint32 dR1 = best_rotated_difference(dT1, 12);
			uint32 dQ2 = dQ1 + dR1;

			tmp[1].set( dQ1, step1diffs_it->second.second.first, step1diffs_it->second.second.second );
			tmp[2] = naf(dQ2);
			goodpaths.push_back(tmp);
		}
	}
	cout << "Saving " << goodpaths.size() << " paths..." << flush;
	save_gz(goodpaths, pathsstring("paths1", modi, modn), binary_archive);
	cout << "done." << endl;
	if (goodpaths.size() > 0)
		show_path(goodpaths[0], container.m_diff);
}



progress_display* dostep_progress = 0;
unsigned dostep_index = 0;
struct dostep_thread {
	dostep_thread(vector<differentialpath>& in, path_container_autobalance& out)
		: pathsin(in), container(out)   
	{}
	vector<differentialpath>& pathsin;
	path_container_autobalance& container;
	md5_forward_thread worker;
	void operator()() {
		try {
			while (true) {
				mut.lock();
				unsigned i = dostep_index;
				unsigned iend = i + (unsigned(pathsin.size() - dostep_index)>>4);
				if (iend > i+4) iend = i+4;
				if (iend == i && i < pathsin.size())
					iend = i+1;
				if (iend != i)
					(*dostep_progress) += iend-i;
				dostep_index = iend;
				mut.unlock();
				if (i >= pathsin.size())
					break;
				for (; i < iend; ++i)
					worker.md5_forward_differential_step(pathsin[i], container);
			}
		} catch (std::exception & e) { cerr << "Worker thread: caught exception:" << endl << e.what() << endl; } catch (...) {}
	}       
};
void dostep_threaded(vector<differentialpath>& in, path_container_autobalance& out)
{
	dostep_index = 0;
	std::string tstring = "t=" + boost::lexical_cast<std::string>(out.t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	if (out.estimatefactor)
		dostep_progress = new progress_display(in.size(), true, cout, tstring, "      ", "e     ");
	else
		dostep_progress = new progress_display(in.size(), true, cout, tstring, "      ", "      ");
	boost::thread_group mythreads;
	for (unsigned i = 0; i < out.threads; ++i)
		mythreads.create_thread(dostep_thread(in,out));
	mythreads.join_all();
	if (dostep_progress->expected_count() != dostep_progress->count())
		*dostep_progress += dostep_progress->expected_count() - dostep_progress->count();
	delete dostep_progress;
}







vector<differentialpath> pathscache;
void dostep(path_container_autobalance& container, bool savetocache)
{
	const unsigned t = container.t;
	const unsigned modn = container.modn;
	const unsigned modi = container.modi;

	if (t == 0 && !container.normalt01) {
		dostep0(container);
		return;
	}
	if (t == 1 && !container.normalt01) {
		dostep1(container);
		return;
	}
	
	vector< differentialpath > pathsin, pathstmp, pathsout;
	if (pathscache.size() != 0) {
		pathsin.swap(pathscache);
		random_permutation(pathsin);
	} else if (container.inputfile.size() == 0) {
		for (unsigned k = 0; k < modn; ++k)
		{
			try {
				std::string filename = pathsstring("paths" + boost::lexical_cast<std::string>(t-1), k, modn);
				cout << "Loading " << filename << "..." << flush;
				load_gz(pathstmp, filename, binary_archive);
				random_permutation(pathstmp);
				for (unsigned j = modi; j < pathstmp.size(); j += modn)
					pathsin.push_back(pathstmp[j]);
				cout << "done: " << pathstmp.size() << " (work:" << pathsin.size() << ")." << endl;
			} catch(...) {
				cout << "failed." << endl;
			}
		}
	} else {
		bool failed = false;
		try {
			cout << "Loading " << container.inputfile << "..." << flush;
			load_gz(pathstmp, binary_archive, container.inputfile);
			random_permutation(pathstmp);
			for (unsigned j = modi; j < pathstmp.size(); j += modn)
				pathsin.push_back(pathstmp[j]);
			cout << "done: " << pathsin.size() << "." << endl;
		} catch(...) {
			failed = true;
			cout << "failed." << endl;
		}
		if (failed) {
			try {
				cout << "Loading (text) " << container.inputfile << "..." << flush;
				load_gz(pathstmp, text_archive, container.inputfile);
				random_permutation(pathstmp);
				for (unsigned j = modi; j < pathstmp.size(); j += modn)
					pathsin.push_back(pathstmp[j]);
				cout << "done: " << pathsin.size() << "." << endl;
			} catch(...) {
				cout << "failed." << endl;
			}
		}
	}
	if (container.showinputpaths)
		for (unsigned r = 0; r < pathsin.size(); ++r)
			show_path(pathsin[r], container.m_diff);

	std::string tstring = "t=" + boost::lexical_cast<std::string>(t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	
	if (container.estimatefactor != 0) {
		cout << "Estimating maxcond for upper bound " << unsigned(double(container.ubound)*container.estimatefactor)
			<< " (=" << container.ubound << " * " << container.estimatefactor << ")..." << endl;
		dostep_threaded(pathsin, container);
//		progress_display show_progress(pathsin.size(), true, cout, tstring, "      ", "e     ");
//		for (unsigned k = 0; k < pathsin.size(); ++k,++show_progress)
//			md5_forward_differential_step(pathsin[k], container);
		container.finish_estimate();
		cout << "Found maxcond = " << container.maxcond << endl;
	}

	dostep_threaded(pathsin, container);
//	progress_display show_progress(pathsin.size(), true, cout, tstring, "      ", "      ");
//	for (unsigned k = 0; k < pathsin.size(); ++k,++show_progress)
//		md5_forward_differential_step(pathsin[k], container);

	pathstmp.swap(pathsout);	
	container.export_results(pathsout);
	
	if (pathsout.size() > 0)
		show_path(pathsout[0], container.m_diff);
	else
		throw std::runtime_error("No valid differential paths found!");
	std::string filenameout = pathsstring("paths" + boost::lexical_cast<std::string>(t), modi, modn);
	cout << "Saving " << pathsout.size() << " paths..." << flush;
	if (savetocache)
		pathscache.swap(pathsout);
	else
		save_gz(pathsout, filenameout, binary_archive);
	cout << "done." << endl;
}
