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

#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1detail.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/progress_display.hpp>

#include "main.hpp"

unsigned int lowpathindex;

void random_permutation(vector<sha1differentialpath>& paths, path_container& container)
{
	// use a pseudo-random permutation fixed by inputs
	seed(paths.size()); addseed(container.t); addseed(container.modn);
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

unsigned lower_eqbits(const sha1differentialpath& _Left, const sha1differentialpath& _Right)
{
	unsigned t = _Left.tend() - 1;
	uint32 LdFt = 0 - _Left[t].getsdr().rotate_left(5).adddiff() - _Left[t-4].getsdr().rotate_left(30).adddiff();
	uint32 LdFtp1 = 0 - _Left[t-3].getsdr().rotate_left(30).adddiff();
	uint32 LdFtp2 = 0 - _Left[t-2].getsdr().rotate_left(30).adddiff();
	uint32 LdFtp3 = 0 - _Left[t-1].getsdr().rotate_left(30).adddiff();
	uint32 LdFtp4 = 0 - _Left[t].getsdr().rotate_left(30).adddiff();
	uint32 LdQtm1 = _Left[t-1].diff();
	uint32 LdQt = _Left[t].diff();
	uint32 RdFt = 0 - _Right[t].getsdr().rotate_left(5).adddiff() - _Right[t-4].getsdr().rotate_left(30).adddiff();
	uint32 RdFtp1 = 0 - _Right[t-3].getsdr().rotate_left(30).adddiff();
	uint32 RdFtp2 = 0 - _Right[t-2].getsdr().rotate_left(30).adddiff();
	uint32 RdFtp3 = 0 - _Right[t-1].getsdr().rotate_left(30).adddiff();
	uint32 RdFtp4 = 0 - _Right[t].getsdr().rotate_left(30).adddiff();
	uint32 RdQtm1 = _Right[t-1].diff();
	uint32 RdQt = _Right[t].diff();
	unsigned b = 0;
	while (b < 32+8) {
		if (b < 32) {
			if (_Left(t-3,(b+2)&31)!=_Right(t-3,(b+2)&31)) break;
			if (_Left(t-2,(b+2)&31)!=_Right(t-2,(b+2)&31)) break;
			if ((LdFt^RdFt)&(1<<b)) break;
			if ((LdQtm1^RdQtm1)&(1<<b)) break;
		}
		if (b >= 2 && b-2 < 32) {
			if ((LdFtp1^RdFtp1)&(1<<(b-2))) break;
			if ((LdQt^RdQt)&(1<<(b-2))) break;
		}
		if (b >= 4 && b-4 < 32 && ((LdFtp2^RdFtp2)&(1<<(b-4)))) break;
		if (b >= 6 && b-6 < 32 && ((LdFtp3^RdFtp3)&(1<<(b-6)))) break;
		if (b >= 8 && b-8 < 32 && ((LdFtp4^RdFtp4)&(1<<(b-8)))) break;
		++b;
	}
	return b;
}





progress_display* dostep_progress = 0;
unsigned dostep_index = 0;
struct dostep_thread {
	dostep_thread(vector<sha1differentialpath>& inlow, vector<sha1differentialpath>& inhigh, path_container& out, bool randomize)
		: pathsinlow(inlow), pathsinhigh(inhigh), container(out), random(randomize), worker(0)
	{}
	~dostep_thread()
	{ if (worker) delete worker; }
	vector<sha1differentialpath>& pathsinlow;
	vector<sha1differentialpath>& pathsinhigh;
	path_container& container;
	sha1_connect_thread* worker;
	sha1differentialpath prelowpath;
	vector<sha1differentialpath> lowpaths;
	bool random;
	void operator()() {
		try {
			worker = new sha1_connect_thread;
			while (true) {
				mut.lock();
				if (dostep_index >= pathsinlow.size())
				{
					mut.unlock();
					break;
				}
				// obtain a stripe of consecutive lowpaths
				static const unsigned stripe = (1<<5); // must be a power of 2
				unsigned count = stripe;
				if (dostep_index+count >= pathsinlow.size())
					count = pathsinlow.size() - dostep_index;
				dostep_index += count;
				(*dostep_progress) += count;
				lowpaths.clear();
				if (random) {
					unsigned kr;
					do {
						kr = (xrng128()%pathsinlow.size());
						kr &= ~(stripe-1); // set kr to beginning of stripe (set necessary low bits to zero)
					} while (pathsinlow[kr].path.size() == 0);
					for (unsigned i = 0; i < count; ++i)
					{
						lowpaths.push_back(pathsinlow[kr+i]);
						pathsinlow[kr+i].clear();
						worker->lowpathindex = kr;
					}
				} else {
					unsigned kr = dostep_index - count;
					for (unsigned i = 0; i < count; ++i)
					{
						lowpaths.push_back(pathsinlow[kr+i]);
						pathsinlow[kr+i].clear();
						worker->lowpathindex = kr;
					}
				}
				mut.unlock();

				for (unsigned i = 0; i < lowpaths.size(); ++i,++worker->lowpathindex)
				{
					if (prelowpath.path.size())
						worker->sha1_connect(lowpaths[i], pathsinhigh, container, lower_eqbits(lowpaths[i],prelowpath));
					else
						worker->sha1_connect(lowpaths[i], pathsinhigh, container, 0);
					prelowpath.swap(lowpaths[i]);
				}
			}
		} 
		catch (std::exception & e) { cerr << "Worker thread: caught exception:" << endl << e.what() << endl; } 
		catch (...) { cerr << "Worker thread: caught unknown exception:" << endl; }
	}       
};
void dostep_threaded(vector<sha1differentialpath>& inlow, vector<sha1differentialpath>& inhigh, path_container& out, bool randomize = true)
{
	mut.lock();
	dostep_index = 0;
	cout << "Starting threads..." << flush;
	boost::thread_group mythreads;
	for (unsigned i = 0; i < out.threads; ++i) {
		cout << i << flush;
		mythreads.create_thread(*new dostep_thread(inlow,inhigh,out,randomize));
		cout << i << flush;
	}
	cout << endl;
	std::string tstring = "t=" + boost::lexical_cast<std::string>(out.t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	dostep_progress = new progress_display(inlow.size(), true, cout, tstring, "      ", "      ");
	mut.unlock();
	mythreads.join_all();
	delete dostep_progress;
}





void load_upper_paths(vector<sha1differentialpath>& pathsinhigh, path_container& container)
{
        const unsigned t = container.t;
        const unsigned modn = container.modn;
        const unsigned modi = container.modi;
        const unsigned mode = container.splitmode; 
        // 0 = split upperpath sorted, 1 = split upperpath random, 2 = split lowerpath sorted, 3 = split lowerpath random
        try {
                cout << "Loading " << container.inputfilehigh << "..." << flush;
                if (modn > 1 && mode <= 1) {
                        vector<sha1differentialpath> tmp;
                        load_bz2(tmp, binary_archive, container.inputfilehigh);
                        if (mode == 1) {
                                random_permutation(tmp, container);
                                pathsinhigh.reserve((tmp.size()+modn)/modn); 
                                for (unsigned j = modi; j < tmp.size(); j += modn) {
                                        pathsinhigh.push_back(tmp[j]);
                                }
                                sort(pathsinhigh.begin(), pathsinhigh.end(), diffpathupper_less());
                        } else {
                                sort(tmp.begin(), tmp.end(), diffpathupper_less());
                                unsigned blockstart = (modi*tmp.size())/modn;
                                unsigned blockend = ((modi+1)*tmp.size())/modn;
                                pathsinhigh.reserve(blockend - blockstart);   
                                for (unsigned k = blockstart; k < blockend; ++k)
                                        pathsinhigh.push_back(tmp[k]);
                        }
                } else {
                        load_bz2(pathsinhigh, binary_archive, container.inputfilehigh);
                        sort(pathsinhigh.begin(), pathsinhigh.end(), diffpathupper_less());
                }
                cout << "done: " << pathsinhigh.size() << "." << endl;
        } catch(exception& e) {
                cout << "failed.\n(" << e.what() << ")" << endl;
        } catch(...) {
                cout << "failed." << endl;
        }
}

void load_lower_paths(vector<sha1differentialpath>& pathsinlow, path_container& container)
{
        const unsigned t = container.t;
        const unsigned modn = container.modn;
        const unsigned modi = container.modi;
        const unsigned mode = container.splitmode;
        try {
                cout << "Loading " << container.inputfilelow << "..." << flush;
                if (modn > 1 && mode >= 2) {
                        vector< sha1differentialpath > pathstmp2;
                        load_bz2(pathstmp2, binary_archive, container.inputfilelow);
                        if (mode == 3) {
                                random_permutation(pathstmp2, container);
                                pathsinlow.reserve((pathstmp2.size()/modn)+1);
                                for (unsigned j = modi; j < pathstmp2.size(); j += modn)
                                        pathsinlow.push_back(pathstmp2[j]);
                                sort(pathsinlow.begin(), pathsinlow.end(), diffpathlower_less());
                        } else {
                                sort(pathstmp2.begin(), pathstmp2.end(), diffpathlower_less());
                                unsigned blockstart = (modi*pathstmp2.size())/modn;
                                unsigned blockend = ((modi+1)*pathstmp2.size())/modn;
                                pathsinlow.reserve(blockend - blockstart);
                                for (unsigned k = blockstart; k < blockend; ++k)
                                        pathsinlow.push_back(pathstmp2[k]);
                        }                       
                } else {                        
                        load_bz2(pathsinlow, binary_archive, container.inputfilelow);
                        sort(pathsinlow.begin(), pathsinlow.end(), diffpathlower_less());
                }
                cout << "done: " << pathsinlow.size() << "." << endl;
        } catch(exception& e) {
                cout << "failed.\n(" << e.what() << ")" << endl;
        } catch(...) {
                cout << "failed." << endl;
        }
}



void dostep(path_container& container)
{
        const unsigned t = container.t;
        const unsigned modn = container.modn;
        const unsigned modi = container.modi;
        const unsigned mode = container.splitmode;
        // 0 = split upperpath sorted, 1 = split upperpath random, 2 = split lowerpath sorted, 3 = split lowerpath random
        vector< sha1differentialpath > pathsinhigh, pathsinlow, pathsout;
        if (mode <= 1) {
                load_upper_paths(pathsinhigh, container);
                load_lower_paths(pathsinlow, container); 
        } else {
                load_lower_paths(pathsinlow, container);
                load_upper_paths(pathsinhigh, container);
        }
        if (container.showinputpaths) {
                for (unsigned r = 0; r < pathsinlow.size(); ++r)
                        show_path(pathsinlow[r]);
                for (unsigned r = 0; r < pathsinhigh.size(); ++r)
                        show_path(pathsinhigh[r]);
        } else {
                if (pathsinlow.size())
                        show_path(pathsinlow[0]);
                if (pathsinhigh.size())
                        show_path(pathsinhigh[0]);
        }
        unsigned badcount = 0;
        for (unsigned i = 0; i < pathsinlow.size(); ++i)
                if (!test_path(pathsinlow[i]))
                        ++badcount;
        if (badcount) {
                cout << "Bad lowerpaths: " << badcount << " out of " << pathsinlow.size() << endl;
                for (unsigned i = 0; i < pathsinlow.size(); ++i)
                        if (!test_path(pathsinlow[i]) || pathsinlow[i].path.size()==0) {
                                show_path(pathsinlow[i]);
                                break;
                        }
        }
        badcount = 0;
        for (unsigned i = 0; i < pathsinhigh.size(); ++i)
                if (!test_path(pathsinhigh[i]))
                        ++badcount;
        if (badcount) {
                cout << "Bad upperpaths: " << badcount << " out of " << pathsinhigh.size() << endl;
                for (unsigned i = 0; i < pathsinhigh.size(); ++i)
                        if (!test_path(pathsinhigh[i]) || pathsinhigh[i].path.size()==0) {
                                show_path(pathsinhigh[i]);
                                break;
                        }
        }

	vector< sha1differentialpath > pathsredo, pathsredohigh, pathsredolow;
	if (container.inputfileredo.size()) {
		cout << "Loading redo inputs: " << flush;
		try {
			load_bz2(pathsredo, binary_archive, container.inputfileredo);
			cout << "done: " << pathsredo.size() << "." << endl;
		} catch (exception& e) {
			cout << "failed.\n(" << e.what() << ")" << endl;
			pathsredo.clear();
		} catch (...) {
			cout << "failed." << endl;
			pathsredo.clear();
		}
	}
	if (pathsredo.size() > 0) {
		cout << "Filtering..." << endl;
		progress_display pd(pathsinlow.size()+pathsinhigh.size());
		int lowtend = pathsinlow[0].tend()-1, hightbegin = pathsinhigh[0].tbegin();
		vector<uint32> redolowtend(pathsredo.size()), redohightbegin(pathsredo.size());
		if (lowtend+1 != hightbegin) {
			cout << "lowtend+1 != hightend: " << lowtend << " " << hightbegin << endl << endl;
		}
		for (unsigned i = 0; i < pathsredo.size(); ++i) {
			redolowtend[i] = pathsredo[i][lowtend].diff();
			redohightbegin[i] = pathsredo[i][hightbegin].diff();
		}
		for (unsigned k = 0; k < pathsinlow.size(); ++k,++pd) {
			bool redolow = false;
			uint32 tmp = pathsinlow[k][ pathsinlow[k].tend()-1 ].diff();
			for (unsigned j = 0; j < pathsredo.size(); ++j) {
				if (pathsinlow[k].tend()-1 != lowtend)
					throw std::runtime_error("tend()-1 != lowtend");
				if (tmp != redolowtend[j])
					continue;
//				if (pathsredo[j].tbegin() > pathsinlow[k].tbegin() || pathsredo[j].tend() < pathsinlow[k].tend())
//					continue;
				bool ok = true;
/*				for (int t = pathsinlow[k].tend()-2; t >= pathsinlow[k].tbegin(); --t)
					if (pathsinlow[k][t].diff() != pathsredo[j][t].diff()) {
						ok = false;
						break;
					}
*/				if (ok) {
					redolow = true;
					break;
				}
			}
			if (redolow) {
				pathsredolow.push_back(pathsinlow[k]);
			}
		}
		for (unsigned k = 0; k < pathsinhigh.size(); ++k,++pd) {
			bool redohigh = false;
			uint32 tmp = pathsinhigh[k][ pathsinhigh[k].tbegin() ].diff();
			for (unsigned j = 0; j < pathsredo.size(); ++j) {
				if (pathsinhigh[k].tbegin() != hightbegin) 
					throw std::runtime_error("pathsinhigh[k].tbegin() != hightbegin");
				if (tmp != redohightbegin[j])
					continue;
//				if (pathsredo[j].tbegin() > pathsinhigh[k].tbegin() || pathsredo[j].tend() < pathsinhigh[k].tend())
//					throw std::runtime_error("blaaat");
//					continue;
				bool ok = true;
/*				for (int t = pathsinhigh[k].tbegin()+1; t < pathsinhigh[k].tend(); ++t)
					if (pathsinhigh[k][t].diff() != pathsredo[j][t].diff()) {
						ok = false;
						break;
					}
*/				if (ok) {
					redohigh = true;
					break;
				}
			}
			if (redohigh) {
				pathsredohigh.push_back(pathsinhigh[k]);
			}
		}
		cout << "# Lowerpaths to redo: " << pathsredolow.size() << endl;
		cout << "# Upperpaths to redo: " << pathsredohigh.size() << endl;
		cout << "Processing selected lowerpaths..." << endl;
		dostep_threaded(pathsredolow, pathsinhigh, container, false);
		container.save_bestpaths();
		
		cout << "Processing selected upperpaths..." << endl;
		dostep_threaded(pathsinlow, pathsredohigh, container, false);
		container.save_bestpaths();
		return;
	}
	if (container.loworder == 0)
		dostep_threaded(pathsinlow, pathsinhigh, container, true);
	else
	{
		if (container.loworder == 2)
			sort(pathsinlow.begin(), pathsinlow.end(), diffpathlower_less_weight());
		dostep_threaded(pathsinlow, pathsinhigh, container, false);
	}
	container.save_bestpaths();
}
