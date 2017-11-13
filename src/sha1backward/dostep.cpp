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



void random_permutation(vector<sha1differentialpath>& paths, path_container_autobalance& container)
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

bool beginpathless(const sha1differentialpath& l, const sha1differentialpath& r)
{
	for (int i = 0; i < 3; ++i)
		if (l[l.tbegin()+i] != r[r.tbegin()+i])
			return l[l.tbegin()+i] < r[r.tbegin()+i];
	for (int i = 3; i < 5; ++i)
		if (l[l.tbegin()+i].getsdr() != r[r.tbegin()+i].getsdr())
			return l[l.tbegin()+i].getsdr() < r[r.tbegin()+i].getsdr();
	return false;
}

bool beginpatheq(const sha1differentialpath& l, const sha1differentialpath& r)
{
	for (int i = 0; i < 3; ++i)
		if (l[l.tbegin()+i] != r[r.tbegin()+i])
			return false;
	for (int i = 3; i < 5; ++i)
		if (l[l.tbegin()+i].getsdr() != r[r.tbegin()+i].getsdr())
			return false;
	return true;
}

void best_filter(vector<sha1differentialpath>& pathsout)
{
	vector<sha1differentialpath> ret;
	std::cout << "Filtering best paths (insize=" << pathsout.size() << " sorting..." << std::flush;
	std::sort(pathsout.begin(), pathsout.end(), beginpathless);
	std::cout << "filtering..." << std::flush;
	std::size_t i = 0;
	while (i < pathsout.size())
	{
		unsigned nrcondi = pathsout[i].nrcond();
		std::size_t e = i + 1;
		while (e < pathsout.size() && beginpatheq(pathsout[i],pathsout[e]))
		{
			unsigned nrconde = pathsout[e].nrcond();
			if (nrconde < nrcondi)
			{
				std::swap(pathsout[i], pathsout[e]);
				nrcondi = nrconde;
			}
			++e;
		}
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
		ret.push_back(pathsout[i]);
#else
		ret.push_back(std::move(pathsout[i]));
#endif
		i = e;
	}
	pathsout = std::move(ret);
	std::cout << " outsize=" << pathsout.size() << ")" << std::endl;
}


vector<sha1differentialpath> pathscache;
void dostep(path_container_autobalance& container, bool savetocache)
{
	const unsigned t = container.t;
	const unsigned modn = container.modn;
	const unsigned modi = container.modi;

	vector< sha1differentialpath > pathsin, pathstmp, pathsout;
	if (pathscache.size() != 0) {
		pathsin.swap(pathscache);
		random_permutation(pathsin, container);
	} else if (container.newinputpath) {
		sha1differentialpath path;
		path.offset = - int(t) + 4;
		path.path.resize(6);
		path.me.resize(6);
		pathsin.push_back( path );
		cout << "Generated 1 new path." << endl;
	} else if (container.inputfile.size() == 0) {
		for (unsigned k = 0; k < modn; ++k)
		{
			try {
				std::string filename = pathsstring("paths" + boost::lexical_cast<std::string>(t+1), k, modn);
				cout << "Loading " << filename << "..." << flush;
				load_bz2(pathstmp, filename, binary_archive);
				random_permutation(pathstmp, container);
				for (unsigned j = modi; j < pathstmp.size(); j += modn)
					pathsin.push_back(pathstmp[j]);
				cout << "done: " << pathstmp.size() << " (work:" << pathsin.size() << ")." << endl;
			} catch(...) {
				cout << "failed." << endl;
			}
		}
	} else {
		try {
			cout << "Loading " << container.inputfile << "..." << flush;
			load_bz2(pathstmp, binary_archive, container.inputfile);
			random_permutation(pathstmp, container);
			for (unsigned j = modi; j < pathstmp.size(); j += modn)
				pathsin.push_back(pathstmp[j]);
			cout << "done: " << pathsin.size() << "." << endl;
		} catch(...) {
			cout << "failed." << endl;
		}
	}
	if (container.showinputpaths)
		for (unsigned r = 0; r < pathsin.size(); ++r)
			show_path(pathsin[r]);

	std::string tstring = "t=" + boost::lexical_cast<std::string>(t) + ": ";
	if (tstring.size() == 5) tstring += " ";
	
	if (container.estimatefactor != 0) {
		cout << "Estimating maxcond for upper bound " << unsigned(double(container.ubound)*container.estimatefactor)
			<< " (=" << container.ubound << " * " << container.estimatefactor << ")..." << endl;
		progress_display show_progress(pathsin.size(), true, cout, tstring, "      ", "e     ");
		for (unsigned k = 0; k < pathsin.size(); ++k,++show_progress) {
			if (container.t < 79 && container.expandprevmessagediff) {
				sdr mmask = pathsin[k].getme(t+1);
				uint32 deltam = mmask.adddiff();
				uint32 dqtm3 = pathsin[k][t-3].getsdr().rotate_left(30).adddiff();

				uint32 addmask = (~mmask.mask)+1; 
				uint32 andmask = mmask.mask & 0x7FFFFFF;
				mmask.sign = 0;
				do {
					mmask.sign += addmask; mmask.sign &= andmask;
					pathsin[k].getme(t+1) = mmask;
					pathsin[k][t-3] = naf(dqtm3 + deltam - mmask.adddiff()).rotate_right(30);
					sha1_backward_differential_step(pathsin[k], container);
				} while (mmask.sign != 0);
			} else
				sha1_backward_differential_step(pathsin[k], container);
		}			
		container.finish_estimate();
		cout << "Found maxcond = " << container.maxcond << endl;
	}

	progress_display show_progress(pathsin.size(), true, cout, tstring, "      ", "      ");
	for (unsigned k = 0; k < pathsin.size(); ++k,++show_progress) {
		if (container.t < 79 && container.expandprevmessagediff) {
			sdr mmask = pathsin[k].getme(t+1);
			uint32 deltam = mmask.adddiff();
			uint32 dqtm3 = pathsin[k][t-3].getsdr().rotate_left(30).adddiff();

			uint32 addmask = (~mmask.mask)+1; 
			uint32 andmask = mmask.mask & 0x7FFFFFF;
			mmask.sign = 0;
			do {
				mmask.sign += addmask; mmask.sign &= andmask;
				pathsin[k].getme(t+1) = mmask;
				pathsin[k][t-3] = naf(dqtm3 + deltam - mmask.adddiff()).rotate_right(30);
				sha1_backward_differential_step(pathsin[k], container);
			} while (mmask.sign != 0);
		} else
			sha1_backward_differential_step(pathsin[k], container);
	}

	pathstmp.swap(pathsout);	
	container.export_results(pathsout);
	
	best_filter(pathsout);
	
	if (pathsout.size() > 0)
		show_path(pathsout[0]);
	std::string filenameout = pathsstring("paths" + boost::lexical_cast<std::string>(t), modi, modn);
	cout << "Saving " << pathsout.size() << " paths..." << flush;
	if (savetocache)
		pathsout.swap(pathscache);
	else
		save_bz2(pathsout, filenameout, binary_archive);
	cout << "done." << endl;
}
