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

	inline uint32 sha1_step_round1_bw(int t, uint32 Q[])
	{
		const int offset = 4;
		uint32 Ft = sha1_f1(Q[offset+t-1], rotate_left(Q[offset+t-2],30), rotate_left(Q[offset+t-3],30));
		return Q[offset+t+1] - ( Ft + sha1_ac[0] + rotate_left(Q[offset+t],5) + rotate_left(Q[offset+t-4],30) );
	}


struct filter_t
{
	const int Qoffset = 4;

	uint32 m[80];
	uint32 Q[85];
	
	uint32 Qcondmask[85];
	uint32 Qset1mask[85];
	uint32 Qprevmask[85];
	
	uint32 mcondmask[80];
	uint32 mset1mask[80];

	sha1differentialpath& path;
	parameters_type& parameters;
	filter_t(sha1differentialpath& _path, parameters_type& _parameters)
		: path(_path), parameters(_parameters)
	{
		for (int i = 0; i < 80; ++i) m[i] = 0;
		for (int i = 0; i < 80; ++i) mcondmask[i] = 0;
		for (int i = 0; i < 80; ++i) mset1mask[i] = 0;
		for (int i = 0; i < 85; ++i) Q[i] = 0;
		for (int i = 0; i < 85; ++i) Qcondmask[i] = 0;
		for (int i = 0; i < 85; ++i) Qset1mask[i] = 0;
		for (int i = 0; i < 85; ++i) Qprevmask[i] = 0;

//		show_path(path);
		for (int t = -4; t <= 20 && t < path.tend(); ++t)
		{
//			cout << "Q" << t << "\t";
			for (unsigned b = 0; b < 32; ++b)
			{
				bitcondition bc = path(t,b);

				if (bc != bc_constant)
				switch (bc)
				{
					case bc_one: case bc_minus:
						Qset1mask[Qoffset+t] |= 1<<b;
					case bc_zero: case bc_plus:
						Qcondmask[Qoffset+t] |= 1<<b;
						break;
					case bc_prevn:
						Qset1mask[Qoffset+t] |= 1<<b;
					case bc_prev:
						Qprevmask[Qoffset+t] |= 1<<b;
						Qcondmask[Qoffset+t] |= 1<<b;
						break;
					default:
						cout << "========= " << bc << " =======" << endl;
						throw;
				}
//				cout << bc;
								
			}
			
			if (t >= 0 && t < 20)
			{
				sdr mt = path.getme(t);
				mcondmask[t] = mt.mask;
				mset1mask[t] = mt.set1conditions();				
				
//				cout << mt;
			}
			
//			cout << endl;
		}
	}
	
	typedef vector<uint32> sol_t;
	
	void gensols(vector<sol_t>& sols, int solt)
	{
		sols.clear();
		if (solt == -5)
		{
			sols.resize(1);
			return;
		}

		vector<sol_t> soltm1;
		gensols(soltm1, solt-1);

		sol_t sol;		
		for (auto& s : soltm1)
		{
			sol = s;
			sol.push_back( 0 );

			uint32 set1 = Qset1mask[Qoffset+solt];
			if (solt > -4)
				set1 ^= (Qprevmask[Qoffset+solt] & sol[Qoffset+solt-1]);
				
			const uint32 mask = ~Qcondmask[Qoffset+solt];
			uint32 v = 0;
			do
			{
				--v;
				v &= mask;

				sol[Qoffset+solt] = v ^ set1;
				
				if (solt > 0)
				{
					uint32 m = sha1_step_round1_bw(solt-1, &sol[0]);
					if (((m^mset1mask[solt-1]) & mcondmask[solt-1]) != 0)
						continue;
				}
				
				sols.push_back(sol);			
				
			} while (v != 0);
			
		}		
		
//		cout << "# sols @ t=" << solt << ": " << sols.size() << endl;
	}
	
};


bool filterpath(sha1differentialpath& path, parameters_type& parameters)
{
	cleanup_path(path);

	filter_t filter(path, parameters);
	
	vector< filter_t::sol_t > sols;
	filter.gensols(sols, parameters.filtert);
	if (0 && !sols.empty())
	{
		show_path(path);
		cout << sols.size() << endl;
	}
	return ! sols.empty();
}

int filterfeasible(parameters_type& parameters)
{
	vector<sha1differentialpath> vecpath;
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
	
	uint64 cnt = 0, cnt2 = 0;
	for (std::size_t i = 0; i < vecpath.size();)
	{
		if (hw(++cnt) == 1)
			cout << "cnt = " << cnt << endl;
		if (filterpath(vecpath[i], parameters))
		{
			++i;
			if (hw(++cnt2) == 1)
			cout << "cnt2 = " << cnt2 << endl;
		
		}
		else
		{
			if (i != vecpath.size()-1)
				std::swap(vecpath[i], vecpath.back());
			vecpath.pop_back();
		}
	}
	cout << "Differential paths remaining: " << vecpath.size() << endl;
	if (!parameters.outfile1.empty())
	{
		cout << "Saving " << parameters.outfile1 << "..." << flush;
		save_bz2(vecpath, binary_archive, parameters.outfile1);
		cout << "done." << endl;
	}
	return 0;
}






int selectpath(parameters_type& parameters)
{
	vector<sha1differentialpath> vecpath;
	cout << "Loading " << parameters.infile1 << "..." << flush;
	try {
		load_bz2(vecpath, binary_archive, parameters.infile1);
		cout << "done (loaded " << vecpath.size() << " paths)." << endl;
	} catch (...) {
		cout << "failed." << endl;
		return 1;
	}
	if (parameters.seli < vecpath.size())
	{
		sha1differentialpath path = vecpath[parameters.seli];
		if (parameters.invert)
		{
			for (int t = path.tbegin(); t < path.tend(); ++t)
			{
				if (t-4 >= path.tbegin() && t+1 < path.tend())
					path.getme(t) = -path.getme(t);
				for (unsigned b = 0; b < 32; ++b)
				{
					bitcondition bc = path(t,b);
					if (bc == bc_plus)
						path.setbitcondition(t,b,bc_minus);
					else if (bc == bc_minus)
						path.setbitcondition(t,b,bc_plus);
				}
			}
		}
		show_path(path);
		cout << "Saving " << parameters.outfile1 << "..." << flush;
		save_bz2(path, binary_archive, parameters.outfile1);
		cout << "done." << endl;
		return 0;
	}
	cout << "Index is out of bounds." << endl;
	return 1;
}
