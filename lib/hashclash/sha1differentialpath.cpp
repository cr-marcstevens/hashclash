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
#include <stdexcept>

#include "rng.hpp"
#include "booleanfunction.hpp"
#include "sha1differentialpath.hpp"
#include "sha1detail.hpp"

using namespace std;

namespace hashclash {

	uint32 bf_simplify(uint32 a, uint32 b, uint32 c) {
		return 0;
	}
	booleanfunction SHA1_F1_data(sha1_f1, "SHA1_F1"), SHA1_F2_data(sha1_f2, "SHA1_F2")
		, SHA1_F3_data(sha1_f3, "SHA1_F3"), SHA1_F4_data(sha1_f4, "SHA1_F4"), BF_simplify(bf_simplify, "BF_simplify");

	void show_path(const sha1differentialpath& path, ostream& o)
	{
		for (int t = path.tbegin(); t < path.tend(); ++t)
		{
			o << "Q" << t << ":\t" << path[t];
			if (t-4 >= path.tbegin() && t+1 < path.tend() && t >= 0 && t < 80) 
			{
				o << path.getme(t);
				vector<unsigned> ambiguous, impossible;
				booleanfunction* F = 0;
				if (t < 20) F = & SHA1_F1_data;
				else if (t < 40) F = & SHA1_F2_data;
				else if (t < 60) F = & SHA1_F3_data;
				else F = & SHA1_F4_data;
				uint32 dF = 0;
				for (unsigned b = 0; b < 32; ++b)
				{
					bitcondition qtm1b = path(t-1,b); if (qtm1b == bc_prev || qtm1b == bc_prevn) qtm1b = bc_constant;
					bf_outcome outcome = F->outcome(qtm1b, path(t-2,((b+2)&31)), path(t-3,((b+2)&31)));
					if (outcome.size() == 1) {
						if (outcome[0] == bc_plus) 			dF += 1<<b;
						else if (outcome[0] == bc_minus)	dF -= 1<<b;
					} else {
						if (outcome.size() == 0) 
							impossible.push_back(b);
						else {
							if (b == 31 && outcome.size() == 2 && outcome(0,31)==outcome(1,31)) {
								dF += 1<<31;
								continue;
							}
							ambiguous.push_back(b);
						}
					}
				}
				uint32 dQtp1 = dF + path.getme(t).adddiff() + path[t].getsdr().rotate_left(5).adddiff() + path[t-4].getsdr().rotate_left(30).adddiff();
				if (dQtp1 == path[t+1].diff())
					o << " ok";
				else
					o << " bad(" << naf(dQtp1 - path[t+1].diff()) << ")";
				if (ambiguous.size() > 0)
				{
					o << " amb:" << ambiguous[0];
					for (unsigned i = 1; i < ambiguous.size(); ++i)
						o << "," << ambiguous[i];
				}
				if (impossible.size() > 0)
				{
					o << " imp:" << impossible[0];
					for (unsigned i = 1; i < impossible.size(); ++i)
						o << "," << impossible[i];
				}
			} else
				o << sdr();
			o << endl;
		}
	}

	bool test_path(const sha1differentialpath& path)
	{
		for (int t = path.tbegin(); t < path.tend(); ++t)
		{
			if (t-4 >= path.tbegin() && t+1 < path.tend() && t >= 0 && t < 80) 
			{
				vector<unsigned> ambiguous, impossible;
				booleanfunction* F = 0;
				if (t < 20) F = & SHA1_F1_data;
				else if (t < 40) F = & SHA1_F2_data;
				else if (t < 60) F = & SHA1_F3_data;
				else F = & SHA1_F4_data;
				uint32 dF = 0;
				for (unsigned b = 0; b < 32; ++b)
				{
					bitcondition qtm1b = path(t-1,b); if (qtm1b == bc_prev || qtm1b == bc_prevn) qtm1b = bc_constant;
					bf_outcome outcome = F->outcome(qtm1b, path(t-2,((b+2)&31)), path(t-3,((b+2)&31)));
					if (outcome.size() == 1) {
						if (outcome[0] == bc_plus) 			dF += 1<<b;
						else if (outcome[0] == bc_minus)	dF -= 1<<b;
					} else {
						if (outcome.size() == 0) 
							impossible.push_back(b);
						else {
							if (b == 31 && outcome.size() == 2 && !outcome.constant())
							{
								dF += 1<<31;
								continue;
							}
							ambiguous.push_back(b);
						}
					}
				}
				uint32 dQtp1 = dF + path.getme(t).adddiff() + path[t].getsdr().rotate_left(5).adddiff() + path[t-4].getsdr().rotate_left(30).adddiff();				
				if (dQtp1 != path[t+1].diff() || ambiguous.size() > 0 || impossible.size() > 0)
					return false;
			}
		}
		return true;
	}

	double deep_analysis_path(const sha1differentialpath& sha1path, const uint32 dQt[80], unsigned tbegin, unsigned tend)
	{

		return 0;
	}

	void cleanup_path(sha1differentialpath& path)
	{
		sha1differentialpath backup;		
		bf_conditions backcond;
		bf_outcome outcome;
		for (int t = path.tbegin()+2; t < path.tend(); ++t)
		{
			for (unsigned b = 0; b < 32; ++b) {
				backcond = BF_simplify.forwardconditions(path(t,b), path(t-1,b), path(t-2,b), bc_constant);
				switch (path(t,b)) {
				case bc_next:
				case bc_nextn:
				case bc_next2:
				case bc_next2n:
					backcond.first = path(t,b);
				}
				switch (path(t-1,b)) {
				case bc_next2:
				case bc_next2n:
					backcond.second = path(t-1,b);
				}
				path.setbitcondition(t,b,backcond.first);
				path.setbitcondition(t-1,b,backcond.second);
				path.setbitcondition(t-2,b,backcond.third);
			}
		}
		for (int t = path.tend()-1; t >= path.tbegin()+2; --t)
		{
			for (unsigned b = 0; b < 32; ++b) {
				backcond = BF_simplify.backwardconditions(path(t,b), path(t-1,b), path(t-2,b), bc_constant);
				switch (path(t-2,b)) {
				case bc_prev2:
				case bc_prev2n:
					backcond.third = path(t-2,b);
				}
				path.setbitcondition(t,b,backcond.first);
				path.setbitcondition(t-1,b,backcond.second);
				path.setbitcondition(t-2,b,backcond.third);
			}
		}
	}

} // namespace hashclash
