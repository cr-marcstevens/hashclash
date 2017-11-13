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

#include <iostream>

#include "booleanfunction.hpp"

using namespace std;

namespace hashclash {

	booleanfunction MD5_F_data(md5_ff, "MD5_FF"), MD5_G_data(md5_gg, "MD5_GG")
		, MD5_H_data(md5_hh, "MD5_HH"), MD5_I_data(md5_ii, "MD5_II");

	inline bool set_containedin_set(const set<uint32>& needle, const set<uint32>& haystack)
	{
		for (set<uint32>::const_iterator cit = needle.begin(); cit != needle.end(); ++cit)
			if (haystack.find(*cit) == haystack.end())
				return false;
		return true;
	}

	inline void join_sets(set<uint32>& bigset, const set<uint32>& smallset)
	{
		set<uint32>::const_iterator cit = smallset.begin();
		set<uint32>::const_iterator citend = smallset.end();
		for (; cit != citend; ++cit)
			bigset.insert(*cit);
	}

	void booleanfunction::find_booleanfunction_outcomes(const bf_conditions& cond, bf_outcome& outcomes, set<uint32>& values)
	{
		outcomes.c = 0;
		values.clear();

		for (uint32 c0 = 0; c0 <= 1; ++c0)
		for (uint32 c1 = 0; c1 <= 1; ++c1)
		for (uint32 b0 = 0; b0 <= 1; ++b0)
		for (uint32 b1 = 0; b1 <= 1; ++b1)
		for (uint32 a0 = 0; a0 <= 1; ++a0)
		for (uint32 a1 = 0; a1 <= 1; ++a1)
		{
			switch (cond.first) {
				case bc_next:
				case bc_nextn:
				case bc_next2:
				case bc_next2n:
				case bc_or2b:
				case bc_constant:
					if (!(a1 == a0)) continue; break;
				case bc_plus:
					if (!(a1 == a0+1)) continue; break;
				case bc_minus:
					if (!(a1 == a0-1)) continue; break;
				case bc_zero:
					if (!(a0 == 0 && a1 == 0)) continue; break;
				case bc_one:
					if (!(a0 == 1 && a1 == 1)) continue; break;
				case bc_prev:
					if (!(a1 == a0 && a0 == b0)) continue; break;
				case bc_prevn:
					if (!(a1 == a0 && a0 != b0)) continue; break;
				case bc_prev2:
					if (!(a1 == a0 && a0 == c0)) continue; break;
				case bc_prev2n:
					if (!(a1 == a0 && a0 != c0)) continue; break;
				case bc_or2:
					if (!(a1 == a0 && (a0 == 1 || c0 == 0))) continue; break;
				default:
					throw;
			}
			switch (cond.second) {
				case bc_or2b:
				case bc_next2:
				case bc_next2n:
				case bc_or2:
				case bc_prev2:
				case bc_prev2n:
				case bc_constant:
					if (!(b1 == b0)) continue; break;
				case bc_plus:
					if (!(b1 == b0+1)) continue; break;
				case bc_minus:
					if (!(b1 == b0-1)) continue; break;
				case bc_zero:
					if (!(b0 == 0 && b1 == 0)) continue; break;
				case bc_one:
					if (!(b0 == 1 && b1 == 1)) continue; break;
				case bc_prev:
					if (!(b1 == b0 && b0 == c0)) continue; break;
				case bc_prevn:
					if (!(b1 == b0 && b0 != c0)) continue; break;
				case bc_next:
					if (!(b1 == b0 && b0 == a0)) continue; break;
				case bc_nextn:
					if (!(b1 == b0 && b0 != a0)) continue; break;
				default:
					throw;
			}
			switch (cond.third) {
				case bc_or2:
				case bc_prev:
				case bc_prevn:
				case bc_prev2:
				case bc_prev2n:
				case bc_constant:
					if (!(c1 == c0)) continue; break;
				case bc_plus:
					if (!(c1 == c0+1)) continue; break;
				case bc_minus:
					if (!(c1 == c0-1)) continue; break;
				case bc_zero:
					if (!(c0 == 0 && c1 == 0)) continue; break;
				case bc_one:
					if (!(c0 == 1 && c1 == 1)) continue; break;
				case bc_next:
					if (!(c1 == c0 && c0 == b0)) continue; break;
				case bc_nextn:
					if (!(c1 == c0 && c0 != b0)) continue; break;
				case bc_next2:
					if (!(c1 == c0 && c0 == a0)) continue; break;
				case bc_next2n:
					if (!(c1 == c0 && c0 != a0)) continue; break;
				case bc_or2b:
					if (!(c1 == c0 && (a0 == 1 || c0 == 0))) continue; break;
				default:
					throw;
			}

			values.insert(a0 | (a1<<1) | (b0<<2) | (b1<<3) | (c0<<4) | (c1<<5));

			uint32 v0 = (f(a0,b0,c0)&1);
			uint32 v1 = (f(a1,b1,c1)&1);
			if (v0 == v1)
				outcomes.c |= bf_outcome::fconstant;				
			else if (v1 == 1)
				outcomes.c |= bf_outcome::fplus;
			else
				outcomes.c |= bf_outcome::fminus;			
		}
	}

	bf_conditions booleanfunction::preferred_conditions(const bf_conditions& conds, std::vector<bf_conditions>& vec_conds)
	{
		if (vec_conds.size() == 0)
			throw std::out_of_range("preferred_conditions(): no conditions");

		// eliminate non-compatible conditions
		// compatible: U_def \subseteq U_abc
		set<uint32>& values = conds_to_values[conds];
		unsigned j = 0;
		unsigned s = 0;
		while (j < vec_conds.size())
		{
			set<uint32>& values2 = conds_to_values[vec_conds[j]];
			if (!set_containedin_set(values2,values))
			{
				vec_conds[j] = vec_conds[vec_conds.size()-1];
				vec_conds.pop_back();
			} else {
				if (values2.size() > s)
					s = unsigned(values2.size());
				++j;
			}
		}
		// eliminate non-maximal conditions: 
		// maximal: |U_def|=max_ghi |U_ghi|
		j = 0;
		while (j < vec_conds.size())
			if (conds_to_values[vec_conds[j]].size() != s)
			{
				vec_conds[j] = vec_conds[vec_conds.size()-1];
				vec_conds.pop_back();
			} else
				++j;

		if (vec_conds.size() == 0)
			throw std::out_of_range("preferred_conditions(): no conditions");

		unsigned bestw = 1<<31;
		unsigned bestindex = 0;
		for (unsigned i = 0; i < vec_conds.size(); ++i)
		{
			unsigned w = 0;
			bitcondition cond0 = vec_conds[i].first;
			bitcondition cond1 = vec_conds[i].second;
			bitcondition cond2 = vec_conds[i].third;
			if (cond0 != bc_constant) w+=16;
			if (cond1 != bc_constant) w+=16;
			if (cond2 != bc_constant) w+=16;
			if (!isdirect(cond0)) w+=4;
			if (!isdirect(cond1)) w+=4;
			if (!isdirect(cond2)) w+=4;
			if (isindirect2(cond0)) w+=1;
			if (isindirect2(cond1)) w+=1;
			if (isindirect2(cond2)) w+=1;
			if (w < bestw)
			{
				bestw = w;
				bestindex = i;
			}
		}
		return vec_conds[bestindex];
	}


	booleanfunction::booleanfunction(const boost::function<uint32(uint32,uint32,uint32)>& F, const std::string& description)
		: f(F), f_description(description)
	{
		outcome_table.resize(bc_max * 16 * 16);
		forward_table.resize(3 * 16 * 16 * 16);
		backward_table.resize(3 * 16 * 16 * 16);

		map<bitcondition, vector<bf_conditions> > 
			outcome_to_conds_forward, outcome_to_conds_backward;

		bf_conditions conds, conds2;
		set<uint32> values, values2;
		bf_outcome outcomes; 
		for (unsigned cond0 = 0; cond0 < bc_max; ++cond0)
		if (isbackward(bitcondition(cond0)))
			for (unsigned cond1 = 0; cond1 < bc_max; ++cond1)
			if (!isindirect2(bitcondition(cond1)))
				for (unsigned cond2 = 0; cond2 < bc_max; ++cond2)
				if (isforward(bitcondition(cond2)))
				{
					conds.first = bitcondition(cond0);
					conds.second = bitcondition(cond1);
					conds.third = bitcondition(cond2);
					find_booleanfunction_outcomes(conds, outcomes, values);
					outcome_table[(cond0<<8) + (cond1<<4) + cond2] = outcomes;
					conds_to_values[conds] = values;
					if (outcomes.size() == 1)
					{
						if (isforward(conds.first) && isforward(conds.second) && isforward(conds.third))
							outcome_to_conds_forward[outcomes[0]].push_back(conds);
						if (isbackward(conds.first) && isbackward(conds.second) && isbackward(conds.third))
							outcome_to_conds_backward[outcomes[0]].push_back(conds);
					}
				}
	
		vector<bf_conditions> vec_conds;
		for (unsigned cond0 = 0; cond0 < bc_max; ++cond0)
		if (isbackward(bitcondition(cond0)))
			for (unsigned cond1 = 0; cond1 < bc_max; ++cond1)
			if (!isindirect2(bitcondition(cond1)))
				for (unsigned cond2 = 0; cond2 < bc_max; ++cond2)
				if (isforward(bitcondition(cond2)))
				{
					conds.first = bitcondition(cond0);
					conds.second = bitcondition(cond1);
					conds.third = bitcondition(cond2);
					outcomes = outcome_table[(cond0<<8) + (cond1<<4) + cond2];
					values = conds_to_values[conds];
					for (unsigned i = 0; i < outcomes.size(); ++i)
					{
						bitcondition bfdiff = outcomes[i];
						vec_conds = outcome_to_conds_forward[bfdiff];
						conds2 = preferred_conditions(conds, vec_conds);
						forward_table[(bfdiff<<12)+(cond0<<8) + (cond1<<4) + cond2] = conds2;

						vec_conds = outcome_to_conds_backward[bfdiff];
						conds2 = preferred_conditions(conds, vec_conds);
						backward_table[(bfdiff<<12)+(cond0<<8) + (cond1<<4) + cond2] = conds2;
					}
				}
		for (unsigned cond0a = 0; cond0a < bc_max; ++cond0a)
			for (unsigned cond1a = 0; cond1a < bc_max; ++cond1a)
				for (unsigned cond2a = 0; cond2a < bc_max; ++cond2a)
				{
					unsigned cond0 = cond0a, cond1 = cond1a, cond2 = cond2a;
					if (!isbackward(bitcondition(cond0)))
						cond0 = bc_constant;
					if (isindirect2(bitcondition(cond1)))
						cond1 = bc_constant;
					if (!isforward(bitcondition(cond2)))
						cond2 = bc_constant;
					if (cond0==cond0a && cond1==cond1a && cond2==cond2a)
						continue;

					outcomes = outcome_table[(cond0<<8) + (cond1<<4) + cond2];
					outcome_table[(cond0a<<8) + (cond1a<<4) + cond2a] = outcomes;
					for (unsigned i = 0; i < outcomes.size(); ++i)
					{
						bitcondition bfdiff = outcomes[i];
						forward_table[(bfdiff<<12)+(cond0a<<8) + (cond1a<<4) + cond2a] 
							= forward_table[(bfdiff<<12)+(cond0<<8) + (cond1<<4) + cond2];
						backward_table[(bfdiff<<12)+(cond0a<<8) + (cond1a<<4) + cond2a] 
							= backward_table[(bfdiff<<12)+(cond0<<8) + (cond1<<4) + cond2];
					}
				}


	}
} // namespace hashclash
