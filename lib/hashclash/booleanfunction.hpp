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

#ifndef HASHCLASH_BOOLEANFUNCTION_HPP
#define HASHCLASH_BOOLEANFUNCTION_HPP

#include <vector>
#include <string>
#include <set>
#include <map>

#include <boost/function.hpp>

#include "types.hpp"
#include "sdr.hpp"
#include "md5detail.hpp"
#include "conditions.hpp"

namespace hashclash {

	typedef triple<bitcondition, bitcondition, bitcondition> bf_conditions;

	struct bf_outcome {
		unsigned char c;

		static const unsigned char fconstant = 1;
		static const unsigned char fplus = 2;
		static const unsigned char fminus = 4;
		bf_outcome(): c(0) {}
		bf_outcome(unsigned char k): c(k) {}
		bf_outcome& operator=(const bf_outcome& r)
		{ c = r.c; return *this; }
		bf_outcome& operator=(unsigned char k)
		{ c = k; return *this; }

		bool constant() const
		{ return (c&fconstant) != 0; }
		bool plus() const
		{ return (c&fplus) != 0; }
		bool minus() const
		{ return (c&fminus) != 0; }
		unsigned size() const
		{
			unsigned w = 0;
			if (constant()) ++w;
			if (plus()) ++w;
			if (minus()) ++w;
			return w;
		}
		uint32 operator()(unsigned index, unsigned b) const
		{
			if (index == 0) {
				if (constant()) return 0;
				if (plus()) return 0+(1<<b);
				return 0-(1<<b);
			} else if (index == 1) {
				if (constant())	{
					if (plus()) return 0+(1<<b);
					return 0-(1<<b);
				} else 
					return 0-(1<<b);
			} else
				return 0-(1<<b);
		}
		bitcondition operator[](unsigned index) const
		{
			if (index == 0) {
				if (constant()) return bc_constant;
				if (plus()) return bc_plus;
				return bc_minus;
			} else if (index == 1) {
				if (constant())	{
					if (plus()) return bc_plus;
					return bc_minus;
				} else 
					return bc_minus;
			} else 
				return bc_minus;
		}
	};
	
	class booleanfunction {
	public:
		booleanfunction(const boost::function<uint32(uint32,uint32,uint32)>& F, const std::string& description = "");
		
		const bf_outcome& outcome(bitcondition input1, bitcondition input2, bitcondition input3)
		{ return outcome_table[(input1<<8) + (input2<<4) + input3]; }
		const bf_outcome& outcome(const bf_conditions& c)
		{ return outcome(c.first, c.second, c.third); }

		const bf_conditions& forwardconditions(bitcondition input1, bitcondition input2, bitcondition input3, bitcondition outcome)
		{ return forward_table[(outcome<<12)+(input1<<8)+(input2<<4)+input3]; }
		const bf_conditions& forwardconditions(const bf_conditions& c, bitcondition outcome)
		{ return forwardconditions(c.first, c.second, c.third, outcome); }

		const bf_conditions& backwardconditions(bitcondition input1, bitcondition input2, bitcondition input3, bitcondition outcome)
		{ return backward_table[(outcome<<12)+(input1<<8)+(input2<<4)+input3]; }
		const bf_conditions& backwardconditions(const bf_conditions& c, bitcondition outcome)
		{ return backwardconditions(c.first, c.second, c.third, outcome); }

		uint32 F(uint32 input1, uint32 input2, uint32 input3) const
		{ return f(input1, input2, input3); }

		const std::string& description() const
		{ return f_description; }
		
	private:
		boost::function<uint32(uint32,uint32,uint32)> f;
		std::string f_description;

		std::vector<bf_outcome> outcome_table;
		std::vector<bf_conditions> forward_table;
		std::vector<bf_conditions> backward_table;

		std::map<bf_conditions, std::set<uint32> > conds_to_values;

		void find_booleanfunction_outcomes(const bf_conditions& cond, bf_outcome& outcomes, std::set<uint32>& values);
		bf_conditions preferred_conditions(const bf_conditions& cond, std::vector<bf_conditions>& vec_conds);
	};

	extern booleanfunction MD5_F_data, MD5_G_data, MD5_H_data, MD5_I_data;

} // namespace hashclash

#endif //HASHCLASH_BOOLEANFUNCTION_HPP
