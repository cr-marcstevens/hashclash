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
#include "differentialpath.hpp"

using namespace std;

namespace hashclash {

	double check_rotation(uint32 dR, uint32 dT, unsigned n, const wordconditions& Qt, const wordconditions& Qtp1, unsigned loopcount)
	{
		uint32 Q0set0 = Qt.set0(), Q0set1 = Qt.set1();
		uint32 Q1set0 = Qtp1.set0(), Q1set1 = Qtp1.set1();
		uint32 prev = Qtp1.prev() | Qt.next()
			, prevn = Qtp1.prevn() | Qt.nextn();

		unsigned okcount = 0;
		for (unsigned count = 0; count < loopcount; ++count)
		{
			uint32 Q0a = (xrng64() & Q0set0) | Q0set1;
			uint32 Q1a = ((xrng64() & Q1set0 & ~prev) | Q1set1 | prevn | (Q0a & prev)) ^ (Q0a & prevn);
			uint32 Ra = Q1a - Q0a;
			uint32 Ta = rotate_right(Ra, n);
			uint32 Tb = Ta + dT;
			uint32 Rb = rotate_left(Tb, n);
			if (Rb-Ra == dR) ++okcount;
		}
		return double(okcount)/double(loopcount);
	}

	void show_path(const differentialpath& path, const uint32 blockdiff[], ostream& o)
	{
		for (int t = path.tbegin(); t < path.tend(); ++t)
		{
			o << "Q" << t << ":\t" << path[t];
			if (t-3 >= path.tbegin() && t+1 < path.tend() && t >= 0 && t < 64) 
			{
				vector<unsigned> ambiguous, impossible;
				booleanfunction* F = 0;
				if (t < 16) F = & MD5_F_data;
				else if (t < 32) F = & MD5_G_data;
				else if (t < 48) F = & MD5_H_data;
				else F = & MD5_I_data;
				uint32 dF = 0;
				for (unsigned b = 0; b < 32; ++b)
				{
					bf_outcome outcome = F->outcome(path(t,b), path(t-1,b), path(t-2,b));
					if (outcome.size() == 1) {
						if (outcome[0] == bc_plus) 			dF += 1<<b;
						else if (outcome[0] == bc_minus)	dF -= 1<<b;
					} else {
						if (outcome.size() == 0) 
							impossible.push_back(b);
						else
							ambiguous.push_back(b);
					}
				}
				uint32 dQtm3 = path[t-3].diff();
				uint32 dR = path[t+1].diff() - path[t].diff();
				uint32 dT2 = naf(dR).rotate_right(md5_rc[t]).adddiff();
				uint32 dT1 = dQtm3 + dF + blockdiff[md5_wt[t]];
				uint32 dR2 = naf(dT1).rotate_left(md5_rc[t]).adddiff();
				double p = check_rotation(dR, dT1, md5_rc[t], path[t], path[t+1]);
				if (dT1 == dT2 || dR == dR2 || p>0)
					o << " ok p=" << p;
				else
					o << " bad p=" << p;
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
			}
			o << endl;
		}
	}

	double test_path(const differentialpath& path, const uint32 blockdiff[])
	{
		double totp = 1;
		for (int t = path.tbegin(); t < path.tend(); ++t)
		{			
			if (t-3 >= path.tbegin() && t+1 < path.tend() && t >= 0 && t < 64) 
			{
				booleanfunction* F = 0;
				if (t < 16) F = & MD5_F_data;
				else if (t < 32) F = & MD5_G_data;
				else if (t < 48) F = & MD5_H_data;
				else F = & MD5_I_data;
				uint32 dF = 0;
				for (unsigned b = 0; b < 32; ++b)
				{
					bf_outcome outcome = F->outcome(path(t,b), path(t-1,b), path(t-2,b));
					if (outcome.size() == 1) {
						if (outcome[0] == bc_plus) 			dF += 1<<b;
						else if (outcome[0] == bc_minus)	dF -= 1<<b;
					} else
						return 0;
				}
				uint32 dQtm3 = path[t-3].diff();
				uint32 dR = path[t+1].diff() - path[t].diff();
				uint32 dT = dQtm3 + dF + blockdiff[md5_wt[t]];
				totp *= check_rotation(dR, dT, md5_rc[t], path[t], path[t+1]);
				if (totp == 0) return 0;
			}
		}
		return totp;
	}

	unsigned totaltunnelstrength(const differentialpath& path)
	{
		unsigned totalstrength = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			// best tunnel first
			bitcondition q5b = path(5,b);
			bitcondition q4b = path(4,b);
			bitcondition q3b = path(3,b);
			bitcondition q14b = path(14,b);
			bitcondition q15b = path(15,b);
			bitcondition q16b = path(16,b);
			if (1)
			if (q14b == bc_constant && q3b == bc_constant && q15b != bc_prev && q15b != bc_prevn && q4b != bc_prev && q4b != bc_prevn) {
				if (q4b == bc_constant || q4b == bc_zero || q4b == bc_plus)
				if (q5b == bc_constant || q5b == bc_prevn || q5b == bc_one || q5b == bc_minus)
				if (q15b == bc_constant || q15b == bc_zero || q15b == bc_plus)
				if (q16b == bc_constant || q16b == bc_prev || q16b == bc_zero || q16b == bc_plus)
				{
					++totalstrength;
					q3b = bc_zero; // disable further tunnels based on Q3
					q4b = bc_zero; //  " on Q4
					q5b = bc_one;  //  " on Q5
				}
			}

			bitcondition q8b = path(8,b);
			bitcondition q12b = path(12,b);
			bitcondition q13b = path(13,b);			
			bitcondition q9b = path(9,b);
			bitcondition q10b = path(10,b);
			bitcondition q11b = path(11,b);
			if (q8b == bc_constant && q9b != bc_prev && q9b != bc_prevn
				&& q12b == bc_constant && q13b != bc_prev && q13b != bc_prevn
				&& (q10b == bc_constant || q10b == bc_minus || q10b == bc_one))
			{
				// Q8Q12m15 tunnel
				++totalstrength;
				// Q9m10 tunnel is only other possible
				if (q9b == bc_constant
					&& (q10b == bc_constant || q10b == bc_one || q10b == bc_minus)
					&& (q11b == bc_constant || q11b == bc_one || q11b == bc_minus || q11b == bc_prev))
					++totalstrength;
			} else {
				if (q9b == bc_constant 
					&& (q10b == bc_constant || q10b == bc_zero || q10b == bc_plus)
					&& (q11b == bc_constant || q11b == bc_one || q11b == bc_minus || q11b == bc_prevn))
					++totalstrength;
				else if (q9b == bc_constant
					&& (q10b == bc_constant || q10b == bc_one || q10b == bc_minus)
					&& (q11b == bc_constant || q11b == bc_one || q11b == bc_minus || q11b == bc_prev))
					++totalstrength;
				else if (q10b == bc_constant
					&& (q11b == bc_constant || q11b == bc_zero || q11b == bc_plus))
					++totalstrength;
			}

			if (q14b == bc_constant && q3b == bc_constant
					&& q15b != bc_prev && q15b != bc_prevn
					&& q4b != bc_prev && q4b != bc_prevn)
			{
				if (q15b == bc_constant || q16b == bc_constant || q16b == bc_prev
					|| ((q15b == bc_one || q15b == bc_minus) && (q16b == bc_one || q16b == bc_minus))
					|| ((q15b == bc_zero || q15b == bc_plus) && (q16b == bc_zero || q16b == bc_plus))
					)
					++totalstrength;
			}

			bitcondition q6b = path(6,b);
			if (q4b == bc_constant 
				&& (q5b == bc_constant || q5b == bc_zero || q5b == bc_plus)
				&& (q6b == bc_constant || q6b == bc_one || q6b == bc_minus || q6b == bc_prevn))
				++totalstrength;
			else if (q4b == bc_constant
				&& (q5b == bc_constant || q5b == bc_one || q5b == bc_minus)
				&& (q6b == bc_constant || q6b == bc_one || q6b == bc_minus || q6b == bc_prev))
				++totalstrength;
			else if (q5b == bc_constant
				&& (q6b == bc_constant || q6b == bc_zero || q6b == bc_plus))
				++totalstrength;

		}
		uint32 Q0Q1eq = (path[0].prev()|~(path[0].set1()^path[-1].set1()));
		uint32 Q1free = ~path[1].mask();
		totalstrength += (hw(Q0Q1eq & Q1free & ~path[2].mask())>>1) 
						+ hw(Q0Q1eq & Q1free & path[2].prev());
		return totalstrength;
	}

	bool check_rotation_fast(uint32 dR, uint32 dT, unsigned n, const wordconditions& Qt, const wordconditions& Qtp1, unsigned loopcount)
	{
		uint32 Q0set0 = Qt.set0(), Q0set1 = Qt.set1();
		uint32 Q1set0 = Qtp1.set0(), Q1set1 = Qtp1.set1();
		uint32 prev = Qtp1.prev() | Qt.next()
			, prevn = Qtp1.prevn() | Qt.nextn();

		unsigned okcount = 0;
		for (unsigned count = 0; count < loopcount; ++count)
		{
			uint32 Q0a = (xrng64() & Q0set0) | Q0set1;
			uint32 Q1a = ((xrng64() & Q1set0 & ~prev) | Q1set1 | prevn | (Q0a & prev)) ^ (Q0a & prevn);
			uint32 Ra = Q1a - Q0a;
			uint32 Ta = rotate_right(Ra, n);
			uint32 Tb = Ta + dT;
			uint32 Rb = rotate_left(Tb, n);
			if (Rb-Ra == dR) {
				++okcount;
				if (okcount >= (loopcount>>6))
					return true;
			}
		}
		return false;
	}

	bool test_path_fast(const differentialpath& path, const uint32 blockdiff[])
	{
		for (int t = path.tbegin(); t < path.tend(); ++t)
		{			
			if (t-3 >= path.tbegin() && t+1 < path.tend() && t >= 0 && t < 64) 
			{
				vector<unsigned> ambiguous, impossible;
				booleanfunction* F = 0;
				if (t < 16) F = & MD5_F_data;
				else if (t < 32) F = & MD5_G_data;
				else if (t < 48) F = & MD5_H_data;
				else F = & MD5_I_data;
				uint32 dF = 0;
				for (unsigned b = 0; b < 32; ++b)
				{
					bf_outcome outcome = F->outcome(path(t,b), path(t-1,b), path(t-2,b));
					if (outcome.size() == 1) {
						if (outcome[0] == bc_plus) 			dF += 1<<b;
						else if (outcome[0] == bc_minus)	dF -= 1<<b;
					} else
						return false;
				}
				uint32 dQtm3 = path[t-3].diff();
				uint32 dR = path[t+1].diff() - path[t].diff();
				uint32 dT = dQtm3 + dF + blockdiff[md5_wt[t]];
				if (false == check_rotation_fast(dR, dT, md5_rc[t], path[t], path[t+1]))
					return false;
			}
		}
		return true;
	}

	void cleanup(differentialpath& path)
	{
		differentialpath backup;		
		bf_conditions backcond;
		bf_outcome outcome;
		backup = path;
		for (int t = path.tbegin()+3; t <= path.tend()-2; ++t)
		{
			booleanfunction* F = 0;
			if (t < 16) F = & MD5_F_data;
			else if (t < 32) F = & MD5_G_data;
			else if (t < 48) F = & MD5_H_data;
			else F = & MD5_I_data;
			for (unsigned b = 0; b < 32; ++b)
			{
				outcome = F->outcome(backup(t,b), backup(t-1,b), backup(t-2,b));
				if (outcome.size() != 1) throw std::runtime_error("hashclash::cleanup(differentialpath&): path is ambiguous!");
				backcond = F->forwardconditions(path(t,b), path(t-1,b), path(t-2,b), outcome[0]);
				path.setbitcondition(t, b, backcond.first);
				path.setbitcondition(t-1, b, backcond.second);
				path.setbitcondition(t-2, b, backcond.third);
			}
		}
		for (int t = path.tend()-2; t >= path.tbegin()+3; --t)
		{
			booleanfunction* F = 0;
			if (t < 16) F = & MD5_F_data;
			else if (t < 32) F = & MD5_G_data;
			else if (t < 48) F = & MD5_H_data;
			else F = & MD5_I_data;
			for (unsigned b = 0; b < 32; ++b)
			{
				outcome = F->outcome(backup(t,b), backup(t-1,b), backup(t-2,b));
				if (outcome.size() != 1) throw std::runtime_error("hashclash::cleanup(differentialpath&): path is ambiguous!");
				backcond = F->backwardconditions(path(t,b), path(t-1,b), path(t-2,b), outcome[0]);
				path.setbitcondition(t, b, backcond.first);
				path.setbitcondition(t-1, b, backcond.second);
				path.setbitcondition(t-2, b, backcond.third);
			}
		}
	}

















	/*** enhance path functions ***/

	typedef triple<wordconditions,wordconditions,unsigned> soltriple;
	struct solutions_less
		: public std::binary_function<soltriple, soltriple, bool>
	{
		bool operator()(const soltriple& _Left, const soltriple& _Right) const
		{ return _Left.third < _Right.third; }
	};

	bool isstronger(const soltriple& lh, const soltriple& rh)
	{
		uint32 lht0 = ~lh.first.set0(), lht1 = lh.first.set1();
		uint32 lhtf = ~lh.first.mask() | lh.first.prev() | lh.first.prevn();
		uint32 lhtt0 = ~lh.second.set0(), lhtt1 = lh.second.set1();
		uint32 lhttp = lh.second.prev(), lhttpn = lh.second.prevn();
		uint32 lhttf = ~lh.second.mask();

		uint32 rht0 = ~rh.first.set0(), rht1 = rh.first.set1();
		uint32 rhtf = ~rh.first.mask() | rh.first.prev() | rh.first.prevn();
		uint32 rhtt0 = ~rh.second.set0(), rhtt1 = rh.second.set1();
		uint32 rhttp = rh.second.prev(), rhttpn = rh.second.prevn();
		uint32 rhttf = ~rh.second.mask();

		// rh not free => lh not free
		if (rhttf != (rhttf | lhttf)) return false;
		if (rhtf != (rhtf | lhtf)) return false;
		
		// rh=0 => lh=0
		if (rht0 != (rht0 & lht0)) return false;
		if (rht1 != (rht1 & lht1)) return false;
		// rh=1 => lh=1
		if (rhtt0 != (rhtt0 & lhtt0)) return false;
		if (rhtt1 != (rhtt1 & lhtt1)) return false;
		
		// rh=^ => lh=0,1,^
		if (rhttp != (rhttp & (lhtt0 | lhtt1 | lhttp))) return false;
		// rh=! => lh=0,1,!
		if (rhttpn != (rhttpn & (lhtt0 | lhtt1 | lhttpn))) return false;

		// rh=lh => false
		if (lh == rh)
			return false;

		return true;
	}

	bool errorinrotation32(const wordconditions& Qt, const wordconditions& Qtp1, 
				uint32 dT, uint32 dR, unsigned rc, uint32 loopcount = (1<<12))
	{
		vector< pair<uint32,double> > rotdiffs;
		rotate_difference(dT, rc, rotdiffs);
		for (unsigned i = 0; i < rotdiffs.size(); ++i)
			if (naf(rotdiffs[i].first-dR).get(0) != 0) {
				if (0 < check_rotation(rotdiffs[i].first, dT, rc, Qt, Qtp1, loopcount))
					return true;
			}
		return false;	
	}

	bool errorinrotationRC(const wordconditions& Qt, const wordconditions& Qtp1, 
				uint32 dT, uint32 dR, unsigned rc, uint32 loopcount = (1<<12))
	{
		vector< pair<uint32,double> > rotdiffs;
		rotate_difference(dT, rc, rotdiffs);
		for (unsigned i = 0; i < rotdiffs.size(); ++i)
			if (naf(rotdiffs[i].first-dR).get(rc) != 0) {
				if (0 < check_rotation(rotdiffs[i].first, dT, rc, Qt, Qtp1, loopcount))
					return true;
			}
		return false;	
	}

	void findsolutions32(const wordconditions& Qt, const wordconditions& Qtp1,
				uint32 dT, uint32 dR, unsigned rc, 
				vector< triple<wordconditions,wordconditions,unsigned> >& solutions,
				uint32 loopcount = (1<<12))
	{
		solutions.clear();

		// first determine bounds on the range of bits we are going to look at
		vector< vector< triple<bitcondition,bitcondition,unsigned> > > vnbc(32);
		double prot = check_rotation(dR, dT, rc, Qt, Qtp1, loopcount);
		int bmin = int(-(log(prot)/log(double(2))))+3;
		
		int bit = 31;
		uint64 cnt = 1;
	    while (solutions.size() == 0)
	    {		
		// we start at bit 31 and work downwards
		// for every free bit in Q_t,Q_{t+1} we store possible extra conditions
		bool Qttp1free = false;
		int bitfree = 0; bit = 31;
		cnt = 1;
		while (bit >= int(rc) && bitfree <= bmin+2 && (bitfree <= bmin || !Qttp1free))
		{
			if (Qt[bit] == bc_constant || Qt[bit] == bc_prev || Qt[bit] == bc_prevn) {
				if (Qtp1[bit] == bc_constant) {
					bitfree+=2;
					Qttp1free = true;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_prevn, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_one, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_zero, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_prev, 1) );
					vnbc[bit].push_back( make_triple(bc_one, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_zero, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_one, bc_zero, 2) );
					vnbc[bit].push_back( make_triple(bc_one, bc_one, 2) );
					vnbc[bit].push_back( make_triple(bc_zero, bc_zero, 2) );
					vnbc[bit].push_back( make_triple(bc_zero, bc_one, 2) );
				} else {
					++bitfree;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(bc_one, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_zero, Qtp1[bit], 1) );
				}
			} else {
				if (Qtp1[bit] == bc_constant) {
					++bitfree;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_one, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_zero, 1) );
				}			
			}
			if (vnbc[bit].size())
				cnt *= vnbc[bit].size();
			--bit;
		}
		
		break;
		if (solutions.size() == 0) {
			++bmin;
			loopcount >>= 1;
			if (loopcount <= (1<<9))
				break;
		}
	    }

		// now we try every combination of these extra conditions
		// and determine which actually are solutions: "no" errors at bit 31 in rotation
		for (uint64 k = 1; k < cnt; ++k)
		{
			wordconditions newQt = Qt, newQtp1 = Qtp1;
			uint64 kk = k;
			unsigned extracond = 0;
			for (int i = 31; i > bit; --i)
			{
				if (vnbc[i].size()) {
					unsigned index = unsigned(kk % vnbc[i].size());
					kk /= vnbc[i].size();
					newQt.set(i, vnbc[i][index].first);
					newQtp1.set(i, vnbc[i][index].second);
					extracond += vnbc[i][index].third;
				}
			}
			if (!errorinrotation32(newQt, newQtp1, dT, dR, rc, loopcount))
				solutions.push_back(make_triple(newQt, newQtp1, extracond));
		}

		// we keep only the basic solutions:
		// we remove all solutions that consist of another solution with added conditions
		unsigned i=0;
		while (i < solutions.size())
		{
			unsigned j = 0;
			while (j < solutions.size())
			{
				if (j != i && isstronger(solutions[i],solutions[j])) {
					swap(solutions[i], solutions[solutions.size()-1]);
					solutions.pop_back();
					break;
				} else
					++j;
			}
			if (j >= solutions.size())
				++i;
		}

		// we sort by number of added conditions so that we try ones with the fewest first
		sort(solutions.begin(), solutions.end(), solutions_less());
	}

	void findsolutionsRC(const wordconditions& Qt, const wordconditions& Qtp1,
				uint32 dT, uint32 dR, unsigned rc, 
				vector< triple<wordconditions,wordconditions,unsigned> >& solutions,
				uint32 loopcount = (1<<12))
	{
		solutions.clear();

		// first determine bounds on the range of bits we are going to look at
		vector< vector< triple<bitcondition,bitcondition,unsigned> > > vnbc(32);
		double prot = check_rotation(dR, dT, rc, Qt, Qtp1, loopcount);
		int bmin = int(-(log(prot)/log(double(2))))+3;
		
		int bit; uint64 cnt;
	    while (solutions.size() == 0)
	    {	
		// we start at bit 31 and work downwards
		// for every free bit in Q_t,Q_{t+1} we store possible extra conditions
		bool Qttp1free = false;
		bit = int(rc)-1;
		int bitfree = 0;
		cnt = 1;
		while (bit >= 0 && bitfree <= bmin+2 && (bitfree <= bmin || !Qttp1free))
		{
			if (Qt[bit] == bc_constant || Qt[bit] == bc_prev || Qt[bit] == bc_prevn) {
				if (Qtp1[bit] == bc_constant) {
					bitfree+=2;
					Qttp1free = true;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_prevn, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_one, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_zero, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_prev, 1) );
					vnbc[bit].push_back( make_triple(bc_one, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_zero, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_one, bc_zero, 2) );
					vnbc[bit].push_back( make_triple(bc_one, bc_one, 2) );
					vnbc[bit].push_back( make_triple(bc_zero, bc_zero, 2) );
					vnbc[bit].push_back( make_triple(bc_zero, bc_one, 2) );
				} else {
					++bitfree;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(bc_one, Qtp1[bit], 1) );
					vnbc[bit].push_back( make_triple(bc_zero, Qtp1[bit], 1) );
				}
			} else {
				if (Qtp1[bit] == bc_constant) {
					++bitfree;
					vnbc[bit].push_back( make_triple(Qt[bit], Qtp1[bit], 0) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_one, 1) );
					vnbc[bit].push_back( make_triple(Qt[bit], bc_zero, 1) );
				}			
			}
			if (vnbc[bit].size())
				cnt *= vnbc[bit].size();
			--bit;
		}

		// now we try every combination of these extra conditions
		// and determine which actually are solutions: "no" errors at bit 31 in rotation
		for (uint64 k = 1; k < cnt; ++k)
		{
			wordconditions newQt = Qt, newQtp1 = Qtp1;
			uint64 kk = k;
			unsigned extracond = 0;
			for (int i = 31; i > bit; --i)
			{
				if (vnbc[i].size()) {
					unsigned index = unsigned(kk % vnbc[i].size());
					kk /= vnbc[i].size();
					newQt.set(i, vnbc[i][index].first);
					newQtp1.set(i, vnbc[i][index].second);
					extracond += vnbc[i][index].third;
				}
			}
			if (!errorinrotationRC(newQt, newQtp1, dT, dR, rc, loopcount))
				solutions.push_back(make_triple(newQt, newQtp1, extracond));
		}
		
		break;
		if (solutions.size() == 0) {
			++bmin;
			loopcount >>= 1;
			if (loopcount <= (1<<9))
				break;
		}
	    }

		// we keep only the basic solutions:
		// we remove all solutions that consist of another solution with added conditions
		unsigned i=0;
		while (i < solutions.size())
		{
			unsigned j = 0;
			while (j < solutions.size())
			{
				if (j != i && isstronger(solutions[i],solutions[j])) {
					swap(solutions[i], solutions[solutions.size()-1]);
					solutions.pop_back();
					break;
				} else
					++j;
			}
			if (j >= solutions.size())
				++i;
		}

		// we sort by number of added conditions so that we try ones with the fewest first
		sort(solutions.begin(), solutions.end(), solutions_less());
	}

	void apply_condition(differentialpath& newpath, int t, unsigned b, bitcondition bc)
	{
		// first insert bitcondition
		if (bc == newpath(t,b) || bc == bc_constant) {
			// do nothing
		} else if (newpath(t,b) == bc_constant) {
			newpath.setbitcondition(t, b, bc);
		} else if (newpath(t,b) == bc_prev || newpath(t,b) == bc_prevn) {
			if ((bc == bc_prevn || bc == bc_prev) && bc != newpath(t,b)) 
			{
				cerr << "t=" << t << ",b=" << b << ",bc=" << bc << ",path(t,b)=" << newpath(t,b) << endl;
				throw std::runtime_error("hashclash::apply_condition(...): trying to overwrite condition with new incompatible condition!");
			}
			if (bc == bc_one || bc == bc_zero) {
				int k = t;
				bitcondition bck = bc;
				while (newpath(k,b) != bc_constant && newpath(k,b) != bck) {
					if (newpath(k,b) == bc_prevn) {
						if (bck == bc_one)
							bck = bc_zero;
						else
							bck = bc_one;
					}
					if (newpath(k,b) == bc_one || newpath(k,b) == bc_zero 
						|| newpath(k,b) == bc_plus || newpath(k,b) == bc_minus)
					{
						cerr << "t=" << t << ",b=" << b << ",bc=" << bc << ",path(t,b)=" << newpath(t,b) << endl;
						throw std::runtime_error("hashclash::apply_condition(...): trying to overwrite a 01+- condition!");
					}
					--k;
				}
				newpath.setbitcondition(k,b,bck);
			}	
		} else {
			cerr << "t=" << t << ",b=" << b << ",bc=" << bc << ",path(t,b)=" << newpath(t,b) << endl;
			throw std::runtime_error("hashclash::apply_condition(...): trying to overwrite a non- .^! condition!");
		}

		// now work from begin to end to update prev, prevn's
		for (int k = newpath.tbegin()+1; k < newpath.tend(); ++k)
		{
			if (newpath(k,b) == bc_prev) {
				if (newpath(k-1,b) == bc_one || newpath(k-1,b) == bc_zero)
					newpath.setbitcondition(k,b,newpath(k-1,b));
			} else if (newpath(k,b) == bc_prevn) {
				if (newpath(k-1,b) == bc_one)
					newpath.setbitcondition(k,b,bc_zero);
				else if (newpath(k-1,b) == bc_zero)
					newpath.setbitcondition(k,b,bc_one);
			}
		} 

	}

	void enhancerot(differentialpath& path, int t, const uint32 blockdiff[])
	{
		differentialpath newpath;	

		double prot;
		bool dev32, devrc;
		vector< triple<wordconditions,wordconditions,unsigned> > solutions32, solutionsRC;

		booleanfunction* F = 0;
		if (t < 16) F = & MD5_F_data;
		else if (t < 32) F = & MD5_G_data;
		else if (t < 48) F = & MD5_H_data;
		else F = & MD5_I_data;
		uint32 dF = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			bf_outcome outcome = F->outcome(path(t,b), path(t-1,b), path(t-2,b));
			if (outcome.size()) {
				if (outcome[0] == bc_plus) 			dF += 1<<b;
				else if (outcome[0] == bc_minus)	dF -= 1<<b;
			}
		}
		uint32 dQtm3 = path[t-3].diff();
		uint32 dR = path[t+1].diff() - path[t].diff();
		uint32 dT = dQtm3 + dF + blockdiff[md5_wt[t]];
		prot = check_rotation(dR, dT, md5_rc[t], path[t], path[t+1], 1<<10);

		if (prot == 1 || prot == 0)
			return;

		dev32 = errorinrotation32(path[t], path[t+1], dT, dR, md5_rc[t], 1<<10);
		devrc = errorinrotationRC(path[t], path[t+1], dT, dR, md5_rc[t], 1<<10);
		if (dev32)
			findsolutions32(path[t], path[t+1], dT, dR, md5_rc[t], solutions32, 1<<10);
		if (devrc)
			findsolutionsRC(path[t], path[t+1], dT, dR, md5_rc[t], solutionsRC, 1<<10);
		if (0 == solutions32.size())
			solutions32.push_back( make_triple(path[t], path[t+1], unsigned(0)) );
		if (0 == solutionsRC.size())
			solutionsRC.push_back( make_triple(path[t], path[t+1], unsigned(0)) );

		double pbestsol = 0;
		unsigned addcond = 64;
		differentialpath bestpath = path;

		uint64 max32RC = solutions32.size() * solutionsRC.size();
		uint64 cnt32RC = 0;
		for (unsigned k = 0; k <= addcond+2 && cnt32RC < max32RC; ++k)
		for (unsigned i = 0; i < solutions32.size(); ++i)
		for (unsigned j = 0; j < solutionsRC.size(); ++j)
		{
			unsigned nrcond = solutions32[i].third + solutionsRC[j].third;
			if (k != nrcond) continue;
			++cnt32RC;

			try {
				newpath = path;
				for (unsigned b = md5_rc[t]; b < 32; ++b)
					if (solutions32[i].first[b] != path(t,b))
						apply_condition(newpath, t, b, solutions32[i].first[b]);
				for (unsigned b = md5_rc[t]; b < 32; ++b)
					if (solutions32[i].second[b] != path(t+1,b))
						apply_condition(newpath, t+1, b, solutions32[i].second[b]);

				for (unsigned b = 0; b < md5_rc[t]; ++b)
					if (solutionsRC[j].first[b] != path(t,b))
						apply_condition(newpath, t, b, solutionsRC[j].first[b]);
				for (unsigned b = 0; b < md5_rc[t]; ++b)
					if (solutionsRC[j].second[b] != path(t+1,b))
						apply_condition(newpath, t+1, b, solutionsRC[j].second[b]);
			} catch (...) {
				cerr << "t=" << t << endl;
				cerr << "path[t] extra conditions:" << endl;
				cerr << path[t] << endl;
				cerr << solutions32[i].first << endl;
				cerr << solutionsRC[j].first << endl;
				cerr << "path[t+1] extra conditions:" << endl;
				cerr << path[t+1] << endl;
				cerr << solutions32[i].second << endl;
				cerr << solutionsRC[j].second << endl;
				cerr << endl;
				show_path(path, blockdiff);
				cerr << endl;
				show_path(newpath, blockdiff);
				cerr << endl;
				throw;
			}
			double pnewpath = test_path(newpath, blockdiff) / double(1<<nrcond);

			if (pnewpath > pbestsol) {
				pbestsol = pnewpath;
				bestpath = newpath;
				addcond = nrcond;
			}
		}
		path = bestpath;		
	}



	void enhancepath(differentialpath& path, const uint32 blockdiff[])
	{
		for (int t = 2; t < 15; ++t)
			enhancerot(path, t, blockdiff);
	}

} // namespace hashclash
