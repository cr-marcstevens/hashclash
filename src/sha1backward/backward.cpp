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
#include <string>
#include <map>

#include <boost/lexical_cast.hpp>

#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1detail.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>

#include "main.hpp"

using namespace hashclash;
using namespace std;

void sha1_backward_differential_step(const sha1differentialpath& path, path_container_autobalance& outpaths)
{
	const int t = int(outpaths.t);
	const unsigned& maxcond = outpaths.maxcond;
	const unsigned& maxsdrs = outpaths.maxsdrs;
	const unsigned& maxweight = outpaths.maxweight;
	const unsigned& minweight = outpaths.minweight;

	booleanfunction* F = 0;
	if (t < 20) F = &SHA1_F1_data; // this implementation only works properly for sha1_f1
	else if (t < 40) F = &SHA1_F2_data;
	else if (t < 60) F = &SHA1_F3_data;
	else if (t < 80) F = &SHA1_F4_data;

	static sha1differentialpath newpath;
	static vector<sdr> sdrs;
	static bitcondition Qtm1b[32], Qtm2b[32], Qtm3b[32];
	static vector<unsigned> bval;
	static uint32 fdiv[32];
	static bf_outcome foutcomes[32];	
//	static std::vector<std::pair<uint32,double> > rotateddiff;
	static vector<sdr> deltam;
	static vector<unsigned> Qtm1prev, Qtm1prevn;

	newpath = path;
	newpath[t-4].clear();

	unsigned totprecond = 0;
	unsigned totcond = 0;
	for (int k = t; k <= outpaths.tend && k <= 80; ++k)
		if (k < newpath.tend())
			totprecond += newpath[k].hw();
	totcond = totprecond + newpath[t-1].hw() + newpath[t-2].hw();	
	
	unsigned minextracond = 0;
	Qtm1prev.clear();
	Qtm1prevn.clear();
	bf_outcome outconstant, outplus, outminus;
	for (unsigned b = 0; b < 32; ++b)
	{
		Qtm1b[b] = newpath(t-1,b);
		if (Qtm1b[b] == bc_prev) {
			Qtm1prev.push_back(b);
			Qtm1b[b] = bc_constant;
			newpath[t-1].set(b, bc_constant);
		} else if (Qtm1b[b] == bc_prevn) {
			Qtm1prevn.push_back(b);
			Qtm1b[b] = bc_constant;
			newpath.setbitcondition(t-1,b, bc_constant);
		}
		Qtm2b[b] = newpath(t-2,(b+2)&31);
		if (b < 31) {
			outconstant = F->outcome( Qtm1b[b], Qtm2b[b], bc_constant );
			outplus = F->outcome( Qtm1b[b], Qtm2b[b], bc_plus );
			outminus = F->outcome( Qtm1b[b], Qtm2b[b], bc_minus );
		} else {
			outconstant = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], bc_constant );
			outplus = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], bc_plus );
			outminus = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], bc_minus );
		}
		if (outconstant.size() > 1 && outplus.size() > 1 && outminus.size() > 1)
			++minextracond;
	}
	sdr sdrQtm3 = newpath[t-3].getsdr();
	uint32 Qtm3_hwnaf = hwnaf(sdrQtm3.adddiff());
	if (totcond + Qtm3_hwnaf + minextracond > maxcond) return;

	/* table deltam */
	sdr mmask = outpaths.m_mask[t];
	deltam.clear();
	deltam.reserve(1<<mmask.hw());
	uint32 addmask = (~mmask.mask)+1; 
	uint32 andmask = mmask.mask & 0x7FFFFFFF;
	mmask.sign = 0;
	do {
		mmask.sign += addmask; mmask.sign &= andmask;
		deltam.push_back(mmask);
	} while (mmask.sign != 0 && !outpaths.onemessagediff);
#if 1
        unsigned minhw = 32;
        for (unsigned i = 0; i < deltam.size(); ++i)
                if (hwnaf(deltam[i].adddiff()) < minhw)
                        minhw = hwnaf(deltam[i].adddiff());
        unsigned ii = 0;
        while (ii < deltam.size()) {
                if (hwnaf(deltam[ii].adddiff()) > minhw) {
                        swap(deltam[ii], deltam.back());
                        deltam.pop_back();
                } else
                        ++ii;
        }
#endif

	/* table sdrs for delta Qtm1 */
	unsigned w = Qtm3_hwnaf+1;	
	if (w < minweight) w = minweight;
        unsigned mincount = 0;
        if (minweight > 0)
                mincount = count_sdrs(sdrQtm3,minweight-1,30);
	while (
			(w < 32) 
			&& (w+1 <= maxweight) 
			&& (totcond + w + 1 + minextracond <= maxcond) 
			&& (count_sdrs(sdrQtm3, w+1,30)-mincount <= maxsdrs)
		)
		++w;
	table_sdrs(sdrs, sdrQtm3, w, 30);

	wordconditions pathQtm1 = newpath[t-1];
	wordconditions& newpathQtm1 = newpath[t-1];
	wordconditions& newpathQtm2 = newpath[t-2];
	wordconditions& newpathQtm3 = newpath[t-3];
	wordconditions& newpathQtm4 = newpath[t-4];

	uint32 dQtm4pc = newpath[t+1].diff() - newpath[t].getsdr().rotate_left(5).adddiff();

	vector<sdr>::const_iterator cit = sdrs.begin(), citend = sdrs.end();
	for (; cit != citend; ++cit)
	{
		sdrQtm3 = *cit;
		unsigned hwQtm3 = sdrQtm3.hw();
		if (hwQtm3 < minweight) continue;
		if (totcond + hwQtm3 + minextracond > maxcond) 
			continue;
		newpathQtm3 = sdrQtm3;

		uint32 cnt = 1;
		uint32 dF_fixed = 0;
		unsigned maxextracond = 0;
		bval.clear();
		for (unsigned b = 0; b < 32; ++b)
		{
			Qtm3b[b] = newpathQtm3.get((b+2)&31);
			if (b < 31)
				foutcomes[b] = F->outcome( Qtm1b[b], Qtm2b[b], Qtm3b[b] );
			else
				foutcomes[b] = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], Qtm3b[b] );
			unsigned fsize = foutcomes[b].size();
			if (fsize > 1)
			{
				fdiv[b] = cnt;
				if (fsize == 2) cnt <<= 1;
				else if (fsize == 3) { cnt += cnt<<1; ++ maxextracond; }
				bval.push_back(b);
			} else
				dF_fixed += foutcomes[b](0,b);
		}
		if (totcond + hwQtm3 + bval.size() > maxcond) 
			continue;
		if (outpaths.estimatefactor != 0) 
		{
			outpaths.estimate(totcond + hwQtm3 + bval.size() + ((maxextracond+1)>>1), cnt*deltam.size());
			continue;
		}

		newpathQtm2 = path[t-2];

		std::reverse(bval.begin(), bval.end());
		bf_conditions newconditions;
		for (uint32 k = 0; k < cnt; ++k)
		{
			uint32 m = k;
			uint32 dF = dF_fixed;
			newpathQtm1 = pathQtm1;
			for (unsigned j = 0; j < bval.size(); ++j)
			{
				const unsigned b = bval[j];
				unsigned i = 0;
				while (m >= fdiv[b])
				{
					m -= fdiv[b];
					++i;
				}
				if (b < 31)
					newconditions = F->backwardconditions( Qtm1b[b], Qtm2b[b], Qtm3b[b], foutcomes[b][i] );
				else {
					newconditions = msb_bf_backwardconditions(*F, Qtm1b[b], Qtm2b[b], Qtm3b[b], foutcomes[b][i] );
					if (msb_bf_outcome(*F, newconditions).size() > 1 || msb_bf_outcome(*F, newconditions)[0] != foutcomes[b][i]) {
						cout << endl << "[" << Qtm1b[b] << Qtm2b[b] << Qtm3b[b] <<"](" << foutcomes[b][i] <<")=>[" << newconditions.first << newconditions.second << newconditions.third << "]" << endl;
						bf_outcome tmp = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], Qtm3b[b]);
						for (unsigned i = 0; i < tmp.size(); ++i)
							cout << tmp[i] << flush;
						cout << endl;
						tmp = F->outcome(Qtm1b[b], Qtm2b[b], Qtm3b[b]);
						for (unsigned i = 0; i < tmp.size(); ++i)
							cout << tmp[i] << flush;
						cout << endl;
						newconditions = F->backwardconditions( Qtm1b[b], Qtm2b[b], Qtm3b[b], bc_constant );
						cout << "(.)=>[" << newconditions.first << newconditions.second << newconditions.third << "]" << endl;
						throw;
					}
				}
				newpathQtm1.set(b, newconditions.first);
				newpathQtm2.set((b+2)&31, newconditions.second);
				newpathQtm3.set((b+2)&31, newconditions.third);

				dF += foutcomes[b](i,b);
			}

			bool contradiction = false;
			for (unsigned i = 0; i < Qtm1prev.size(); ++i) {
				const unsigned b = Qtm1prev[i];
				const bitcondition bcQtm1b = newpathQtm1[b];
				switch (bcQtm1b) {
				case bc_constant:
					newpathQtm1.set(b, bc_prev);
					break;
				case bc_one:
				case bc_zero:
					if (newpathQtm2[b] == bc_prev) {
						newpathQtm3.set(b, bcQtm1b);
					} else if (newpathQtm2[b] == bc_prevn) {
						if (bcQtm1b == bc_one)
							newpathQtm3.set(b, bc_zero);
						else
							newpathQtm3.set(b, bc_one);
					} else if (newpathQtm2[b] != bc_constant && newpathQtm2[b] != bcQtm1b) {
						//cerr << "[^" << newpathQtm1[b] << newpathQtm2[b] << newpathQtm3[b] << "]" << flush;
						contradiction = true;
					}
					newpathQtm2.set(b, bcQtm1b);
					break;
				default:
					cerr << "[^" << newpathQtm1[b] << newpathQtm2[b] << newpathQtm3[b] << "]" << flush;
					contradiction = true;
				}
			}
			for (unsigned i = 0; i < Qtm1prevn.size(); ++i) {
				const unsigned b = Qtm1prevn[i];
				const bitcondition bcQtm1b = newpathQtm1[b];
				switch (bcQtm1b) {
				case bc_constant:
					newpathQtm1.set(b, bc_prevn);
					break;
				case bc_one:
				case bc_zero:
					if (newpathQtm2[b] == bc_prev) {
						if (bcQtm1b == bc_one)
							newpathQtm3.set(b, bc_zero);
						else
							newpathQtm3.set(b, bc_one);
					} else if (newpathQtm2[b] == bc_prevn) {
						newpathQtm3.set(b, bcQtm1b);
					} else if (newpathQtm2[b] != bc_constant && newpathQtm2[b] == bcQtm1b) {
						//cerr << "[^" << newpathQtm1[b] << newpathQtm2[b] << newpathQtm3[b] << "]" << flush;
						contradiction = true;
					}
					if (bcQtm1b == bc_one)
						newpathQtm2.set(b, bc_zero);
					else
						newpathQtm2.set(b, bc_one);
					break;
				default:
					cerr << "[^" << newpathQtm1[b] << newpathQtm2[b] << newpathQtm3[b] << "]" << flush;
					contradiction = true;
				}
			}
			if (contradiction)
				continue;

			unsigned ncond = totprecond + newpathQtm1.hw() + newpathQtm2.hw() + newpathQtm3.hw();
			if (ncond > maxcond) return;
			for (unsigned i = 0; i < deltam.size(); ++i)
			{
				newpath.getme(t) = deltam[i];
				uint32 dQtm4 = dQtm4pc - dF - deltam[i].adddiff();
				newpathQtm4 = naf(dQtm4).rotate_right(30);
				unsigned ncond2 = ncond;
				if (outpaths.includenaf) {
					if (outpaths.halfnafweight)
						ncond2 += (newpathQtm4.hw()>>1);
					else
						ncond2 += newpathQtm4.hw();
				}
				if (ncond2 <= maxcond)
					outpaths.push_back(newpath, ncond2);
			}
		} // for cnt
	} // for sdrs	
}
