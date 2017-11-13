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

void sha1_forward_differential_step(const sha1differentialpath& path, path_container_autobalance& outpaths)
{
	const int t = int(outpaths.t);
	const unsigned& maxcond = outpaths.maxcond;
	const unsigned& maxsdrs = outpaths.maxsdrs;
	const unsigned& maxweight = outpaths.maxweight;
	const unsigned& minweight = outpaths.minweight;

	booleanfunction* F = 0;
	if (t < 20) F = &SHA1_F1_data;
	else if (t < 40) F = &SHA1_F2_data;
	else if (t < 60) F = &SHA1_F3_data;
	else if (t < 80) F = &SHA1_F4_data;

	static sha1differentialpath newpath;
	static vector<sdr> sdrs;
	static bitcondition Qtm1b[32], Qtm2b[32], Qtm3b[32];
	static vector<unsigned> bval;
	static uint32 fdiv[32];
	static bf_outcome foutcomes[32];
	static vector<sdr> deltam;

	newpath = path;
	unsigned totprecond = 0;
	unsigned totcond = 0;
	for (int k = outpaths.tbegin; k < t-3; ++k)
		if (k >= newpath.tbegin())
			totprecond += newpath[k].hw();
	totcond = totprecond + newpath[t-2].hw() + newpath[t-3].hw();
	
	unsigned minextracond = 0;
	bf_outcome outconstant, outplus, outminus;
	for (unsigned b = 0; b < 32; ++b)
	{
		Qtm2b[b] = newpath(t-2,(b+2)&31);
		Qtm3b[b] = newpath(t-3,(b+2)&31);
		if (b < 31) {
			outconstant = F->outcome( bc_constant, Qtm2b[b], Qtm3b[b] );
			outplus = F->outcome( bc_plus, Qtm2b[b], Qtm3b[b] );
			outminus = F->outcome( bc_minus, Qtm2b[b], Qtm3b[b] );
		} else {
			outconstant = msb_bf_outcome(*F, bc_constant, Qtm2b[b], Qtm3b[b] );
			outplus = msb_bf_outcome(*F, bc_plus, Qtm2b[b], Qtm3b[b] );
			outminus = msb_bf_outcome(*F, bc_minus, Qtm2b[b], Qtm3b[b] );
		}
		if (outconstant.size() > 1 && outplus.size() > 1 && outminus.size() > 1)
			++minextracond;
	}
	unsigned Qtm1_hwnaf = hwnaf(newpath[t-1].diff());
	unsigned Qt_hwnaf = hwnaf(newpath[t].diff());
	if (totcond + Qt_hwnaf + Qtm1_hwnaf + minextracond > maxcond) return;

	/* table deltam */
	sdr mmask = outpaths.m_mask[t];
	deltam.clear();
	deltam.reserve(1<<mmask.hw());
	uint32 addmask = (~mmask.mask)+1; 
	uint32 andmask = mmask.mask & ~0x80000000;
	mmask.sign = 0;
	do {
		mmask.sign += addmask; mmask.sign &= andmask;
		deltam.push_back(mmask);
	} while (mmask.sign != 0 && !outpaths.onemessagediff);

	/* table sdrs for delta Qtm1 */
	unsigned w = Qtm1_hwnaf+1;
	sdr sdrQtm1 = newpath[t-1].getsdr();
	if (w < minweight) w = minweight;
        unsigned mincount = 0;
        if (minweight > 0)
                mincount = count_sdrs(sdrQtm1,minweight-1,5);
	while (
			(w < 32) 
			&& (w+1 <= maxweight) 
			&& (totcond + w + 1 + minextracond <= maxcond) 
			&& (count_sdrs(sdrQtm1, w+1,5)-mincount <= maxsdrs)
		)
		++w;
	table_sdrs(sdrs, sdrQtm1, w, 5);
	// we have to keep Q_t-1 intact if t=0 or 1
	if (t == 0 || t == 1) {
		sdrs.clear();
		sdrs.push_back(path[t-1].getsdr());
	}

	wordconditions& newpathQtp1 = newpath[t+1];
	wordconditions& newpathQt = newpath[t];
	wordconditions& newpathQtm1 = newpath[t-1];
	wordconditions& newpathQtm2 = newpath[t-2];
	wordconditions& newpathQtm3 = newpath[t-3];

	uint32 dQtp1pc = newpath[t].getsdr().rotate_left(5).adddiff() + newpath[t-4].getsdr().rotate_left(30).adddiff();
	vector<sdr>::const_iterator cit = sdrs.begin(), citend = sdrs.end();
	for (; cit != citend; ++cit)
	{
		sdrQtm1 = *cit;
		unsigned hwQtm1 = sdrQtm1.hw();
		if (hwQtm1 < minweight) continue;
		if (hwQtm1 > maxweight) continue;
		if (totcond + Qt_hwnaf + hwQtm1 + minextracond > maxcond) 
			continue;

		// we have to keep Q_t-1 intact if t=0 or 1
		if (t == 0 || t == 1)
			newpathQtm1 = path[t-1];
		else
			newpathQtm1 = sdrQtm1;

		uint32 cnt = 1;
		uint32 dF_fixed = 0;
		unsigned maxextracond = 0;
		bval.clear();
		for (unsigned b = 0; b < 32; ++b)
		{
			Qtm1b[b] = newpathQtm1.get(b);
			if (b == 31)
				foutcomes[b] = msb_bf_outcome(*F, Qtm1b[b], Qtm2b[b], Qtm3b[b]);
			else
				foutcomes[b] = F->outcome( Qtm1b[b], Qtm2b[b], Qtm3b[b] );
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
		if (totcond + Qt_hwnaf + hwQtm1 + bval.size() > maxcond) 
			continue;
		if (outpaths.estimatefactor != 0) 
		{
			outpaths.estimate(totcond + Qt_hwnaf + hwQtm1 + bval.size() + ((maxextracond+1)>>1), cnt*deltam.size());
			continue;
		}

		newpathQtm2 = path[int(t)-2];
		newpathQtm3 = path[int(t)-3];

		std::reverse(bval.begin(), bval.end());
		bf_conditions newconditions;
		for (uint32 k = 0; k < cnt; ++k)
		{
			uint32 m = k;
			uint32 dF = dF_fixed;
			for (unsigned j = 0; j < bval.size(); ++j)
			{
				const unsigned b = bval[j];
				unsigned i = 0;
				while (m >= fdiv[b])
				{
					m -= fdiv[b];
					++i;
				}
				if (b == 31)
					newconditions = msb_bf_forwardconditions(*F, Qtm1b[b], Qtm2b[b], Qtm3b[b], foutcomes[b][i]);
				else
					newconditions = F->forwardconditions( Qtm1b[b], Qtm2b[b], Qtm3b[b], foutcomes[b][i] );
				newpathQtm1.set(b, newconditions.first);
				newpathQtm2.set((b+2)&31, newconditions.second);
				newpathQtm3.set((b+2)&31, newconditions.third);
				dF += foutcomes[b](i,b);
			}
			unsigned ncond = totprecond + newpathQtm3.hw() + newpathQtm2.hw() + newpathQtm1.hw() + Qt_hwnaf;
			if (ncond > maxcond) continue;
			for (unsigned i = 0; i < deltam.size(); ++i)
			{
				newpath.getme(t) = deltam[i];
				uint32 dQtp1 = dQtp1pc + dF + deltam[i].adddiff();
				newpathQtp1 = naf(dQtp1);
				unsigned ncond2 = ncond;
				if (outpaths.includenaf) {
					if (outpaths.halfnafweight)
						ncond2 += (newpathQtp1.hw()>>1);
					else
						ncond2 += newpathQtp1.hw();
				}
				if (ncond2 <= maxcond)
					outpaths.push_back(newpath, ncond2);
			}
		} // for cnt
	} // for sdrs	

}
