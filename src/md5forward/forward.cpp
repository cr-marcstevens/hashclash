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

#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>

#include "main.hpp"

using namespace hashclash;
using namespace std;

void md5_forward_thread::md5_forward_differential_step(const differentialpath& path, path_container_autobalance& outpaths)
{
	const unsigned t = outpaths.t;
	const unsigned& maxcond = outpaths.maxcond;
	const unsigned maxsdrs = outpaths.maxsdrs;
	const unsigned maxweight = outpaths.maxweight;
	const unsigned minweight = outpaths.minweight;

	const uint32 m_diff_t = outpaths.m_diff[md5_wt[t]];

	booleanfunction* F = 0;
	if (t < 16) F = &MD5_F_data;
	else if (t < 32) F = &MD5_G_data;
	else if (t < 48) F = &MD5_H_data;
	else if (t < 64) F = &MD5_I_data;

	newpath = path;
	unsigned totprecond = 0;
	unsigned totcond = 0;
	for (int k = max(newpath.tbegin(), outpaths.tbegin); k < int(t)-2; ++k)
		totprecond += newpath[k].hw();
	totcond = totprecond + newpath[int(t)-2].hw() + newpath[int(t)-1].hw();
	
	unsigned minextracond = 0;
	for (unsigned b = 0; b < 32; ++b)
	{
		Qtm1b[b] = newpath(int(t)-1,b);
		Qtm2b[b] = newpath(int(t)-2,b);
		bf_outcome outconstant = F->outcome( bc_constant, Qtm1b[b], Qtm2b[b] );
		bf_outcome outplus = F->outcome( bc_plus, Qtm1b[b], Qtm2b[b] );
		bf_outcome outminus = F->outcome( bc_minus, Qtm1b[b], Qtm2b[b] );
		if (outconstant.size() > 1 && outplus.size() > 1 && outminus.size() > 1)
			++minextracond;
	}
	uint32 Qtdiff = newpath[t].diff();
	uint32 Qtm3diff = newpath[int(t)-3].diff();
	unsigned Qt_hwnaf = hwnaf(Qtdiff);
	if (totcond + Qt_hwnaf + minextracond > maxcond) return;
	wordconditions& newpathQtp1 = newpath[t+1];
	wordconditions& newpathQt = newpath[t];
	wordconditions& newpathQtm1 = newpath[int(t)-1];
	wordconditions& newpathQtm2 = newpath[int(t)-2];

	unsigned w = Qt_hwnaf+1;
	if (w < minweight) w = minweight;
        unsigned mincount = 0;
        if (minweight > 0)
                mincount = count_sdrs(Qtdiff,minweight-1);
	while (
			(w < 32) 
			&& (w+1 <= maxweight) 
			&& (totcond + w + 1 + minextracond <= maxcond) 
			&& (count_sdrs(Qtdiff, w+1)-mincount <= maxsdrs)
		)
		++w;
	table_sdrs(sdrs, Qtdiff, w);

	// we have to keep Q_0 intact if t=0
	if (t == 0) {
		sdrs.clear();
		sdrs.push_back(path[0].getsdr());
	}

	vector<sdr>::const_iterator cit = sdrs.begin(), citend = sdrs.end();
	for (; cit != citend; ++cit)
	{
		sdr sdrQt = *cit;
		unsigned hwQt = sdrQt.hw();
		if (hwQt < minweight) continue;
		if (totcond + hwQt + minextracond > maxcond) 
			continue;

		// we have to keep Q_0 intact if t=0
		if (t == 0)
			newpathQt = path[0];
		else
			newpathQt = sdrQt;

		uint32 cnt = 1;
		uint32 dF_fixed = 0;
		unsigned maxextracond = 0;
		bval.clear();
		for (unsigned b = 0; b < 32; ++b)
		{
			Qtb[b] = newpathQt.get(b);
			foutcomes[b] = F->outcome( Qtb[b], Qtm1b[b], Qtm2b[b] );
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
		if (totcond + hwQt + bval.size() > maxcond) 
			continue;
		if (outpaths.estimatefactor != 0) 
		{
			outpaths.estimate(totcond + hwQt + bval.size() + ((maxextracond+1)>>1), cnt);
			continue;
		}

		newpathQtm1 = path[int(t)-1];
		newpathQtm2 = path[int(t)-2];

		std::reverse(bval.begin(), bval.end());			
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
				const bf_conditions& newconditions = F->forwardconditions( Qtb[b], Qtm1b[b], Qtm2b[b], foutcomes[b][i] );
				newpathQt.set(b, newconditions.first);
				newpathQtm1.set(b, newconditions.second);
				newpathQtm2.set(b, newconditions.third);
				dF += foutcomes[b](i,b);
			}
			uint32 dT = dF + Qtm3diff + m_diff_t;
			newpathQtp1 = naf(Qtdiff + best_rotated_difference(dT,md5_rc[t]));
			unsigned ncond = totprecond + newpathQt.hw() + newpathQtm1.hw() + newpathQtm2.hw();
			if (outpaths.includenaf)
			{
				if (outpaths.halfnafweight)
					ncond += (newpathQtp1.hw()>>1);
				else
					ncond += newpathQtp1.hw();
			}
			if (ncond <= maxcond)
				outpaths.push_back(newpath, ncond);
		} // for cnt
	} // for sdrs	
}
