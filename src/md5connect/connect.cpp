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

#include "main.hpp"

#include <hashclash/timer.hpp>
#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>



using namespace hashclash;
using namespace std;

booleanfunction* md5bf(unsigned t)
{
	if (t < 16) return &MD5_F_data;
	else if (t < 32) return &MD5_G_data;
	else if (t < 48) return &MD5_H_data;
	return &MD5_I_data;
}

//unsigned t;
//booleanfunction* Ft;
//booleanfunction* Ftp1;
//booleanfunction* Ftp2;
//booleanfunction* Ftp3;
//vector<uint32> dFt, dFtp1, dFtp2, dFtp3;
//uint32 dmt, dmtp1, dmtp2, dmtp3;
//uint32 dQtp1, dQtp2, dQtp3, dQtp4;
//differentialpath newpath;

//bool lastdFp1, lastdFp2, lastdFp3;
void md5_connect_thread::connectbits(const connect_bitdata& in, vector<connect_bitdata>& out,
				 unsigned b, const differentialpath& lower, const differentialpath& upper,
				 vector<byteconditions>* newconds/* = 0*/
				 )
{
//	static connect_bitdata result;
//	static bf_outcome bfo0, bfo1, bfo2, bfo3;
//	static bf_conditions bfc0, bfc1, bfc2, bfc3;
//	static bitcondition Qt, Qtp1;
//	static byteconditions newcond;
	bitcondition Qtp2 = upper[t+2][b];
	bitcondition Qtp3 = upper[t+3][b];
	bitcondition Qtm1 = lower[t-1][b];
	bitcondition Qtm2 = lower[t-2][b];
	unsigned dQtbit = (in.dQt>>b)&1;
	unsigned dQtp1bit = (in.dQtp1>>b)&1;

	/**** carry dQt ****/
	for (unsigned dqt = 0; dqt <= dQtbit; ++dqt)
	{
		if (dqt == 0) {
			if (dQtbit == 0) {
				Qt = bc_constant;
				result.dQt = in.dQt;
			} else {
				Qt = bc_minus;
				result.dQt = in.dQt + (1<<b);
			}
		} else {
			Qt = bc_plus;
			result.dQt = in.dQt - (1<<b);
		}

		/**** bitconditions Ft ****/
		bfo0 = Ft->outcome(Qt, Qtm1, Qtm2);
		for (unsigned bfo0index = 0; bfo0index < bfo0.size(); ++bfo0index)
		{			
			result.dFt = in.dFt - bfo0(bfo0index, b);
			if (result.dFt & (1<<b)) continue;			
			bfc0 = Ft->forwardconditions(Qt, Qtm1, Qtm2, bfo0[bfo0index]);

			/**** carry dQtp1 ****/
			for (unsigned dqtp1 = 0; dqtp1 <= dQtp1bit; ++dqtp1)
			{
				if (dqtp1 == 0) {
					if (dQtp1bit == 0) {
						Qtp1 = bc_constant;
						result.dQtp1 = in.dQtp1;
					} else {
						Qtp1 = bc_minus;
						result.dQtp1 = in.dQtp1 + (1<<b);
					}
				} else {
					Qtp1 = bc_plus;
					result.dQtp1 = in.dQtp1 - (1<<b);
				}

				/**** bitconditions Ftp1 ****/				
				bfo1 = Ftp1->outcome(Qtp1, bfc0.first, bfc0.second);
				for (unsigned bfo1index = 0; bfo1index < bfo1.size(); ++bfo1index)
				{					
					result.dFtp1 = in.dFtp1 - bfo1(bfo1index, b);
					if (result.dFtp1 & (1<<b)) continue;					
					bfc1 = Ftp1->forwardconditions(Qtp1, bfc0.first, bfc0.second, bfo1[bfo1index]);

					/**** bitconditions Ftp2 ****/
					bfo2 = Ftp2->outcome(Qtp2, bfc1.first, bfc1.second);
					for (unsigned bfo2index = 0; bfo2index < bfo2.size(); ++bfo2index)
					{						
						result.dFtp2 = in.dFtp2 - bfo2(bfo2index, b);
						if (result.dFtp2 & (1<<b)) continue;						
						bfc2 = Ftp2->forwardconditions(Qtp2, bfc1.first, bfc1.second, bfo2[bfo2index]);

						/**** bitconditions Ftp3 ****/
						bfo3 = Ftp3->outcome(Qtp3, bfc2.first, bfc2.second);
						for (unsigned bfo3index = 0; bfo3index < bfo3.size(); ++bfo3index)
						{						
							result.dFtp3 = in.dFtp3 - bfo3(bfo3index, b);
							if (result.dFtp3 & (1<<b)) continue;
							
							out.push_back(result);							
							if (newconds)
							{
								bfc3 = Ftp3->forwardconditions(Qtp3, bfc2.first, bfc2.second, bfo3[bfo3index]);
								newcond.set(
									bfc0.third, bfc1.third, bfc2.third, 
									bfc3.third, bfc3.second, bfc3.first);
								newconds->push_back(newcond);
							}
						} // dFtp3
					} // dFtp2
				} // dFtp1
			} // dQtp1
		} // dFt
	} // dQt
}

void md5_connect_thread::connectbits2(const connect_bitdata& in,
				 unsigned b, const differentialpath& lower, const differentialpath& upper)
{
//	static bf_outcome bfo0, bfo1, bfo2, bfo3;
//	static bf_conditions bfc0, bfc1, bfc2, bfc3;
	bitcondition Qt, Qtp1;
	bitcondition Qtp2 = upper[t+2][b];
	bitcondition Qtp3 = upper[t+3][b];
	bitcondition Qtm1 = lower[t-1][b];
	bitcondition Qtm2 = lower[t-2][b];
	unsigned dQtbit = (in.dQt>>b)&1;
	unsigned dQtp1bit = (in.dQtp1>>b)&1;
	/**** carry dQt ****/
	for (unsigned dqt = 0; dqt <= dQtbit; ++dqt)
	{
		if (dqt == 0) {
			if (dQtbit == 0) {
				Qt = bc_constant;
			} else {
				Qt = bc_minus;
			}
		} else {
			Qt = bc_plus;
		}

		/**** bitconditions Ft ****/
		bfo0 = Ft->outcome(Qt, Qtm1, Qtm2);
		for (unsigned bfo0index = 0; bfo0index < bfo0.size(); ++bfo0index)
		{			
			if ((in.dFt - bfo0(bfo0index, b)) & (1<<b)) continue;			
			bfc0 = Ft->forwardconditions(Qt, Qtm1, Qtm2, bfo0[bfo0index]);

			/**** carry dQtp1 ****/
			lastdFp1 = true;
			for (unsigned dqtp1 = 0; dqtp1 <= dQtp1bit; ++dqtp1)
			{
				if (dqtp1 == 0) {
					if (dQtp1bit == 0) {
						Qtp1 = bc_constant;
					} else {
						Qtp1 = bc_minus;
					}
				} else {
					Qtp1 = bc_plus;
				}

				/**** bitconditions Ftp1 ****/				
				bfo1 = Ftp1->outcome(Qtp1, bfc0.first, bfc0.second);
				for (unsigned bfo1index = 0; bfo1index < bfo1.size(); ++bfo1index)
				{					
					if ((in.dFtp1 - bfo1(bfo1index, b)) & (1<<b)) continue;					
					bfc1 = Ftp1->forwardconditions(Qtp1, bfc0.first, bfc0.second, bfo1[bfo1index]);

					/**** bitconditions Ftp2 ****/
					lastdFp2 = true;
					bfo2 = Ftp2->outcome(Qtp2, bfc1.first, bfc1.second);
					for (unsigned bfo2index = 0; bfo2index < bfo2.size(); ++bfo2index)
					{						
						if ((in.dFtp2 - bfo2(bfo2index, b)) & (1<<b)) continue;						

						lastdFp3 = true;
					} // dFtp2
				} // dFtp1
			} // dQtp1
		} // dFt
	} // dQt
}

unsigned md5_connect_thread::md5_connect_bits(const vector<differentialpath>& lowers, unsigned index,
						  const differentialpath& upper, path_container& container)
{
	const differentialpath& lower = lowers[index];
//	static connect_bitdata startbit0;
//	static vector<connect_bitdata> bitdataresults[33];


	// first verify whether lower and upper paths can be connected
	startbit0.dFt = dFt[index];
	startbit0.dFtp1 = dFtp1[index];
	startbit0.dFtp2 = dFtp2[index];
	startbit0.dFtp3 = dFtp3[index];
	startbit0.dQt = lowdQt[index]; //lower[t].diff();
	startbit0.dQtp1 = dQtp1; // upper[t+1].diff();
	
	bitdataresults[0].clear();
	bitdataresults[0].push_back(startbit0);	

	for (unsigned b = 0; b < 32; ++b)
	{
		bitdataresults[b+1].clear();
		for (unsigned k = 0; k < bitdataresults[b].size(); ++k)
			connectbits(bitdataresults[b][k], bitdataresults[b+1], b, lower, upper);
		if (bitdataresults[b+1].size() == 0)
		{
			// they cannot be connected
			// check whether dF[t+1][b], dF[t+2][b] or dF[t+3][b] were involved or not
			lastdFp1 = false;
			lastdFp2 = false;
			lastdFp3 = false;
			for (unsigned k = 0; k < bitdataresults[b].size(); ++k)
				connectbits2(bitdataresults[b][k], b, lower, upper);
			return b;
		}

		// remove duplicates
		sort(bitdataresults[b+1].begin(), bitdataresults[b+1].end());
		bitdataresults[b+1].erase( unique(bitdataresults[b+1].begin(), bitdataresults[b+1].end())
									, bitdataresults[b+1].end());
	}
	// they can be connected
	// we rebuild the list and now keeping track
	// which bitconds lead which connect_bitdata to which connect_bitdata
//	static vector<connect_bitdata> bitdatastart[32];
//	static vector<connect_bitdata> bitdataend[32];
//	static vector<byteconditions>  bitdatanewcond[32];
	for (unsigned b = 0; b < 32; ++b)
	{
		bitdatastart[b].clear();
		bitdataend[b].clear();
		bitdatanewcond[b].clear();

		for (unsigned k = 0; k < bitdataresults[b].size(); ++k)
		{
			connectbits(bitdataresults[b][k], bitdataend[b], b, lower, upper, &bitdatanewcond[b]);
			bitdatastart[b].resize(bitdataend[b].size(), bitdataresults[b][k]);
		}
		if (bitdatastart[b].size() != bitdataend[b].size()
			|| bitdataend[b].size() != bitdatanewcond[b].size())
			throw;
	}

	//static differentialpath newpath2;
	newpath2 = upper;
	for (int k = lower.tbegin(); k < lower.tend(); ++k)
		newpath2[k] = lower[k];

#if 1
	vector< vector<unsigned> > tunnelstrength(32);
	vector< vector<unsigned> > colbitconditions(32);
	vector< vector<unsigned> > cumtunnelstrength, cumbitconditions;
	unsigned maxtunnelstrength = 0, minconditions = 9999;
	if (container.Qcondstart > t+3) {
		if (tmppath.path.size() == 0) {
			tmppath.get(-3); tmppath.get(16);
			for (int i = -3; i <= 16; ++i)
				for (unsigned b = 0; b < 32; ++b)
					tmppath.setbitcondition(i,b,bc_one);
		}	
		for (unsigned b = 0; b < 32; ++b) {
			for (int i = newpath2.tbegin(); i < newpath2.tend(); ++i)
				tmppath.setbitcondition(i, b, newpath2(i,b));
			unsigned bc = 0;
			for (int i = newpath2.tbegin(); i < newpath2.tend(); ++i)
				if ((i < t-2 || i > t+3) && newpath2(i,b) != bc_constant)
					++bc;
			for (unsigned i = 0; i < bitdatanewcond[b].size(); ++i) {
				tmppath.setbitcondition(t-2,b,bitdatanewcond[b][i][0]);
				tmppath.setbitcondition(t-1,b,bitdatanewcond[b][i][1]);
				tmppath.setbitcondition(t,b,bitdatanewcond[b][i][2]);
				tmppath.setbitcondition(t+1,b,bitdatanewcond[b][i][3]);
				tmppath.setbitcondition(t+2,b,bitdatanewcond[b][i][4]);
				tmppath.setbitcondition(t+3,b,bitdatanewcond[b][i][5]);
				tunnelstrength[b].push_back(totaltunnelstrength(tmppath));
				unsigned bcc = bc;
				if (bitdatanewcond[b][i][0] != bc_constant) ++bcc;
				if (bitdatanewcond[b][i][1] != bc_constant) ++bcc;
				if (bitdatanewcond[b][i][2] != bc_constant) ++bcc;
				if (bitdatanewcond[b][i][3] != bc_constant) ++bcc;
				if (bitdatanewcond[b][i][4] != bc_constant) ++bcc;
				if (bitdatanewcond[b][i][5] != bc_constant) ++bcc;
				colbitconditions[b].push_back(bcc);
			}
			for (int i = newpath.tbegin(); i < newpath.tend(); ++i)
				tmppath.setbitcondition(i, b, bc_one);
		}
		cumtunnelstrength = tunnelstrength;
		cumbitconditions = colbitconditions;
		for (int b = 0; b < 31; ++b) {
			for (unsigned i = 0; i < bitdatanewcond[b].size(); ++i) {
				unsigned cts = cumtunnelstrength[b][i];
				pair< vector<connect_bitdata>::const_iterator, vector<connect_bitdata>::const_iterator >
					bounds = equal_range(bitdatastart[b+1].begin(), bitdatastart[b+1].end(), bitdataend[b][i]);
				for (; bounds.first != bounds.second; ++bounds.first) {
					unsigned j = bounds.first - bitdatastart[b+1].begin();
					if (tunnelstrength[b+1][j] + cts >= cumtunnelstrength[b+1][j]) {
						if (tunnelstrength[b+1][j] + cts > cumtunnelstrength[b+1][j] || colbitconditions[b+1][j] + cumbitconditions[b][i] < cumbitconditions[b+1][j])
							cumbitconditions[b+1][j] = colbitconditions[b+1][j] + cumbitconditions[b][i];
						cumtunnelstrength[b+1][j] = tunnelstrength[b+1][j] + cts;
					}
				}
			}
		}
		for (unsigned i = 0; i < cumtunnelstrength[31].size(); ++i)
			if (cumtunnelstrength[31][i] >= maxtunnelstrength) {
				if (cumtunnelstrength[31][i] > maxtunnelstrength || cumbitconditions[31][i] < minconditions)
					minconditions = cumbitconditions[31][i];
				maxtunnelstrength = cumtunnelstrength[31][i];
			}
		int maxcompl = maxtunnelstrength;
		for (int k = newpath2.tbegin(); k < newpath2.tend() && k < 64; ++k)
			if (k >= container.Qcondstart)
				maxcompl -= int(newpath2[k].hw());
		if (maxcompl < container.bestmaxcomp) return 32;
		if (maxtunnelstrength < container.bestmaxtunnel) return 32;
		if (minconditions > container.bestpathcond) return 32;
	}
#endif

	//static unsigned bindex[32];
	unsigned bit = 31;
	bindex[bit] = 0;
	vector<unsigned> tsv(33,0), bcv(33,0);
	tsv[bit] = tunnelstrength[31][0];
	bcv[bit] = colbitconditions[31][0];
	while (bit <= 31)
	{
		if (bindex[bit] < bitdataend[bit].size())
		{	
#if 1
			unsigned ts = tsv[bit] = tsv[bit+1] + tunnelstrength[bit][bindex[bit]];
			unsigned bc = bcv[bit] = bcv[bit+1] + colbitconditions[bit][bindex[bit]];
			// ts == tsv[bit+1] after --bit;
#endif
			--bit;
			if (bit != 0) {
				bindex[bit] = bitdataend[bit].size();
				for (unsigned k = 0; k < bitdataend[bit].size(); ++k)
					if (bitdataend[bit][k] == bitdatastart[bit+1][bindex[bit+1]])
					{
#if 1
						if (cumtunnelstrength[bit][k] + ts < container.bestmaxtunnel) continue;
						if (cumtunnelstrength[bit][k] + ts == container.bestmaxtunnel && cumbitconditions[bit][k] + bc > container.bestpathcond) continue;
#endif
						bindex[bit] = k;
						break;
					}
			} else {
				for (unsigned k = 0; k < bitdataend[bit].size(); ++k)
					if (bitdataend[bit][k] == bitdatastart[bit+1][bindex[bit+1]])
					{ // we have a full path !!!
						bindex[0] = k;
#if 1
						if (cumtunnelstrength[bit][k] + ts < container.bestmaxtunnel) continue;
						if (cumtunnelstrength[bit][k] + ts == container.bestmaxtunnel && cumbitconditions[bit][k] + bc > container.bestpathcond) continue;
#endif						
						for (unsigned b = 0; b < 32; ++b)
						{
							newpath2.setbitcondition(t-2,b,bitdatanewcond[b][bindex[b]][0]);
							newpath2.setbitcondition(t-1,b,bitdatanewcond[b][bindex[b]][1]);
							newpath2.setbitcondition(t+0,b,bitdatanewcond[b][bindex[b]][2]);
							newpath2.setbitcondition(t+1,b,bitdatanewcond[b][bindex[b]][3]);
							newpath2.setbitcondition(t+2,b,bitdatanewcond[b][bindex[b]][4]);
							newpath2.setbitcondition(t+3,b,bitdatanewcond[b][bindex[b]][5]);
						}
						container.push_back(newpath2);
					}
				bindex[bit] = bitdataend[bit].size();
			}
		} else
		{
			++bit;
			if (bit <= 31)
			{
#if 1
				unsigned ts = tsv[bit+1];
				unsigned bc = bcv[bit+1];
#endif
				++bindex[bit];
				if (bit < 31)
				{
					while (bindex[bit] < bitdataend[bit].size()
						&& (((bitdataend[bit][bindex[bit]]) != (bitdatastart[bit+1][bindex[bit+1]]))
#if 1
							|| (cumtunnelstrength[bit][bindex[bit]] + ts < container.bestmaxtunnel)
							|| (cumtunnelstrength[bit][bindex[bit]] + ts == container.bestmaxtunnel && cumbitconditions[bit][bindex[bit]] + bc > container.bestpathcond)
#endif
						))
						++bindex[bit];
				} else {
					while (bindex[bit] < bitdataend[bit].size()
#if 1
						&& ((cumtunnelstrength[bit][bindex[bit]] + ts < container.bestmaxtunnel)
							|| (cumtunnelstrength[bit][bindex[bit]] + ts == container.bestmaxtunnel && cumbitconditions[bit][bindex[bit]] + bc > container.bestpathcond)
						)
#endif
						)
						++bindex[bit];
				}
			}
		}
	}


	return 32;
}

void md5_connect_thread::md5_connect(const vector<differentialpath>& lowerpaths
				 , const differentialpath& upperpath
				 , path_container& container)
{
//	static timer sw(true);
//	static vector<unsigned> countb(33, 0);
//	static vector<unsigned> countbaborted(33, 0);
//	static vector<unsigned> countbdepth(33, 0);
//	static uint64 count, countall;
//	static vector<unsigned char> isgood;
//	static vector<int> lowerpathsmaxtunnel;
	if (lowerpathsmaxtunnel.size() == 0) {
		lowerpathsmaxtunnel.resize(lowerpaths.size());
		for (unsigned i = 0; i < lowerpaths.size(); ++i)
		{
			differentialpath tmp = lowerpaths[i];
		try {
		cleanup(tmp);
		} catch (std::exception& e) {
			cerr << "hashclash::cleanup(differentialpath&): unknown exception!:" << endl << e.what() << endl;
			show_path(tmp, container.m_diff);
		}catch (...) {
			cerr << "hashclash::cleanup(differentialpath&): unknown exception!:" << endl;
			show_path(tmp, container.m_diff);
		}

			//cleanup(tmp);
			lowerpathsmaxtunnel[i] = totaltunnelstrength(tmp);
		}
	}
        int uppercompl = 0;
        for (int k = upperpath.tbegin(); k < upperpath.tend() && k < 64; ++k)
        	if (k >= container.Qcondstart)
                	uppercompl -= int(upperpath[k].hw());
	if (container.showstats && sw.time() > 60)
	{
		cout << endl << count << "\t" << countall << "\t" 
			<< double(countall)/double(count) << endl;
		for (unsigned b = 0; b <= 32; ++b)
		{
			cout << b << ":\t" << countb[b] << "\t" << countbaborted[b] << 
				"\t" << double(countbaborted[b])/double(countb[b]) <<
				"\t" << double(countbdepth[b])/double(countb[b]) << endl;
		}
		sw.start();
	}
	t = container.t;
	Ft = md5bf(t);
	Ftp1 = md5bf(t+1);
	Ftp2 = md5bf(t+2);
	Ftp3 = md5bf(t+3);
	dmt = container.m_diff[md5_wt[t]];
	dmtp1 = container.m_diff[md5_wt[t+1]];
	dmtp2 = container.m_diff[md5_wt[t+2]];
	dmtp3 = container.m_diff[md5_wt[t+3]];

	unsigned bequal = 0;
	if (0 != newpath.path.size())
	{
		while (bequal < 32
			&& newpath[t+1][bequal] == upperpath[t+1][bequal]
			&& newpath[t+2][bequal] == upperpath[t+2][bequal]
			&& newpath[t+3][bequal] == upperpath[t+3][bequal]
			)
			++bequal;		
	}
	newpath = upperpath;

	isgood.resize(lowerpaths.size(),0);
	dFt.resize(lowerpaths.size());
	dFtp1.resize(lowerpaths.size());
	dFtp2.resize(lowerpaths.size());
	dFtp3.resize(lowerpaths.size());
	dQtp1 = upperpath[t+1].diff();
	dQtp2 = upperpath[t+2].diff();
	dQtp3 = upperpath[t+3].diff();
	dQtp4 = upperpath[t+4].diff();
	uint32 dTtp1 = best_rotated_difference(dQtp2 - dQtp1, 32-md5_rc[t+1]) - dmtp1;
	uint32 dTtp2 = best_rotated_difference(dQtp3 - dQtp2, 32-md5_rc[t+2]) - dmtp2;
	uint32 dTtp3 = best_rotated_difference(dQtp4 - dQtp3, 32-md5_rc[t+3]) - dmtp3;
	for (unsigned i = 0; i < lowerpaths.size(); ++i)
	{
		uint32 dQt = lowdQt[i];
		uint32 dTt = best_rotated_difference(dQtp1 - dQt, 32-md5_rc[t]);
		uint32 dFti = dTt - dmt - lowdQtm3[i];
		uint32 dFtp1i = dTtp1 - lowdQtm2[i];
		uint32 dFtp2i = dTtp2 - lowdQtm1[i];
		uint32 dFtp3i = dTtp3 - dQt;
		if (isgood[i] != 0)
		{
			if (isgood[i] <= bequal)
			{
				if (0 != ((dFti^dFt[i])<<(32-isgood[i]))
					|| 0 != ((dFtp1i^dFtp1[i])<<(32-isgood[i]))
					|| 0 != ((dFtp2i^dFtp2[i])<<(32-isgood[i]))
					|| 0 != ((dFtp3i^dFtp3[i])<<(32-isgood[i]))
					)
					isgood[i] = 0;
			} else
				isgood[i] = 0;
		} 		
		dFt[i] = dFti;
		dFtp1[i] = dFtp1i;
		dFtp2[i] = dFtp2i;
		dFtp3[i] = dFtp3i;
	}
	
	countall += lowerpaths.size();
	for (unsigned i = 0; i < lowerpaths.size(); ++i)
	{
		if (isgood[i]) continue;
		if (lowerpathsmaxtunnel[i]+uppercompl <= container.bestmaxcomp) continue;
		++count;
		unsigned b = md5_connect_bits(lowerpaths, i, upperpath, container);
		++countb[b];
#ifndef DONT_SKIP_LOWERPATHS
		// skip all lowerpaths that have the same characteristics on the b lower bits
		if (b < 32) {
			isgood[i] = b+1;

			unsigned j = i + 1;
			unsigned jmax = binary_search_lower_paths(lowerpaths, i, b) + 1;
			countbdepth[b] += jmax - j;

			uint32 maskt = 0;
			for (unsigned k = 0; k < b; ++k)
				maskt |= 1<<k;
			uint32 masktp1 = maskt;
			uint32 masktp2 = maskt;
			uint32 masktp3 = maskt;
			maskt |= 1<<b;
			if (lastdFp1) masktp1 |= 1<<b;
			if (lastdFp2) masktp2 |= 1<<b;
			if (lastdFp3) masktp3 |= 1<<b;

			uint32 dFti = dFt[i] & maskt;
			uint32 dFtp1i = dFtp1[i] & masktp1;
			uint32 dFtp2i = dFtp2[i] & masktp2;
			uint32 dFtp3i = dFtp3[i] & masktp3;
			for (; j < jmax; ++j)
			{
				// verify if dF for t,t+1,t+2,t+3 are equal on these bits
				// if so skip this path
				if ( dFti == (dFt[j]&maskt) 
					&& dFtp1i == (dFtp1[j]&masktp1)
					&& dFtp2i == (dFtp2[j]&masktp2)
					&& dFtp3i == (dFtp3[j]&masktp3)
					)
				{
					isgood[j] = b+1;
					++countbaborted[b];
				}
			}
		}
#endif // DONT_SKIP_LOWERPATHS
	}
}
