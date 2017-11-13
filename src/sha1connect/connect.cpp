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
#include <hashclash/saveload_bz2.hpp>
#include <hashclash/sha1detail.hpp>
#include <hashclash/sha1differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>


using namespace hashclash;
using namespace std;

booleanfunction* sha1bf(unsigned t)
{
	if (t < 20) return &SHA1_F1_data;
	else if (t < 40) return &SHA1_F2_data;
	else if (t < 60) return &SHA1_F3_data;
	return &SHA1_F4_data;
}

#define CHECK_TUNNEL_BITCONDITION(tc,pc) \
	if (tc != bc_constant) { \
		if (tc == bc_one && (pc == bc_zero || pc == bc_plus)) continue; \
		if (tc == bc_zero && (pc == bc_one || pc == bc_minus)) continue; \
	}
	
inline void connect_helper(unsigned val, unsigned diff, bitcondition& cond, uint32& out, uint32 in, unsigned b)
{
	if (val == 0) {
		if (diff) {
			cond = bc_minus;
			out = in + (1<<b);
		} else {
			cond = bc_constant;
			out = in;
		}
	} else {
		cond = bc_plus;
		out = in - (1<<b);
	}
}
inline void connect_helper(unsigned val, unsigned diff, uint32& out, uint32 in, unsigned b)
{
	if (val == 0) {
		if (diff)
			out = in + (1<<b);
		else
			out = in;
	} else
		out = in - (1<<b);
}
inline bitcondition connect_helper(unsigned val, unsigned diff)
{
	if (val == 0) {
		if (diff)
			return bc_minus;
		else
			return bc_constant;
	} else
		return bc_plus;
}


template<int steps, bool storefdata, bool storenewconds, bool compmincond>
void sha1_connect_thread::connectbits_01234(const pair<connect_bitdata,unsigned>& in2, map<connect_bitdata,unsigned>& out, 
				 unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, 
				 vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& outdata, vector<unsigned>& outmincond)
{
	const connect_bitdata& in = in2.first;
	const unsigned b0 = b;
	const unsigned b1 = (b+30)&31;
	const unsigned b2 = (b+28)&31;
	const unsigned b3 = (b+26)&31;
	const unsigned b4 = (b+24)&31;
/*	 connect_bitdata result;
	 bf_outcome bfo0, bfo1, bfo2, bfo3, bfo4;
	 bf_conditions bfc0, bfc1, bfc2, bfc3, bfc4;
	 bitcondition Qtm3b, Qtm2b, Qtm1b, Qtb, Qtp1b, Qtp2b, Qtp3b;
	 uint32 dft, dftp1, dftp2, dftp3, dftp4;
	 pair<byteconditions,byteconditions> newcond;
	 bitcondition Ttm3b, Ttm2b, Ttm1b, Ttb, Ttp1b, Ttp2b, Ttp3b;
*/	Ttm3b = tbc(t-3,(b0+2)&31);
	Ttm2b = tbc(t-2,(b0+2)&31);
	Ttm1b = tbc(t-1,b0);
	Ttb = tbc(t, b1);
	Ttp1b = tbc(t+1, b2);
	Ttp2b = tbc(t+2, b3);
	Ttp3b = tbc(t+3, b4);
	
	//update the running bitconditions
	result = in;
	result.rqtm2[1] = result.rqtm2[0];
	result.rqtm1[1] = result.rqtm1[0];
	result.rqt[1] = result.rqt[0];
	result.rqtp1[1] = result.rqtp1[0];

	Qtm3b = lower(t-3,(b0+2)&31);
	Qtm2b = lower(t-2,(b0+2)&31);
	Qtp2b = upper(t+2,b3); if (Qtp2b == bc_prev) { Qtp2b = bc_constant; } //if (storenewconds) throw; }
	Qtp3b = upper(t+3,b4); if (Qtp3b == bc_prev) { Qtp3b = bc_constant; } //if (storenewconds) throw; }
	unsigned dQtm1bit = (in.dQtm1>>b0)&1;
	unsigned dQtbit = (in.dQt>>b1)&1;
	unsigned dQtp1bit = (in.dQtp1>>b2)&1;
	unsigned mmtbit = (mmaskt>>b0)&1;
	unsigned mmtp1bit = (mmasktp1>>b1)&1;
	unsigned mmtp2bit = (mmasktp2>>b2)&1;
	unsigned mmtp3bit = (mmasktp3>>b3)&1;
	unsigned mmtp4bit = (mmasktp4>>b4)&1;
	for (unsigned dqtm1 = 0; dqtm1 <= dQtm1bit; ++dqtm1) {
		connect_helper(dqtm1, dQtm1bit, Qtm1b, result.dQtm1, in.dQtm1, b0);
		if (b0 == 1 && result.dQtm1 != dQtm1b2b31) continue;
		if (b0 == 26 && result.dQtm1 != dQtm1b27b31) continue;
		if (b0 == 31 && (Qtm1b == Qtm1b31not || Qtm1b == Qtm1b31not2)) continue;
		if (b0 == 31)
			bfo0 = msb_bf_outcome(*Ft, Qtm1b, Qtm2b, Qtm3b);
		else
			bfo0 = Ft->outcome(Qtm1b, Qtm2b, Qtm3b);			
		for (unsigned bfo0index = 0; bfo0index < bfo0.size(); ++bfo0index) {
			dft = in.dFt - bfo0(bfo0index, b0);
			for (unsigned mmt = 0; mmt <= mmtbit; ++mmt) {
				if (mmt == 1 && b0 == 31) continue;
				connect_helper(mmt, mmtbit, result.dFt, dft, b0);
				if (result.dFt & (1<<b0)) continue;
				if (b0 == 31)
					bfc0 = msb_bf_forwardconditions(*Ft, Qtm1b, Qtm2b, Qtm3b, bfo0[bfo0index]);
				else
					bfc0 = Ft->forwardconditions(Qtm1b, Qtm2b, Qtm3b, bfo0[bfo0index]);
				if ((Qtm1freemask & (1<<b0)) && bfc0.first != bc_constant) continue;
				if ((Qtm2freemask & (1<<((b0+2)&31))) && bfc0.second != bc_constant) continue;
				if ((Qtm3freemask & (1<<((b0+2)&31))) && bfc0.third != bc_constant) continue;
				CHECK_TUNNEL_BITCONDITION(Ttm3b,bfc0.third);
				CHECK_TUNNEL_BITCONDITION(Ttm2b,bfc0.second);
				CHECK_TUNNEL_BITCONDITION(Ttm1b,bfc0.first);
				result.rqtm2[0] = bfc0.second;				
				if (steps == 1) {
					if (storefdata) result.fqtm1[b0] = bfc0.first;
					if (storenewconds) {
						newcond.first.set(bc_constant, bc_constant, bc_constant, bc_constant, bc_constant, bc_constant, bfc0.third);
						newcond.second.set(bc_constant, bc_constant, bc_constant, bc_constant, connect_helper(mmt, mmtbit) );
						newconds.push_back(newcond);
						outdata.push_back(result);
						unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1);
						outmincond.push_back(curcond);
					}
					else
					{
						if (compmincond)
						{
							unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1);
							auto itb = out.insert(make_pair(result,curcond));
							if (!itb.second && curcond < itb.first->second)
								itb.first->second = curcond;
						}
						else
							out[result];
					}
					continue;
				}

				usedFtp1 = true;
				for (unsigned dqt = 0; dqt <= dQtbit; ++dqt) {
					connect_helper(dqt, dQtbit, Qtb, result.dQt, in.dQt, b1);
					if (b1 == 1 && result.dQt != dQtb2b31) continue;
					if (b1 == 26 && result.dQt != dQtb27b31) continue;
					if (b1 == 31 && (Qtb == Qtb31not || Qtb == Qtb31not2)) continue;
					if (b1 == 31) 
						bfo1 = msb_bf_outcome(*Ftp1, Qtb, bfc0.first, in.rqtm2[1]);
					else
						bfo1 = Ftp1->outcome(Qtb, bfc0.first, in.rqtm2[1]);
					for (unsigned bfo1index = 0; bfo1index < bfo1.size(); ++bfo1index) {
						dftp1 = in.dFtp1 - bfo1(bfo1index, b1);
						for (unsigned mmtp1 = 0; mmtp1 <= mmtp1bit; ++mmtp1) {
							connect_helper(mmtp1, mmtp1bit, result.dFtp1, dftp1, b1);
							if (mmtp1 == 1 && b1 == 31) continue;
							if (result.dFtp1 & (1<<b1)) continue;
							if (b1 == 31)
								bfc1 = msb_bf_forwardconditions(*Ftp1, Qtb, bfc0.first, in.rqtm2[1], bfo1[bfo1index]);
							else
								bfc1 = Ftp1->forwardconditions(Qtb, bfc0.first, in.rqtm2[1], bfo1[bfo1index]);
							if ((Qtfreemask & (1<<b1)) && bfc1.first != bc_constant) continue;
							if ((Qtm1freemask & (1<<b0)) && bfc1.second != bc_constant) continue;
							if ((Qtm2freemask & (1<<b0)) && bfc1.third != bc_constant) continue;
							CHECK_TUNNEL_BITCONDITION(Ttm1b,bfc1.second);
							CHECK_TUNNEL_BITCONDITION(Ttb,bfc1.first);
							result.rqtm1[0] = bfc1.second;							
							if (steps == 2) {
								if (storefdata) result.fqt[b1] = bfc1.first;
							
								if (storenewconds) 
								{
									newcond.first.set(bc_constant, bc_constant, bc_constant, bc_constant, bc_constant, bfc1.third, bfc0.third);
									newcond.second.set(bc_constant, bc_constant, bc_constant, connect_helper(mmtp1, mmtp1bit), connect_helper(mmt, mmtbit) );
									newconds.push_back(newcond);
									outdata.push_back(result);
									unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1);
									outmincond.push_back(curcond);
								}
								else
								{
									if (compmincond)
									{
										unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1);
										auto itb = out.insert(make_pair(result,curcond));
										if (!itb.second && curcond < itb.first->second)
											itb.first->second = curcond;
									} 
									else 
										out[result];
								}
								continue;
							}							

							usedFtp2 = true;
							for (unsigned dqtp1 = 0; dqtp1 <= dQtp1bit; ++dqtp1) {
								connect_helper(dqtp1, dQtp1bit, Qtp1b, result.dQtp1, in.dQtp1, b2);
								if (b2 == 1 && result.dQtp1 != dQtp1b2b31) continue;
								if (b2 == 26 && result.dQtp1 != dQtp1b27b31) continue;
								if (b2 == 31 && (Qtp1b == Qtp1b31not || Qtp1b == Qtp1b31not2)) continue;
								if (b2 == 31)
									bfo2 = msb_bf_outcome(*Ftp2, Qtp1b, bfc1.first, in.rqtm1[1]);
								else
									bfo2 = Ftp2->outcome(Qtp1b, bfc1.first, in.rqtm1[1]);
								for (unsigned bfo2index = 0; bfo2index < bfo2.size(); ++bfo2index) {
									dftp2 = in.dFtp2 - bfo2(bfo2index, b2);
									for (unsigned mmtp2 = 0; mmtp2 <= mmtp2bit; ++mmtp2) {
										if (mmtp2 == 1 && b2 == 31) continue;
										connect_helper(mmtp2, mmtp2bit, result.dFtp2, dftp2, b2);
										if (result.dFtp2 & (1<<b2)) continue;
										if (b2 == 31)
											bfc2 = msb_bf_forwardconditions(*Ftp2, Qtp1b, bfc1.first, in.rqtm1[1], bfo2[bfo2index]);
										else
											bfc2 = Ftp2->forwardconditions(Qtp1b, bfc1.first, in.rqtm1[1], bfo2[bfo2index]);
										if ((Qtp1freemask & (1<<b2)) && bfc2.first != bc_constant) continue;
										if ((Qtfreemask & (1<<b1)) && bfc2.second != bc_constant) continue;
										if ((Qtm1freemask & (1<<b1)) && bfc2.third != bc_constant) continue;
										CHECK_TUNNEL_BITCONDITION(Ttb,bfc2.second);
										CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc2.first);
										result.rqt[0] = bfc2.second;
										if (steps == 3) {
											if (storefdata) result.fqtp1[b2] = bfc2.first;
											if (storenewconds) 
											{
												newcond.first.set(bc_constant, bc_constant, bc_constant, bc_constant, bfc2.third, bfc1.third, bfc0.third);
												newcond.second.set(bc_constant, bc_constant, connect_helper(mmtp2, mmtp2bit), 
													connect_helper(mmtp1, mmtp1bit), connect_helper(mmt, mmtbit) );
												newconds.push_back(newcond);
												outdata.push_back(result);
												unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1);
												outmincond.push_back(curcond);
											}
											else
											{
												if (compmincond)
												{
													unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1);
													auto itb = out.insert(make_pair(result,curcond));
													if (!itb.second && curcond < itb.first->second)
														itb.first->second = curcond;
												} 
												else 
													out[result];
											}
											continue;
										}
													
										usedFtp3 = true;
										if (b3 == 31)
											bfo3 = msb_bf_outcome(*Ftp3, Qtp2b, bfc2.first, in.rqt[1]);
										else
											bfo3 = Ftp3->outcome(Qtp2b, bfc2.first, in.rqt[1]);
										for (unsigned bfo3index = 0; bfo3index < bfo3.size(); ++bfo3index) {
											dftp3 = in.dFtp3 - bfo3(bfo3index, b3);
											for (unsigned mmtp3 = 0; mmtp3 <= mmtp3bit; ++mmtp3) {
												if (mmtp3 == 1 && b3 == 31) continue;
												connect_helper(mmtp3, mmtp3bit, result.dFtp3, dftp3, b3);
												if (result.dFtp3 & (1<<b3)) continue;
												if (b3 == 31)
													bfc3 = msb_bf_forwardconditions(*Ftp3, Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
												else
													bfc3 = Ftp3->forwardconditions(Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
												if ((Qtp2freemask & (1<<b3)) && bfc3.first != bc_constant) continue;
												if ((Qtp1freemask & (1<<b2)) && bfc3.second != bc_constant) continue;
												if ((Qtfreemask & (1<<b2)) && bfc3.third != bc_constant) continue;
												CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc3.second);
												CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc3.first);
												result.rqtp1[0] = bfc3.second;
												if (steps == 4) {
													if (storefdata) result.fqtp2[b3] = bfc3.first;
													if (storenewconds) 
													{
														newcond.first.set(bc_constant, bc_constant, bc_constant, bfc3.third, bfc2.third, bfc1.third, bfc0.third);
														newcond.second.set(bc_constant, connect_helper(mmtp3, mmtp3bit),
															connect_helper(mmtp2, mmtp2bit), connect_helper(mmtp1, mmtp1bit), connect_helper(mmt, mmtbit) );
														newconds.push_back(newcond);
														outdata.push_back(result);
														unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1)  + (bfc3.third==bc_constant?0:1);
														outmincond.push_back(curcond);
													}
													else
													{
														if (compmincond)
														{
															unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1)  + (bfc3.third==bc_constant?0:1);
															auto itb = out.insert(make_pair(result,curcond));
															if (!itb.second && curcond < itb.first->second)
																itb.first->second = curcond;
														} 
														else 
															out[result];
													}
													continue;
												}												

												usedFtp4 = true;
												if (b4 == 31)
													bfo4 = msb_bf_outcome(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1]);
												else
													bfo4 = Ftp4->outcome(Qtp3b, bfc3.first, in.rqtp1[1]);
												for (unsigned bfo4index = 0; bfo4index < bfo4.size(); ++bfo4index) {
													dftp4 = in.dFtp4 - bfo4(bfo4index, b4);
													for (unsigned mmtp4 = 0; mmtp4 <= mmtp4bit; ++mmtp4) {
														if (mmtp4 == 1 && b4 == 31) continue;
														connect_helper(mmtp4, mmtp4bit, result.dFtp4, dftp4, b4);
														if (result.dFtp4 & (1<<b4)) continue;
															if (b4 == 31)
																bfc4 = msb_bf_forwardconditions(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
															else
																bfc4 = Ftp4->forwardconditions(Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
															if ((Qtp3freemask & (1<<b4)) && bfc4.first != bc_constant) continue;
															if ((Qtp2freemask & (1<<b3)) && bfc4.second != bc_constant) continue;
															if ((Qtp1freemask & (1<<b3)) && bfc4.third != bc_constant) continue;
															CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc4.second);
															CHECK_TUNNEL_BITCONDITION(Ttp3b,bfc4.first);
														if (storenewconds) 
														{
															newcond.first.set(bfc4.first, bfc4.second, bfc4.third, bfc3.third, bfc2.third, bfc1.third, bfc0.third);
															newcond.second.set( connect_helper(mmtp4, mmtp4bit), connect_helper(mmtp3, mmtp3bit),
																connect_helper(mmtp2, mmtp2bit), connect_helper(mmtp1, mmtp1bit), connect_helper(mmt, mmtbit) );
															newconds.push_back(newcond);
															outdata.push_back(result);
															unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1) + (bfc3.third==bc_constant?0:1)
																	 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
															outmincond.push_back(curcond);
														}
														else
														{
															if (compmincond)
															{
																unsigned curcond = in2.second + (bfc0.third==bc_constant?0:1) + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1) + (bfc3.third==bc_constant?0:1)
																		 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
																auto itb = out.insert(make_pair(result,curcond));
																if (!itb.second && curcond < itb.first->second)
																	itb.first->second = curcond;
															} 
															else 
																out[result];
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

template<bool storenewconds>
void sha1_connect_thread::connectbits_1234(const pair<connect_bitdata,unsigned>& in2, map<connect_bitdata,unsigned>& out, 
				 unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, 
				 vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& outdata, vector<unsigned>& outmincond)
{
	const connect_bitdata& in = in2.first;
	const unsigned b0 = b&31;
	const unsigned b1 = (b+30)&31;
	const unsigned b2 = (b+28)&31;
	const unsigned b3 = (b+26)&31;
	const unsigned b4 = (b+24)&31;
/*	 connect_bitdata result;
	 bf_outcome bfo1, bfo2, bfo3, bfo4;
	 bf_conditions bfc1, bfc2, bfc3, bfc4;
	 bitcondition Qtm2b, Qtb, Qtp1b, Qtp2b, Qtp3b;
	 uint32 dftp1, dftp2, dftp3, dftp4;
	 pair<byteconditions,byteconditions> newcond;
	 bitcondition Ttm1b, Ttb, Ttp1b, Ttp2b, Ttp3b;
*/
	Ttm1b = tbc(t-1,b0);
	Ttb = tbc(t, b1);
	Ttp1b = tbc(t+1, b2);
	Ttp2b = tbc(t+2, b3);
	Ttp3b = tbc(t+3, b4);

	//update the running bitconditions
	result = in;
	result.rqtm2[1] = result.rqtm2[0];
	result.rqtm1[1] = result.rqtm1[0];
	result.rqt[1] = result.rqt[0];
	result.rqtp1[1] = result.rqtp1[0];

	Qtm2b = lower(t-2,(b0+2)&31);
	Qtp2b = upper(t+2,b3); if (Qtp2b == bc_prev) { Qtp2b = bc_constant; } // if (storenewconds) throw; }
	Qtp3b = upper(t+3,b4); if (Qtp3b == bc_prev) { Qtp3b = bc_constant; } // if (storenewconds) throw; }
	unsigned dQtbit = (in.dQt>>b1)&1;
	unsigned dQtp1bit = (in.dQtp1>>b2)&1;
	unsigned mmtp1bit = (mmasktp1>>b1)&1;
	unsigned mmtp2bit = (mmasktp2>>b2)&1;
	unsigned mmtp3bit = (mmasktp3>>b3)&1;
	unsigned mmtp4bit = (mmasktp4>>b4)&1;

	result.rqtm2[0] = bc_constant;

	usedFtp1 = true;
	if (b0 != 0 && b0 != 1) throw;
	for (unsigned dqt = 0; dqt <= dQtbit; ++dqt) {
		connect_helper(dqt, dQtbit, Qtb, result.dQt, in.dQt, b1);
		if (b1 == 1 && result.dQt != dQtb2b31) continue;
		if (b1 == 26 && result.dQt != dQtb27b31) continue;
		if (b1 == 31 && (Qtb == Qtb31not || Qtb == Qtb31not2)) continue;
		if (b1 == 31)
			bfo1 = msb_bf_outcome(*Ftp1, Qtb, in.fqtm1[b0], in.rqtm2[1]);
		else
			bfo1 = Ftp1->outcome(Qtb, in.fqtm1[b0], in.rqtm2[1]);
		for (unsigned bfo1index = 0; bfo1index < bfo1.size(); ++bfo1index) {
			dftp1 = in.dFtp1 - bfo1(bfo1index, b1);
			for (unsigned mmtp1 = 0; mmtp1 <= mmtp1bit; ++mmtp1) {
				if (mmtp1 == 1 && b1 == 31) continue;
				connect_helper(mmtp1, mmtp1bit, result.dFtp1, dftp1, b1);
				if (result.dFtp1 & (1<<b1)) continue;
				if (b1 == 31)
					bfc1 = msb_bf_forwardconditions(*Ftp1, Qtb, in.fqtm1[b0], in.rqtm2[1], bfo1[bfo1index]);
				else
					bfc1 = Ftp1->forwardconditions(Qtb, in.fqtm1[b0], in.rqtm2[1], bfo1[bfo1index]);
				if ((Qtfreemask & (1<<b1)) && bfc1.first != bc_constant) continue;
				if ((Qtm1freemask & (1<<b0)) && bfc1.second != bc_constant) continue;
				if ((Qtm2freemask & (1<<b0)) && bfc1.third != bc_constant) continue;
				CHECK_TUNNEL_BITCONDITION(Ttm1b,bfc1.second);
				CHECK_TUNNEL_BITCONDITION(Ttb,bfc1.first);
				result.rqtm1[0] = bfc1.second;							

				usedFtp2 = true;
				for (unsigned dqtp1 = 0; dqtp1 <= dQtp1bit; ++dqtp1) {
					connect_helper(dqtp1, dQtp1bit, Qtp1b, result.dQtp1, in.dQtp1, b2);
					if (b2 == 1 && result.dQtp1 != dQtp1b2b31) continue;
					if (b2 == 26 && result.dQtp1 != dQtp1b27b31) continue;
					if (b2 == 31 && (Qtp1b == Qtp1b31not || Qtp1b == Qtp1b31not2)) continue;
					if (b2 == 31)
						bfo2 = msb_bf_outcome(*Ftp2, Qtp1b, bfc1.first, in.rqtm1[1]);
					else
						bfo2 = Ftp2->outcome(Qtp1b, bfc1.first, in.rqtm1[1]);
					for (unsigned bfo2index = 0; bfo2index < bfo2.size(); ++bfo2index) {
						dftp2 = in.dFtp2 - bfo2(bfo2index, b2);
						for (unsigned mmtp2 = 0; mmtp2 <= mmtp2bit; ++mmtp2) {
							if (mmtp2 == 1 && b2 == 31) continue;
							connect_helper(mmtp2, mmtp2bit, result.dFtp2, dftp2, b2);
							if (result.dFtp2 & (1<<b2)) continue;
							if (b2 == 31)
								bfc2 = msb_bf_forwardconditions(*Ftp2, Qtp1b, bfc1.first, in.rqtm1[1], bfo2[bfo2index]);
							else
								bfc2 = Ftp2->forwardconditions(Qtp1b, bfc1.first, in.rqtm1[1], bfo2[bfo2index]);
							if ((Qtp1freemask & (1<<b2)) && bfc2.first != bc_constant) continue;
							if ((Qtfreemask & (1<<b1)) && bfc2.second != bc_constant) continue;
							if ((Qtm1freemask & (1<<b1)) && bfc2.third != bc_constant) continue;
							CHECK_TUNNEL_BITCONDITION(Ttb,bfc2.second);
							CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc2.first);
							result.rqt[0] = bfc2.second;
																
							usedFtp3 = true;
							if (b3 == 31)
								bfo3 = msb_bf_outcome(*Ftp3, Qtp2b, bfc2.first, in.rqt[1]);
							else
								bfo3 = Ftp3->outcome(Qtp2b, bfc2.first, in.rqt[1]);
							for (unsigned bfo3index = 0; bfo3index < bfo3.size(); ++bfo3index) {
								dftp3 = in.dFtp3 - bfo3(bfo3index, b3);
								for (unsigned mmtp3 = 0; mmtp3 <= mmtp3bit; ++mmtp3) {
									if (mmtp3 == 1 && b3 == 31) continue;
									connect_helper(mmtp3, mmtp3bit, result.dFtp3, dftp3, b3);
									if (result.dFtp3 & (1<<b3)) continue;
									if (b3 == 31)
										bfc3 = msb_bf_forwardconditions(*Ftp3, Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
									else
										bfc3 = Ftp3->forwardconditions(Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
									if ((Qtp2freemask & (1<<b3)) && bfc3.first != bc_constant) continue;
									if ((Qtp1freemask & (1<<b2)) && bfc3.second != bc_constant) continue;
									if ((Qtfreemask & (1<<b2)) && bfc3.third != bc_constant) continue;
									CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc3.second);
									CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc3.first);
									result.rqtp1[0] = bfc3.second;

									usedFtp4 = true;
									if (b4 == 31)
										bfo4 = msb_bf_outcome(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1]);
									else
										bfo4 = Ftp4->outcome(Qtp3b, bfc3.first, in.rqtp1[1]);
									for (unsigned bfo4index = 0; bfo4index < bfo4.size(); ++bfo4index) {
										dftp4 = in.dFtp4 - bfo4(bfo4index, b4);
										for (unsigned mmtp4 = 0; mmtp4 <= mmtp4bit; ++mmtp4) {
											if (mmtp4 == 1 && b4 == 31) continue;
											connect_helper(mmtp4, mmtp4bit, result.dFtp4, dftp4, b4);
											if (result.dFtp4 & (1<<b4)) continue;
												if (b4 == 31)
													bfc4 = msb_bf_forwardconditions(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
												else
													bfc4 = Ftp4->forwardconditions(Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
												if ((Qtp3freemask & (1<<b4)) && bfc4.first != bc_constant) continue;
												if ((Qtp2freemask & (1<<b3)) && bfc4.second != bc_constant) continue;
												if ((Qtp1freemask & (1<<b3)) && bfc4.third != bc_constant) continue;
												CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc4.second);
												CHECK_TUNNEL_BITCONDITION(Ttp3b,bfc4.first);
											unsigned curcond = in2.second + (bfc1.third==bc_constant?0:1) + (bfc2.third==bc_constant?0:1) + (bfc3.third==bc_constant?0:1)
													 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
											if (storenewconds) 
											{
												newcond.first.set(bfc4.first, bfc4.second, bfc4.third, bfc3.third, bfc2.third, bfc1.third);
												newcond.second.set( connect_helper(mmtp4, mmtp4bit), connect_helper(mmtp3, mmtp3bit),
													connect_helper(mmtp2, mmtp2bit), connect_helper(mmtp1, mmtp1bit));
												newconds.push_back(newcond);
												outdata.push_back(result);
												outmincond.push_back(curcond);
											}
											else
											{
												auto itb = out.insert(make_pair(result,curcond));
												if (!itb.second && curcond < itb.first->second)
													itb.first->second = curcond;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

template<bool storenewconds>
void sha1_connect_thread::connectbits_234(const pair<connect_bitdata,unsigned>& in2, map<connect_bitdata,unsigned>& out, 
				 unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, 
				 vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& outdata, vector<unsigned>& outmincond)
{
	const connect_bitdata& in = in2.first;
	const unsigned b0 = b&31;
	const unsigned b1 = (b+30)&31;
	const unsigned b2 = (b+28)&31;
	const unsigned b3 = (b+26)&31;
	const unsigned b4 = (b+24)&31;
/*	 connect_bitdata result;
	 bf_outcome bfo2, bfo3, bfo4;
	 bf_conditions bfc2, bfc3, bfc4;
	 bitcondition Qtp1b, Qtp2b, Qtp3b;
	 uint32 dftp2, dftp3, dftp4;
	 pair<byteconditions,byteconditions> newcond;
	 bitcondition Ttb, Ttp1b, Ttp2b, Ttp3b;
*/
	Ttb = tbc(t, b1);
	Ttp1b = tbc(t+1, b2);
	Ttp2b = tbc(t+2, b3);
	Ttp3b = tbc(t+3, b4);
	//update the running bitconditions
	result = in;
	result.rqtm2[1] = result.rqtm2[0];
	result.rqtm1[1] = result.rqtm1[0];
	result.rqt[1] = result.rqt[0];
	result.rqtp1[1] = result.rqtp1[0];

	Qtp2b = upper(t+2,b3); if (Qtp2b == bc_prev) { Qtp2b = bc_constant; } // if (storenewconds) throw; }
	Qtp3b = upper(t+3,b4); if (Qtp3b == bc_prev) { Qtp3b = bc_constant; } // if (storenewconds) throw; }
	unsigned dQtp1bit = (in.dQtp1>>b2)&1;
	unsigned mmtp2bit = (mmasktp2>>b2)&1;
	unsigned mmtp3bit = (mmasktp3>>b3)&1;
	unsigned mmtp4bit = (mmasktp4>>b4)&1;

	result.rqtm2[0] = bc_constant;
	result.rqtm1[0] = bc_constant;

	usedFtp2 = true;
	if (b1 != 0 && b1 != 1) throw;
	for (unsigned dqtp1 = 0; dqtp1 <= dQtp1bit; ++dqtp1) {
		connect_helper(dqtp1, dQtp1bit, Qtp1b, result.dQtp1, in.dQtp1, b2);
		if (b2 == 1 && result.dQtp1 != dQtp1b2b31) continue;
		if (b2 == 26 && result.dQtp1 != dQtp1b27b31) continue;
		if (b2 == 31 && (Qtp1b == Qtp1b31not || Qtp1b == Qtp1b31not2)) continue;
		if (b2 == 31)
			bfo2 = msb_bf_outcome(*Ftp2, Qtp1b, in.fqt[b1], in.rqtm1[1]);
		else
			bfo2 = Ftp2->outcome(Qtp1b, in.fqt[b1], in.rqtm1[1]);
		for (unsigned bfo2index = 0; bfo2index < bfo2.size(); ++bfo2index) {
			dftp2 = in.dFtp2 - bfo2(bfo2index, b2);
			for (unsigned mmtp2 = 0; mmtp2 <= mmtp2bit; ++mmtp2) {
				if (mmtp2 == 1 && b2 == 31) continue;
				connect_helper(mmtp2, mmtp2bit, result.dFtp2, dftp2, b2);
				if (result.dFtp2 & (1<<b2)) continue;
				if (b2 == 31)
					bfc2 = msb_bf_forwardconditions(*Ftp2, Qtp1b, in.fqt[b1], in.rqtm1[1], bfo2[bfo2index]);
				else
					bfc2 = Ftp2->forwardconditions(Qtp1b, in.fqt[b1], in.rqtm1[1], bfo2[bfo2index]);
				if ((Qtp1freemask & (1<<b2)) && bfc2.first != bc_constant) continue;
				if ((Qtfreemask & (1<<b1)) && bfc2.second != bc_constant) continue;
				if ((Qtm1freemask & (1<<b1)) && bfc2.third != bc_constant) continue;
				CHECK_TUNNEL_BITCONDITION(Ttb,bfc2.second);
				CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc2.first);
				result.rqt[0] = bfc2.second;
											
				usedFtp3 = true;
				if (b3 == 31)
					bfo3 = msb_bf_outcome(*Ftp3, Qtp2b, bfc2.first, in.rqt[1]);
				else
					bfo3 = Ftp3->outcome(Qtp2b, bfc2.first, in.rqt[1]);
				for (unsigned bfo3index = 0; bfo3index < bfo3.size(); ++bfo3index) {
					dftp3 = in.dFtp3 - bfo3(bfo3index, b3);
					for (unsigned mmtp3 = 0; mmtp3 <= mmtp3bit; ++mmtp3) {
						if (mmtp3 == 1 && b3 == 31) continue;
						connect_helper(mmtp3, mmtp3bit, result.dFtp3, dftp3, b3);
						if (result.dFtp3 & (1<<b3)) continue;
						if (b3 == 31)
							bfc3 = msb_bf_forwardconditions(*Ftp3, Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
						else
							bfc3 = Ftp3->forwardconditions(Qtp2b, bfc2.first, in.rqt[1], bfo3[bfo3index]);
						if ((Qtp2freemask & (1<<b3)) && bfc3.first != bc_constant) continue;
						if ((Qtp1freemask & (1<<b2)) && bfc3.second != bc_constant) continue;
						if ((Qtfreemask & (1<<b2)) && bfc3.third != bc_constant) continue;
						CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc3.second);
						CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc3.first);
						result.rqtp1[0] = bfc3.second;

						usedFtp4 = true;
						if (b4 == 31)
							bfo4 = msb_bf_outcome(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1]);
						else
							bfo4 = Ftp4->outcome(Qtp3b, bfc3.first, in.rqtp1[1]);
						for (unsigned bfo4index = 0; bfo4index < bfo4.size(); ++bfo4index) {
							dftp4 = in.dFtp4 - bfo4(bfo4index, b4);
							for (unsigned mmtp4 = 0; mmtp4 <= mmtp4bit; ++mmtp4) {
								if (mmtp4 == 1 && b4 == 31) continue;
								connect_helper(mmtp4, mmtp4bit, result.dFtp4, dftp4, b4);
								if (result.dFtp4 & (1<<b4)) continue;
									if (b4 == 31)
										bfc4 = msb_bf_forwardconditions(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
									else
										bfc4 = Ftp4->forwardconditions(Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
									if ((Qtp3freemask & (1<<b4)) && bfc4.first != bc_constant) continue;
									if ((Qtp2freemask & (1<<b3)) && bfc4.second != bc_constant) continue;
									if ((Qtp1freemask & (1<<b3)) && bfc4.third != bc_constant) continue;
									CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc4.second);
									CHECK_TUNNEL_BITCONDITION(Ttp3b,bfc4.first);
								unsigned curcond = in2.second + (bfc2.third==bc_constant?0:1) + (bfc3.third==bc_constant?0:1)
										 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
								if (storenewconds) 
								{
									newcond.first.set(bfc4.first, bfc4.second, bfc4.third, bfc3.third, bfc2.third);
									newcond.second.set( connect_helper(mmtp4, mmtp4bit), connect_helper(mmtp3, mmtp3bit),
										connect_helper(mmtp2, mmtp2bit));
									newconds.push_back(newcond);
									outdata.push_back(result);
									outmincond.push_back(curcond);
								}
								else
								{
									auto itb = out.insert(make_pair(result,curcond));
									if (!itb.second && curcond < itb.first->second)
										itb.first->second = curcond;
								}
							}
						}
					}
				}
			}
		}
	}
}

template<bool storenewconds>
void sha1_connect_thread::connectbits_34(const pair<connect_bitdata,unsigned>& in2, map<connect_bitdata,unsigned>& out, 
				 unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, 
				 vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& outdata, vector<unsigned>& outmincond)
{
	const connect_bitdata& in = in2.first;
	const unsigned b0 = b&31;
	const unsigned b1 = (b+30)&31;
	const unsigned b2 = (b+28)&31;
	const unsigned b3 = (b+26)&31;
	const unsigned b4 = (b+24)&31;
/*	 connect_bitdata result;
	 bf_outcome bfo3, bfo4;
	 bf_conditions bfc3, bfc4;
	 bitcondition Qtp2b, Qtp3b;
	 uint32 dftp3, dftp4;
	 pair<byteconditions,byteconditions> newcond;
	 bitcondition Ttp1b, Ttp2b, Ttp3b;
*/
	Ttp1b = tbc(t+1, b2);
	Ttp2b = tbc(t+2, b3);
	Ttp3b = tbc(t+3, b4);
	//update the running bitconditions
	result = in;
	result.rqtm2[1] = result.rqtm2[0];
	result.rqtm1[1] = result.rqtm1[0];
	result.rqt[1] = result.rqt[0];
	result.rqtp1[1] = result.rqtp1[0];

	Qtp2b = upper(t+2,b3); if (Qtp2b == bc_prev) { Qtp2b = bc_constant; } // if (storenewconds) throw; }
	Qtp3b = upper(t+3,b4); if (Qtp3b == bc_prev) { Qtp3b = bc_constant; } // if (storenewconds) throw; }
	unsigned mmtp3bit = (mmasktp3>>b3)&1;
	unsigned mmtp4bit = (mmasktp4>>b4)&1;

	result.rqtm2[0] = bc_constant;
	result.rqtm1[0] = bc_constant;
	result.rqt[0] = bc_constant;

	usedFtp3 = true;
	if (b2 != 0 && b2 != 1) throw;
	if (b3 == 31)
		bfo3 = msb_bf_outcome(*Ftp3, Qtp2b, in.fqtp1[b2], in.rqt[1]);
	else
		bfo3 = Ftp3->outcome(Qtp2b, in.fqtp1[b2], in.rqt[1]);
	for (unsigned bfo3index = 0; bfo3index < bfo3.size(); ++bfo3index) {
		dftp3 = in.dFtp3 - bfo3(bfo3index, b3);
		for (unsigned mmtp3 = 0; mmtp3 <= mmtp3bit; ++mmtp3) {
			if (mmtp3 == 1 && b3 == 31) continue;
			connect_helper(mmtp3, mmtp3bit, result.dFtp3, dftp3, b3);
			if (result.dFtp3 & (1<<b3)) continue;
			if (b3 == 31)
				bfc3 = msb_bf_forwardconditions(*Ftp3, Qtp2b, in.fqtp1[b2], in.rqt[1], bfo3[bfo3index]);
			else
				bfc3 = Ftp3->forwardconditions(Qtp2b, in.fqtp1[b2], in.rqt[1], bfo3[bfo3index]);
			if ((Qtp2freemask & (1<<b3)) && bfc3.first != bc_constant) continue;
			if ((Qtp1freemask & (1<<b2)) && bfc3.second != bc_constant) continue;
			if ((Qtfreemask & (1<<b2)) && bfc3.third != bc_constant) continue;
			CHECK_TUNNEL_BITCONDITION(Ttp1b,bfc3.second);
			CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc3.first);
			result.rqtp1[0] = bfc3.second;

			usedFtp4 = true;
			if (b4 == 31)
				bfo4 = msb_bf_outcome(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1]);
			else
				bfo4 = Ftp4->outcome(Qtp3b, bfc3.first, in.rqtp1[1]);
			for (unsigned bfo4index = 0; bfo4index < bfo4.size(); ++bfo4index) {
				dftp4 = in.dFtp4 - bfo4(bfo4index, b4);
				for (unsigned mmtp4 = 0; mmtp4 <= mmtp4bit; ++mmtp4) {
					if (mmtp4 == 1 && b4 == 31) continue;
					connect_helper(mmtp4, mmtp4bit, result.dFtp4, dftp4, b4);
					if (result.dFtp4 & (1<<b4)) continue;
						if (b4 == 31)
							bfc4 = msb_bf_forwardconditions(*Ftp4, Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
						else
							bfc4 = Ftp4->forwardconditions(Qtp3b, bfc3.first, in.rqtp1[1], bfo4[bfo4index]);
						if ((Qtp3freemask & (1<<b4)) && bfc4.first != bc_constant) continue;
						if ((Qtp2freemask & (1<<b3)) && bfc4.second != bc_constant) continue;
						if ((Qtp1freemask & (1<<b3)) && bfc4.third != bc_constant) continue;
						CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc4.second);
						CHECK_TUNNEL_BITCONDITION(Ttp3b,bfc4.first);
					unsigned curcond = in2.second + (bfc3.third==bc_constant?0:1)
							 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
					if (storenewconds) 
					{
						newcond.first.set(bfc4.first, bfc4.second, bfc4.third, bfc3.third);
						newcond.second.set( connect_helper(mmtp4, mmtp4bit), connect_helper(mmtp3, mmtp3bit));
						newconds.push_back(newcond);
						outdata.push_back(result);
						outmincond.push_back(curcond);
					}
					else
					{
						auto itb = out.insert(make_pair(result,curcond));
						if (!itb.second && curcond < itb.first->second)
							itb.first->second = curcond;
					}
				}
			}
		}
	}
}

template<bool storenewconds>
void sha1_connect_thread::connectbits_4(const pair<connect_bitdata,unsigned>& in2, map<connect_bitdata,unsigned>& out, 
				 unsigned b, const sha1differentialpath& lower, const sha1differentialpath& upper, 
				 vector< pair<byteconditions,byteconditions> >& newconds, vector<connect_bitdata>& outdata, vector<unsigned>& outmincond)
{
	const connect_bitdata& in = in2.first;
	const unsigned b0 = b&31;
	const unsigned b1 = (b+30)&31;
	const unsigned b2 = (b+28)&31;
	const unsigned b3 = (b+26)&31;
	const unsigned b4 = (b+24)&31;
/*	connect_bitdata result;
	bf_outcome bfo4;
	bf_conditions bfc4;
	bitcondition Qtp2b, Qtp3b;
	uint32 dftp4;
	pair<byteconditions,byteconditions> newcond;
	bitcondition Ttp2b, Ttp3b;
*/
	Ttp2b = tbc(t+2, b3);
	Ttp3b = tbc(t+3, b4);
	//update the running bitconditions
	result = in;
	result.rqtm2[1] = result.rqtm2[0];
	result.rqtm1[1] = result.rqtm1[0];
	result.rqt[1] = result.rqt[0];
	result.rqtp1[1] = result.rqtp1[0];

	Qtp2b = upper(t+2,b3); if (Qtp2b == bc_prev) { Qtp2b = bc_constant; } // if (storenewconds) throw; }
	Qtp3b = upper(t+3,b4); if (Qtp3b == bc_prev) { Qtp3b = bc_constant; } // if (storenewconds) throw; }
	unsigned mmtp4bit = (mmasktp4>>b4)&1;

	result.rqtm2[0] = bc_constant;
	result.rqtm1[0] = bc_constant;
	result.rqt[0] = bc_constant;
	result.rqtp1[0] = bc_constant;

	usedFtp4 = true;
	if (b3 != 0 && b3 != 1) throw;
	if (b4 == 31)
		bfo4 = msb_bf_outcome(*Ftp4, Qtp3b, in.fqtp2[b3], in.rqtp1[1]);
	else
		bfo4 = Ftp4->outcome(Qtp3b, in.fqtp2[b3], in.rqtp1[1]);
	for (unsigned bfo4index = 0; bfo4index < bfo4.size(); ++bfo4index) {
		dftp4 = in.dFtp4 - bfo4(bfo4index, b4);
		for (unsigned mmtp4 = 0; mmtp4 <= mmtp4bit; ++mmtp4) {
			if (mmtp4 == 1 && b4 == 31) continue;
			connect_helper(mmtp4, mmtp4bit, result.dFtp4, dftp4, b4);
			if (result.dFtp4 & (1<<b4)) continue;
				if (b4 == 31)
					bfc4 = msb_bf_forwardconditions(*Ftp4, Qtp3b, in.fqtp2[b3], in.rqtp1[1], bfo4[bfo4index]);
				else
					bfc4 = Ftp4->forwardconditions(Qtp3b, in.fqtp2[b3], in.rqtp1[1], bfo4[bfo4index]);
				if ((Qtp3freemask & (1<<b4)) && bfc4.first != bc_constant) continue;
				if ((Qtp2freemask & (1<<b3)) && bfc4.second != bc_constant) continue;
				if ((Qtp1freemask & (1<<b3)) && bfc4.third != bc_constant) continue;
				CHECK_TUNNEL_BITCONDITION(Ttp2b,bfc4.second);
				CHECK_TUNNEL_BITCONDITION(Ttp3b,bfc4.first);
			unsigned curcond = in2.second
					 + (bfc4.first==bc_constant?0:1) + (bfc4.second==bc_constant?0:1) + (bfc4.third==bc_constant?0:1);
			if (storenewconds) 
			{
				newcond.first.set(bfc4.first, bfc4.second, bfc4.third);
				newcond.second.set( connect_helper(mmtp4, mmtp4bit) );
				newconds.push_back(newcond);
				outdata.push_back(result);
				outmincond.push_back(curcond);
			}
			else
			{
				auto itb = out.insert(make_pair(result,curcond));
				if (!itb.second && curcond < itb.first->second)
					itb.first->second = curcond;
			}
		}
	}
}

void set_sdr_bit(sdr& bsdr, unsigned b, bitcondition bc) {
	switch (bc) {
		case bc_constant:
			bsdr.mask &= ~uint32(1<<b);
			bsdr.sign &= bsdr.mask;
			break;
		case bc_plus:
			bsdr.mask |= 1<<b;
			bsdr.sign |= 1<<b;
			break;
		case bc_minus:
			bsdr.mask |= 1<<b;
			bsdr.sign &= ~uint32(1<<b);
			break;
		default:
			throw;
	}
}

unsigned sha1_connect_thread::connect_paths(const sha1differentialpath& lowerpath, const sha1differentialpath& upperpath, connect_bitdata& startdata, path_container& container, bool returnb40)
{
/*	static vector<connect_bitdata> bitdataresults[41];
	static vector<connect_bitdata> bitdatastart[40];
	static vector<connect_bitdata> bitdataend[40];
	static vector< unsigned > bitdatamincond[40];
	static vector< pair<byteconditions,byteconditions> >  bitdatanewcond[40];
	static vector< unsigned > bitdatanewcondhw[40];
	static vector< unsigned > mincond[41];	
*/

	bitdataresults[0].clear();
	bitdataresults[0][startdata];

	unsigned b = 0;
	while (b < 32) {
		bitdataresults[b+1].clear();
		usedFtp1 = usedFtp2 = usedFtp3 = usedFtp4 = false;
		if (b < 2)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<1,false,false,false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 4)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<2,false,false,false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 6)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<3,false,false,false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 8)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<4,false,false,false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 32)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<5,false,false,false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else throw;
		
		if (bitdataresults[b+1].size() == 0) 
			return b;
		// remove results that have more conditions than our current limit
//		for (auto it = bitdataresults[b+1].begin(); it != bitdataresults[b+1].end(); )
//			if (it->second > container.bestpathcond)
//				it = bitdataresults[b+1].erase(it);
//			else
//				++it;
			
		++b;
	}

	unsigned mincond0 = 0;
	for (int k = lowerpath.tbegin(); k < t-3 && k < lowerpath.tend(); ++k)
		mincond0 += lowerpath[k].hw();
	for (int k = t+4; k < upperpath.tend(); ++k)
		mincond0 += upperpath[k].hw();
	bitdataresults[0].clear();
	bitdataresults[0][startdata]=mincond0;

	b = 0;
	while (b < 32+8) {		
		bitdataresults[b+1].clear();
		usedFtp1 = usedFtp2 = usedFtp3 = usedFtp4 = false;
		if (b < 2)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<1,true,false,true>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 4)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<2,true,false,true>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 6)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<3,true,false,true>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 8)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<4,true,false,true>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 32)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_01234<5,true,false,true>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 34)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_1234<false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 36)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_234<false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 38)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_34<false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else if (b < 40)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
				connectbits_4<false>(*it, bitdataresults[b+1], b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
		else throw;
		if (bitdataresults[b+1].size() == 0)
			break;

		// remove results that have more conditions than our current limit
		for (auto it = bitdataresults[b+1].begin(); it != bitdataresults[b+1].end(); )
			if (it->second > container.bestpathcond)
				it = bitdataresults[b+1].erase(it);
			else
				++it;
			
		++b;
	}
	if (b < 32+8) return b;
	if (returnb40) return 40;
	
	unsigned overalmincond = bitdataresults[40].begin()->second;
	
//	mincond[0].clear();
//	mincond[0].push_back(mincond0);

	for (b = 0; b < 40; ++b)
	{
		bitdatastart[b].clear();
		bitdataend[b].clear();
		bitdatamincond[b].clear();
		bitdatanewcond[b].clear();
		bitdatanewcondhw[b].clear();
		map<connect_bitdata,unsigned> tmpmap;
//		mincond[b+1].clear();
//		mincond[b+1].resize(bitdataresults[b+1].size(),262144);
		if (b < 2)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
			{
				tmpmap.clear();
				connectbits_01234<1,true,true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 4)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
			{
				tmpmap.clear();
				connectbits_01234<2,true,true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 6)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it)
			{
				tmpmap.clear();
				connectbits_01234<3,true,true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 8)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_01234<4,true,true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 32)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_01234<5,true,true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 34)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_1234<true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 36)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_234<true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 38)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_34<true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		else if (b < 40)
			for (auto it = bitdataresults[b].begin(); it != bitdataresults[b].end(); ++it) {
				tmpmap.clear();
				connectbits_4<true>(*it, tmpmap, b, lowerpath, upperpath, bitdatanewcond[b], bitdataend[b], bitdatamincond[b]);
				if (bitdatanewcond[b].size() != bitdataend[b].size()) { cerr << "#bitdatanewcond!=#bitdataend @ b=" << b << endl; throw; }
				bitdatastart[b].resize(bitdataend[b].size(), it->first);
			}
		bitdatanewcondhw[b].resize(bitdatanewcond[b].size());
		for (unsigned k = 0; k < bitdatanewcond[b].size(); ++k)
			bitdatanewcondhw[b][k] = bitdatanewcond[b][k].first.hw();
		if (bitdatastart[b].size() != bitdataend[b].size()
			|| bitdataend[b].size() != bitdatanewcond[b].size())
			throw;
	}
//	unsigned overalmincond = mincond[40][0];
	if (container.determinelowestcond) {
		if (overalmincond < container.bestpathcond) {
			container.bestpathcond = overalmincond;
			cout << "Best path: totcond=" << container.bestpathcond << " count=?" << endl;
		}
		return 40;
	}
	if (overalmincond > container.bestpathcond) 
		return 40;
	unsigned bestcond2 = 262144;

	/*static*/ sha1differentialpath newpath2;
	newpath2 = lowerpath;
	for (int k = upperpath.tbegin(); k < upperpath.tend(); ++k) {
		newpath2[k] = upperpath[k];
		newpath2.getme(k) = upperpath.getme(k);
	}
	unsigned bindex[40];
	unsigned bcond[40];
	unsigned bit = 39;
	bindex[bit] = 0;
	bcond[bit] = 0;
	//static uint64 badcnt = 0, okcnt = 0;
	while (bit <= 39)
	{
		if (bindex[bit] < bitdataend[bit].size())
		{			
			// check and update # conds
//			if (bcond[bit] + bitdatamincond[bit][bindex[bit]] > container.bestpathcond) {
//				--bit;
//				bindex[bit] = bitdataend[bit].size();
//				continue;
//			}
			bcond[bit-1] = bcond[bit] + bitdatanewcondhw[bit][bindex[bit]];
			// set bitconditions for committed bit
			if (bit < 32) {
				newpath2.setbitcondition(t-3,(bit+2)&31,bitdatanewcond[bit][bindex[bit]].first[6]);
				set_sdr_bit(newpath2.getme(t),(bit+0)&31,bitdatanewcond[bit][bindex[bit]].second[4]);
			}
			if (bit >= 2 && bit < 32+2) {
				newpath2.setbitcondition(t-2,(bit+0)&31,bitdatanewcond[bit][bindex[bit]].first[5]);
				set_sdr_bit(newpath2.getme(t+1),(bit+30)&31,bitdatanewcond[bit][bindex[bit]].second[3]);
			}
			if (bit >= 4 && bit < 32+4) {
				newpath2.setbitcondition(t-1,(bit+30)&31,bitdatanewcond[bit][bindex[bit]].first[4]);
				set_sdr_bit(newpath2.getme(t+2),(bit+28)&31,bitdatanewcond[bit][bindex[bit]].second[2]);
			}
			if (bit >= 6 && bit < 32+6) {
				newpath2.setbitcondition(t+0,(bit+28)&31,bitdatanewcond[bit][bindex[bit]].first[3]);
				set_sdr_bit(newpath2.getme(t+3),(bit+26)&31,bitdatanewcond[bit][bindex[bit]].second[1]);
			}
			if (bit >= 8 && bit < 32+8) {
				newpath2.setbitcondition(t+3,(bit+24)&31,bitdatanewcond[bit][bindex[bit]].first[0]);
				newpath2.setbitcondition(t+2,(bit+26)&31,bitdatanewcond[bit][bindex[bit]].first[1]);
				newpath2.setbitcondition(t+1,(bit+26)&31,bitdatanewcond[bit][bindex[bit]].first[2]);
				set_sdr_bit(newpath2.getme(t+4),(bit+24)&31,bitdatanewcond[bit][bindex[bit]].second[0]);
			}

			--bit;
			// if bit==0 then create path
			if (bit != 0) {
				bindex[bit] = bitdataend[bit].size();
				for (unsigned k = 0; k < bitdataend[bit].size(); ++k)
					if (bitdataend[bit][k] == bitdatastart[bit+1][bindex[bit+1]])
					{
						bindex[bit] = k;
						break;
					}
			} else {
				for (unsigned k = 0; k < bitdataend[0].size(); ++k)
					if (bitdataend[0][k] == bitdatastart[1][bindex[1]])
					{ // we have a full path !!!
						if (bcond[0] + bitdatanewcondhw[0][k] + mincond0 > container.bestpathcond) continue;
						bindex[0] = k;
						for (unsigned b = 0; b < 1; ++b)
						{
							if (b < 32) {
								newpath2.setbitcondition(t-3,(b+2)&31,bitdatanewcond[b][bindex[b]].first[6]);
								set_sdr_bit(newpath2.getme(t),(b+0)&31,bitdatanewcond[b][bindex[b]].second[4]);
							}
							if (b >= 2 && b < 32+2) {
								newpath2.setbitcondition(t-2,(b+0)&31,bitdatanewcond[b][bindex[b]].first[5]);
								set_sdr_bit(newpath2.getme(t+1),(b+30)&31,bitdatanewcond[b][bindex[b]].second[3]);
							}
							if (b >= 4 && b < 32+4) {
								newpath2.setbitcondition(t-1,(b+30)&31,bitdatanewcond[b][bindex[b]].first[4]);
								set_sdr_bit(newpath2.getme(t+2),(b+28)&31,bitdatanewcond[b][bindex[b]].second[2]);
							}
							if (b >= 6 && b < 32+6) {
								newpath2.setbitcondition(t+0,(b+28)&31,bitdatanewcond[b][bindex[b]].first[3]);
								set_sdr_bit(newpath2.getme(t+3),(b+26)&31,bitdatanewcond[b][bindex[b]].second[1]);
							}
							if (b >= 8 && b < 32+8) {
								newpath2.setbitcondition(t+3,(b+24)&31,bitdatanewcond[b][bindex[b]].first[0]);
								newpath2.setbitcondition(t+2,(b+26)&31,bitdatanewcond[b][bindex[b]].first[1]);
								newpath2.setbitcondition(t+1,(b+26)&31,bitdatanewcond[b][bindex[b]].first[2]);
								set_sdr_bit(newpath2.getme(t+4),(b+24)&31,bitdatanewcond[b][bindex[b]].second[0]);
							}
						}
						container.push_back(newpath2);						
					}
				bindex[bit] = bitdataend[bit].size();
			}
		} else
		{
			++bit;
			if (bit <= 39)
			{
				++bindex[bit];
				if (bit < 39)
				{
					while (bindex[bit] < bitdataend[bit].size()
						&& ((bitdataend[bit][bindex[bit]]) != (bitdatastart[bit+1][bindex[bit+1]])))
						++bindex[bit];
					if (bindex[bit] < bitdataend[bit].size())
						bcond[bit] = bcond[bit+1] + bitdatanewcond[bit+1][bindex[bit+1]].first.hw();
				}
			}
		}
	}
	return 40;
}

void sha1_connect_thread::sha1_connect(const sha1differentialpath& lowerpath
				 , const vector<sha1differentialpath>& upperpaths
				 , path_container& container, unsigned prev_eq_b)
{
	if (lowerpath.path.size() == 0) {
		cout << "lowerpath empty!" << endl;
		return;
	}
	if (dpFt.size() == 0) {
		t = container.t;
		mmaskt = container.m_mask[t];
		mmasktp1 = container.m_mask[t+1];
		mmasktp2 = container.m_mask[t+2];
		mmasktp3 = container.m_mask[t+3];
		mmasktp4 = container.m_mask[t+4];

		Ft = sha1bf(t);
		Ftp1 = sha1bf(t+1);
		Ftp2 = sha1bf(t+2);
		Ftp3 = sha1bf(t+3);
		Ftp4 = sha1bf(t+4);

		dpFt.resize(upperpaths.size());
		dpFtp1.resize(upperpaths.size());
		dpFtp2.resize(upperpaths.size());
		dpFtp3.resize(upperpaths.size());
		dpFtp4.resize(upperpaths.size());
		dQtp1.resize(upperpaths.size());
		for (unsigned i = 0; i < upperpaths.size(); ++i) {
			dpFt[i] = upperpaths[i][t+1].diff();
			dpFtp1[i] = upperpaths[i][t+2].diff() - upperpaths[i][t+1].getsdr().rotate_left(5).adddiff();
			dpFtp2[i] = upperpaths[i][t+3].diff() - upperpaths[i][t+2].getsdr().rotate_left(5).adddiff();
			dpFtp3[i] = upperpaths[i][t+4].diff() - upperpaths[i][t+3].getsdr().rotate_left(5).adddiff();
			dpFtp4[i] = upperpaths[i][t+5].diff() - upperpaths[i][t+4].getsdr().rotate_left(5).adddiff();
			dQtp1[i] = upperpaths[i][t+1].diff();
		}
		prevb.resize(upperpaths.size(),40);

		tbc = container.tunnelconditions;
		tbc.get(t-3); tbc.get(t+3);
		Qtm3freemask = Qtm2freemask = Qtm1freemask = Qtfreemask = Qtp1freemask = Qtp2freemask = Qtp3freemask = 0;
		for (unsigned b = 0; b < 32; ++b) {
			if (container.tunnelconditions(t-3,b) == bc_plus || container.tunnelconditions(t-3,b) == bc_minus) Qtm3freemask |= 1<<b;
			if (container.tunnelconditions(t-2,b) == bc_plus || container.tunnelconditions(t-2,b) == bc_minus) Qtm2freemask |= 1<<b;
			if (container.tunnelconditions(t-1,b) == bc_plus || container.tunnelconditions(t-1,b) == bc_minus) Qtm1freemask |= 1<<b;
			if (container.tunnelconditions(t,b) == bc_plus || container.tunnelconditions(t,b) == bc_minus) Qtfreemask |= 1<<b;
			if (container.tunnelconditions(t+1,b) == bc_plus || container.tunnelconditions(t+1,b) == bc_minus) Qtp1freemask |= 1<<b;
			if (container.tunnelconditions(t+2,b) == bc_plus || container.tunnelconditions(t+2,b) == bc_minus) Qtp2freemask |= 1<<b;
			if (container.tunnelconditions(t+3,b) == bc_plus || container.tunnelconditions(t+3,b) == bc_minus) Qtp3freemask |= 1<<b;
		}
	}
	uint32 dlFt = 0 - lowerpath[t].getsdr().rotate_left(5).adddiff() - lowerpath[t-4].getsdr().rotate_left(30).adddiff();
	uint32 dlFtp1 = 0 - lowerpath[t-3].getsdr().rotate_left(30).adddiff();
	uint32 dlFtp2 = 0 - lowerpath[t-2].getsdr().rotate_left(30).adddiff();
	uint32 dlFtp3 = 0 - lowerpath[t-1].getsdr().rotate_left(30).adddiff();
	uint32 dlFtp4 = 0 - lowerpath[t].getsdr().rotate_left(30).adddiff();

	uint32 dQtm1 = lowerpath[t-1].diff();
	uint32 dQt = lowerpath[t].diff();

	vector<connect_bitdata> bitdataresults[41];
	connect_bitdata startdata;
	startdata.dQtm1 = dQtm1;
	startdata.dQt = dQt;
	for (unsigned j = 0; j < 2; ++j)
	{
		startdata.rqtm2[j] = startdata.rqtm1[j] = startdata.rqt[j] = startdata.rqtp1[j] = bc_constant;
		startdata.fqtm1[j] = startdata.fqt[j] = startdata.fqtp1[j] = startdata.fqtp2[j] = bc_constant;
	}

	Qtm1b31not = Qtb31not = bc_plus;
	Qtm1b31not2 = Qtb31not2 = bc_constant;
	for (unsigned b = 27; b <= 31; ++b) {
		if (lowerpath(t-1,b) == bc_minus) Qtm1b31not = bc_plus;
		if (lowerpath(t-1,b) == bc_plus) Qtm1b31not = bc_minus;
		if (lowerpath(t,b) == bc_minus) Qtb31not = bc_plus;
		if (lowerpath(t,b) == bc_plus) Qtb31not = bc_minus;
	}
	if (lowerpath(t-1,31) == bc_constant) Qtm1b31not2 = Qtm1b31not;
	if (lowerpath(t,31) == bc_constant) Qtb31not2 = Qtb31not;
	sdr sdrQtm1 = lowerpath[t-1].getsdr();
	sdrQtm1.mask &= uint32(0)-uint32(1<<2); sdrQtm1.sign &= sdrQtm1.mask;	
	dQtm1b2b31 = sdrQtm1.adddiff();
	sdrQtm1.mask &= uint32(0)-uint32(1<<27); sdrQtm1.sign &= sdrQtm1.mask;
	dQtm1b27b31 = sdrQtm1.adddiff();
	sdr sdrQt = lowerpath[t].getsdr();
	sdrQt.mask &= uint32(0)-uint32(1<<2); sdrQt.sign &= sdrQt.mask;	
	dQtb2b31 = sdrQt.adddiff();
	sdrQt.mask &= uint32(0)-uint32(1<<27); sdrQt.sign &= sdrQt.mask;
	dQtb27b31 = sdrQt.adddiff();

	bool connected = false;
	unsigned i = 0;
	unsigned testcount = 0;
	unsigned highestb = 0;
	/*static vector<uint64> bcnt(41,0);*/
	while (i < upperpaths.size()) {
		Qtp1b31not = bc_plus;
		Qtp1b31not2 = bc_constant;
		for (unsigned b = 27; b <= 31; ++b) {
			if (upperpaths[i](t+1,b) == bc_minus) Qtp1b31not = bc_plus;
			if (upperpaths[i](t+1,b) == bc_plus) Qtp1b31not = bc_minus;
		}
		if (upperpaths[i](t+1,31) == bc_constant) Qtp1b31not2 = Qtp1b31not;
		sdr sdrQtp1 = upperpaths[i][t+1].getsdr();
		sdrQtp1.mask &= uint32(0)-uint32(1<<2); sdrQtp1.sign &= sdrQtp1.mask;	
		dQtp1b2b31 = sdrQtp1.adddiff();
		sdrQtp1.mask &= uint32(0)-uint32(1<<27); sdrQtp1.sign &= sdrQtp1.mask;
		dQtp1b27b31 = sdrQtp1.adddiff();

		bitdataresults[0].clear();
		startdata.dQtp1 = dQtp1[i];
		startdata.dFt = dpFt[i] + dlFt;
		startdata.dFtp1 = dpFtp1[i] + dlFtp1;
		startdata.dFtp2 = dpFtp2[i] + dlFtp2;
		startdata.dFtp3 = dpFtp3[i] + dlFtp3;
		startdata.dFtp4 = dpFtp4[i] + dlFtp4;
		unsigned b = prevb[i];
		if (b >= prev_eq_b)
			b = connect_paths(lowerpath, upperpaths[i], startdata, container, false);
		prevb[i] = b;
		if (0) { //b == 40) {
			cout << "?" << flush; 
			mmaskt = container.m_mask[t];
			mmasktp1 = container.m_mask[t+1];
			mmasktp2 = container.m_mask[t+2];
			mmasktp3 = container.m_mask[t+3];
			mmasktp4 = container.m_mask[t+4];
			b = connect_paths(lowerpath, upperpaths[i], startdata, container, false);
			mmaskt = mmasktp1 = mmasktp2 = mmasktp3 = mmasktp4 = 0;
		}
		if (b < 32+8) {			
			unsigned j = i+1;
			uint32 ftmask = ~uint32(0); if (b < 32) ftmask >>= (31-b);
			uint32 ftp1mask = ~uint32(0); 
			if (b < 2) ftp1mask = 0; else if (b-2 < 32) { ftp1mask >>=(31-(b-2)); if (!usedFtp1) ftp1mask >>=1; }
			uint32 qtp1mask = ~uint32(0); 
			if (b < 4) qtp1mask = 0; else if (b-4 < 32) { qtp1mask >>=(31-(b-4)); if (!usedFtp2) qtp1mask >>=1; }
			uint32 ftp2mask = ~uint32(0); 
			if (b < 4) ftp2mask = 0; else if (b-4 < 32) { ftp2mask >>=(31-(b-4)); if (!usedFtp2) ftp2mask >>=1; }
			unsigned qtp2bits = 32; 
			if (b < 6) qtp2bits = 0; else if (b-6 < 32) { qtp2bits = b-5; if (!usedFtp3) qtp2bits -= 1; }
			uint32 ftp3mask = ~uint32(0); 
			if (b < 6) ftp3mask = 0; else if (b-6 < 32) { ftp3mask >>=(31-(b-6)); if (!usedFtp3) ftp3mask >>=1; }
			unsigned qtp3bits = 32; 
			if (b < 8) qtp3bits = 0; else if (b-8 < 32) { qtp3bits = b-7; if (!usedFtp4) qtp3bits -= 1; }
			uint32 ftp4mask = ~uint32(0); 
			if (b < 8) ftp4mask = 0; else if (b-8 < 32) { ftp4mask >>=(31-(b-8)); if (!usedFtp4) ftp4mask >>=1; }
			
			while (j < upperpaths.size()) {
				if ((dpFt[i]^dpFt[j])&ftmask) break;
				if ((dpFtp1[i]^dpFtp1[j])&ftp1mask) break;
				if ((dQtp1[i]^dQtp1[j])&qtp1mask) break;
				if ((dpFtp2[i]^dpFtp2[j])&ftp2mask) break;
				if ((dpFtp3[i]^dpFtp3[j])&ftp3mask) break;
				if ((dpFtp4[i]^dpFtp4[j])&ftp4mask) break;
				bool bcok = true;
				for (unsigned k = 0; bcok && k < qtp2bits; ++k)
					if (upperpaths[i](t+2,k) != upperpaths[j](t+2,k))
						bcok = false;
				for (unsigned k = 0; bcok && k < qtp3bits; ++k)
					if (upperpaths[i](t+3,k) != upperpaths[j](t+3,k))
						bcok = false;
				if (!bcok) break;
				prevb[j]=b;
				++j;
			}
			i = j;
		} else {
			++i;
			connected = true;
		}
		if (b > highestb) highestb = b;
		++bcnt[b];
	}
	if (connected)
		cout << "+" << flush;
	/*static timer sw(true);*/
	if (container.showstats && sw.time() > 3600) {
		cout << endl;
		for (unsigned i = 0; i < bcnt.size(); ++i)
			cout << i << ": " << bcnt[i] << endl;
		sw.start();
	}
	/*static timer sw2(true);*/
	if (sw2.time() > 600) {
		container.save_bestpaths();
		sw2.start();
	}
}
