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
#include <unordered_map>
#include <thread>
#include <mutex>
#include <random>
#include <atomic>

#define MD5DETAIL_INLINE_IMPL
#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/progress_display.hpp>

#include <boost/lexical_cast.hpp>

#include "main.hpp"

//#define CPUPERFORMANCE
#ifdef CPUPERFORMANCE
#include <hashclash/cpuperformance.hpp>
//uint64 cpu_step_t[64];
#define UPDATE(s) update_performance_counter __update(cpu_step_t[s]);
#ifdef __GNUC__
#include <sched.h>
#endif
#else
#define UPDATE(s)
#endif

#ifdef DOTESTCOUNTS
#define TESTCOUNT(s) ++testcounts[s];
#else
#define TESTCOUNT(s)
#endif


void textcoll_solver_t::prepare_block1()
{
	uint64_t N = 1ULL<<21;
	
	vector< halfstate_t > in, out;

	// Q17-Q13,m1
	generateQ13Q17(out, N);

	// m15 Q12 ... m12 Q9
	for (int t = 15; t >= 12; --t)
	{
		std::swap(in, out);
		randomize_vector(in);
		extend_step_bw(t, out, in, N);
	}

	std::swap(in, out);
	randomize_vector(in);
	extend_step_fw(17, out, in, N); // m6

	std::swap(in, out);
	randomize_vector(in);
	extend_step_m11(out, in, N);    // m11

	std::swap(in, out);
	randomize_vector(in);
	extend_step_fw(19, out, in, N); // m0

	std::swap(in, out);
	randomize_vector(in);
	extend_step_fw(20, out, in, N); // m5
	for (size_t i = 0; i < out.size(); )
	{
		if (MA.checkword(5, out[i].m[5]) && MA.checkword(5, out[i].m[5] + m_diff[5]))
		{
			++i;
			continue;
		}
		std::swap( out[i], out.back() );
		out.pop_back();
	}
	std::cout << "out.size() after m5+m_diff filtering: " << out.size() << std::endl;

	std::swap(in, out);
	randomize_vector(in);
	extend_step_m10(out, in, N); // m10, m15

	std::swap(in, out);
	randomize_vector(in);
	extend_step_fw(23, out, in, N); // m4
		
	std::cout << "Saving 'Q7Q24.bin.gz'..." << std::flush;
	save_gz(out, "Q7Q24", binary_archive);
		
	std::map<char, size_t> charcount;
	for (int wt : { 0, 1, 4, 5, 6, 10, 11, 12, 13, 14, 15 })
	{
		for (unsigned b = 0; b < 4; ++b)
		{
			charcount.clear();
			for (auto& S: out)
				++charcount[ char((S.m[wt]>>(8*b))&0xFF) ];
			std::cout << "m" << wt << "[" << b << "]=msg[" << (wt*4+b) << "]: ";
			std::vector< std::pair<char,size_t> > charcountsorted(charcount.begin(), charcount.end());
			std::sort(charcountsorted.begin(), charcountsorted.end(), 
				[](const std::pair<char,size_t>& l, const std::pair<char,size_t>& r)
				{
					return l.second > r.second; 
				});
			for (auto& cc : charcountsorted)
				std::cout << cc.first;
			std::cout << std::endl;
			std::cout << "m" << wt << "[" << b << "]=msg[" << (wt*4+b) << "]: ";
			std::cout.precision(4);
			for (auto& cc : charcountsorted)
				std::cout << " " << cc.first << ":" << double(cc.second)/double(out.size());
			std::cout << std::endl << std::endl;
		}
	}
}

void textcoll_solver_t::generateQ13Q17(vec_halfstate_t& out, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nGenerate Q13-Q17,m1:";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m1range = MA.word_range(1);
			const size_t m1rangecount = m1range.count(), Q17count = masked_value_QtQtm1(17, halfstate_t()).count();
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S.Qt(13) = masked_value_Qt(13,S).sample(rng);
				S.Qt(14) = masked_value_QtQtm1(14,S).sample(rng);
				S.Qt(15) = masked_value_QtQtm1(15,S).sample(rng);
				S.Qt(16) = masked_value_QtQtm1(16,S).sample(rng);
				if (!checkrotationQtQtp1(13, S)) continue;
				if (!checkrotationQtQtp1(14, S)) continue;
				if (!checkrotationQtQtp1(15, S)) continue;

				auto Q17 = masked_value_QtQtm1(17,S);
				if (Q17count < m1rangecount)
				{
					S.Qt(17) = Q17.sample(rng);
					S.computeWt(16);
					if (!MA.checkword(1, S.m[1]))
						continue;
				} else {
					S.m[1] = MA.sampleword(1, rng);
					S.computeQtp1(16);
					if (!Q17.check(S.Qt(17)))
					continue;
				}
				if (!checkrotationQtQtp1(16, S)) continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});

	std::cout << "Generate Q13-Q17,m1: attempts: 2^" << log(double(attempts))/log(2.0) << " success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_fw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	unsigned wt = md5_wt[t];
	std::cout << "\n======\nStep " << t << ": (Q" << (t+1) << ",m" << wt << "):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd,t,wt](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto wtrange = MA.word_range(wt);
			const size_t wtrangecount = wtrange.count(), Qtp1count = masked_value_QtQtm1(t+1, halfstate_t()).count();
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				auto Qtp1mv = masked_value_QtQtm1(t+1,S);
				if (Qtp1count < wtrangecount)
				{
					S.Qt(t+1) = Qtp1mv.sample(rng);
					S.computeWt(t);
					if (!MA.checkword(wt, S.m[wt]))
						continue;
				} else {
					S.m[wt] = MA.sampleword(wt, rng);
					S.computeQtp1(t);
					if (!Qtp1mv.check(S.Qt(t+1)))
						continue;
				}
				if (!checkrotationQtQtp1(t, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_bw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	unsigned wt = md5_wt[t];
	std::cout << "\n======\nStep " << t << ": (Q" << (t-3) << ",m" << wt << "):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd,t,wt](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto wtrange = MA.word_range(wt);
			const size_t wtrangecount = wtrange.count(), Qtm3count = masked_value_QtQtp1(t-3, halfstate_t()).count();
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				auto Qtm3mv = masked_value_QtQtp1(t-3,S);
				if (Qtm3count < wtrangecount)
				{
					S.Qt(t-3) = Qtm3mv.sample(rng);
					S.computeWt(t);
					if (!MA.checkword(wt, S.m[wt]))
						continue;
				} else {
					S.m[wt] = MA.sampleword(wt, rng);
					S.computeQtm3(t);
					if (!Qtm3mv.check(S.Qt(t-3)))
						continue;
				}
				if (!checkrotationQtQtp1(t-3, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_m11(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nStep: (m11,Q8,Q19):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m11range = MA.word_range(11);
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				S.m[11] = MA.sampleword(11, rng);
				S.computeQtp1(18);
				auto Q19mv = masked_value_QtQtm1(19, S);
				if (!Q19mv.check(S.Qt(19)))
					continue;
				S.computeQtm3(11);
				auto Q8mv = masked_value_QtQtp1(8, S);
				if (!Q8mv.check(S.Qt(8)))
					continue;
				if (!checkrotationQtQtp1(8, S))
					continue;
				if (!checkrotationQtQtp1(18, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step m11: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_m10(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nStep: (m10,Q7,Q22,Q23):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m11range = MA.word_range(11);
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				S.m[10] = MA.sampleword(10, rng);
				
				S.computeQtm3(10);
				auto Q7mv = masked_value_QtQtp1(7,S);
				if (!Q7mv.check(S.Qt(7)))
					continue;

				S.computeQtp1(21); // m10
				auto Q22mv = masked_value_QtQtm1(22,S);
				if (!Q22mv.check(S.Qt(22)))
					continue;

				S.computeQtp1(22); // m15
				auto Q23mv = masked_value_QtQtm1(23,S);
				if (!Q23mv.check(S.Qt(23)))
					continue;
				if (!checkrotationQtQtp1(7,S))
					continue;
				if (!checkrotationQtQtp1(21,S))
					continue;
				if (!checkrotationQtQtp1(22,S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step m10: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}


void textcoll_solver_t::filltables(const differentialpath& _diffpath)
{
	diffpath = _diffpath;
	
	for (int t = -3; t <= 64; ++t)
	{
		dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
		Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
		Qprev[offset+t] = Qprev2[offset+t] = 0;
	}

	// build tables
	for (int t = diffpath.tbegin(); t < diffpath.tend(); ++t)
	{
		dQ[offset+t] = diffpath[t].diff();
		Qprev[offset+t] = diffpath[t].prev() | diffpath[t].prevn();
		Qprev2[offset+t] = diffpath[t].prev2() | diffpath[t].prev2n();
		Qvaluemask[offset+t] = (~diffpath[t].set0()) | diffpath[t].set1()
							| diffpath[t].prev() | diffpath[t].prevn()
							| diffpath[t].prev2() | diffpath[t].prev2n();
		Qvalue[offset+t] = diffpath[t].set1() | diffpath[t].prevn() | diffpath[t].prev2n();
	}
	for (int t = 0; t < 64; ++t)
	{
		booleanfunction* F = 0;
		if (t < 16) F = & MD5_F_data;
		else if (t < 32) F = & MD5_G_data;
		else if (t < 48) F = & MD5_H_data;
		else F = & MD5_I_data;
		uint32 dF = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			bf_outcome outcome = F->outcome(diffpath(t,b), diffpath(t-1,b), diffpath(t-2,b));
			if (outcome.size()) {
				if (outcome[0] == bc_plus) 			dF += 1<<b;
				else if (outcome[0] == bc_minus)	dF -= 1<<b;
			} else throw std::runtime_error("ambiguous path!!");
		}
		dT[offset+t] = dQ[offset+t-3] + dF + m_diff[md5_wt[t]];
		dR[offset+t] = dQ[offset+t+1] - dQ[offset+t];
	}

#if 0
	// experimentally measure how many Q24 solutions you need on average for a full solution
	counter_exponential_print attempt("attempt"), ok("ok"), partial("partial");
	fullstate_t S1, S2;
	while (true)
	{
		++attempt;
		S1.Qt(21) = masked_value_Qt(21,S1).sample();
		S1.Qt(22) = masked_value_QtQtm1(22,S1).sample();
		S1.Qt(23) = masked_value_QtQtm1(23,S1).sample();
		S1.Qt(24) = masked_value_QtQtm1(24,S1).sample();
		for (int k = 0; k < 16; ++k)
		{
			S1.m[k] = xrng64();
			S2.m[k] = S1.m[k] + m_diff[k];
		}
		for (int t = 21; t <= 24; ++t)
			S2.Qt(t) = S1.Qt(t) + dQt(t);
		for (int t = 24; t < 64; ++t)
		{
			S1.computeQtp1(t);
			S2.computeQtp1(t);
		}
		if (S2.Qt(61)-S1.Qt(61) != uint32_t(1)<<31)
			continue;
		++partial;
		if (S2.Qt(62)-S1.Qt(62) != uint32_t(1)<<31)
			continue;
		if (S2.Qt(63)-S1.Qt(63) != uint32_t(1)<<31)
			continue;
		if (S2.Qt(64)-S1.Qt(64) != uint32_t(1)<<31)
			continue;
		if ( (S1.Qt(64)>>31) != (S1.Qt(63)>>31) )
			continue;
		if ( (S1.Qt(64)>>31) != (S1.Qt(62)>>31) )
			continue;
		++ok;
		std::cout << "p=" << log(double(ok())/double(attempt()))/log(2.0) << " partial";
		std::cout << "p=" << log(double(ok())/double(partial()))/log(2.0) << std::endl;
	}
#endif
}


