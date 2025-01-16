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


void textcoll_solver_t::start_block1()
{
	vector< halfstate_t > Q7Q24states;
	
	std::cout << "Trying to load 'Q7Q24.bin.gz'..." << std::flush;
	try {
		load_gz(Q7Q24states, "Q7Q24", binary_archive);
		std::cout << "done: " << Q7Q24states.size() << "." << std::endl;
	} catch (std::exception& e)
	{
		Q7Q24states.clear();
		std::cout << "failed." << std::endl;
		return;
	}

	randomize_vector(Q7Q24states);
	
	// sort on maximum Q9m9 tunnel	
	std::sort(Q7Q24states.begin(), Q7Q24states.end(), 
		[this](const halfstate_t& l, const halfstate_t& r)
		{
			return Q9m9tunnel(l) > Q9m9tunnel(r);
		});
	std::cout << "Q9m9 tunnel best state: " << Q9m9tunnel(Q7Q24states.front()) << std::endl;

	size_t skipstates = 0;
	for (auto& Q7Q24state : Q7Q24states)
	{
		if (skipstates > 0)
		{
			--skipstates;
			continue;
		}
		completeQ7Q24(Q7Q24state);
	}
}

void textcoll_solver_t::completeQ7Q24(const halfstate_t& Q7Q24state)
{
	fullstate_t S = Q7Q24state;

	// sanity check S
	for (int t = 10; t <= 23; ++t)
	{
		unsigned wt = md5_wt[t];
		uint32 Wtorg = S.m[wt];
		S.computeWt(t);
		if (S.m[wt] != Wtorg)
			throw std::runtime_error("Q7Q24state: incorrect Wt");
		if (!MA.checkword(wt, S.m[wt]))
			throw std::runtime_error("Q7Q24state: Wt alphabet violation");
	}

	std::cout << "MSG: ";
	for (unsigned t = 0; t < 16; ++t)
		for (unsigned b = 0; b < 4; ++b)
			std::cout << char( (S.m[t]>>(b*8))&0xFF);
	std::cout << std::endl;


	// iterate over all m4 values in random order
	if (m4rndrange.empty())
	{
		auto m4range = MA.word_range(4);
		m4rndrange.reserve(m4range.count());
		for (uint32_t m4 : m4range)
			m4rndrange.emplace_back(m4);
	}
	randomize_vector(m4rndrange);

	size_t m4attempts = 0, m4ok = 0;
	for (auto& m4 :	m4rndrange)
	{
		if (++m4attempts > 32 && m4ok == 0)
			return;
		S.m[4] = m4;
		compute_good_m10(S);
		std::cout << "m4=" << word2str(m4) << ": good_m10 size: " << good_m10.size() << std::endl;
		if (good_m10.size() < 6000) //(1ULL<<14))
			continue;
		++m4ok;
		
		std::cout << "MSG: ";
		for (unsigned t = 0; t < 16; ++t)
			for (unsigned b = 0; b < 4; ++b)
				std::cout << char( (S.m[t]>>(b*8))&0xFF);
		std::cout << std::endl;

		uint32_t m13org = S.m[13], Q10org = S.Qt(10);

		std::vector< std::pair<uint32_t,uint32_t> > m13Q10good;
		for (uint32_t m13 : MA.word_range(13))
		{
			S.m[13] = m13;
			S.computeQtm3(13);
			if (!masked_value_QtQtp1(10,S).check(S.Qt(10)))
				continue;
			if (!checkrotationQtQtp1(10,S))
				continue;
			m13Q10good.emplace_back(S.m[13], S.Qt(10));
		}
		std::cout << "m13Q10good size: " << m13Q10good.size() << std::endl;
		uint32_t Q9m9 = ~Qvaluemask[offset+9] & ~Qprev[offset+10] & S.Qt(11);
		std::sort(m13Q10good.begin(), m13Q10good.end(), 
			[Q9m9](const std::pair<uint32_t,uint32_t>& l, const std::pair<uint32_t,uint32_t>& r)
			{
				return hammingweight(Q9m9 & ~l.second)
					> hammingweight(Q9m9 & ~r.second);
			});

		// find m12 that still satisfies Q9 and Q8 conditions
		// forget m12 for which there is already another m12 with same Q8 and Q7
		// i.e. the difference is part of the Q9m9 tunnel

		std::map< std::array<uint32_t, 3>, std::array< uint32_t, 2 > > Q7810m1213;

		progress_display pdQ7810(m13Q10good.size());
		std::atomic_size_t m13readptr(0), m13writeptr(0);

		run_workload(threads, 
			[this,&m13Q10good,&pdQ7810,&m13readptr,&m13writeptr,&S,&Q7810m1213]
			(size_t ji, size_t jn)
			{
				fullstate_t myS = S;
				auto m12range = MA.word_range(12);
				std::map< std::array<uint32_t, 3>, std::array< uint32_t, 2 > > tmp;
				bool stop = false;
				while (!stop)
				{
					
					size_t m13ptr = m13readptr.fetch_add(1);
					if (m13ptr >= m13Q10good.size())
						break;

					auto m13Q10 = m13Q10good[m13ptr];
					myS.m[13] = m13Q10.first;
					myS.Qt(10) = m13Q10.second;

					// compute all results for m12range and store in tmp					
					for (uint32_t m12 : m12range)
					{
						myS.m[12] = m12;
						myS.m[10] = 0; // fixed dummy value

						myS.computeQtm3(12);
						if (!masked_value_QtQtp1(9,myS).check(myS.Qt(9)))
							continue;

						myS.computeQtm3(11);
						if (!masked_value_QtQtp1(8,myS).check(myS.Qt(8)))
							continue;
					
						myS.computeQtm3(10); // virtual Q7-under-fixed-dummy-m10

						if (!checkrotationQtQtp1(9,myS))
							continue;
						if (!checkrotationQtQtp1(8,myS))
							continue;

						tmp.insert( std::make_pair( std::array<uint32_t,3>({myS.Qt(7),myS.Qt(8),myS.Qt(10)}), std::array<uint32_t,2>({ myS.m[12], myS.m[13] })));
					}

					// wait till ready to write
					while (m13ptr != m13writeptr)
						std::this_thread::yield();
					// write results in tmp to Q7810m1213
					{
						lock_t lock(mut);
						++pdQ7810;
						for (auto& kv : tmp)
						{
							auto itf = Q7810m1213.insert(kv);
							if (itf.second && hammingweight(uint64_t(Q7810m1213.size()))==1) std::cout << Q7810m1213.size() << " " << std::flush;

						}
						if (Q7810m1213.size() >= 1ULL<<21)
							stop = true;
					}
					// finished writing, increase write ptr for next thread
					++m13writeptr;
					
					tmp.clear();
				}	
			});
		std::cout << "\nQ7810m1213size: " << Q7810m1213.size() << std::endl;
		

		// iterate over those m12, and find all (m10, m12) satisfying Q7-Q9, Q22-Q23
		static const size_t maxLUT = 1ULL<<29;
		
		progress_display pdm12m13(Q7810m1213.size());


		vecQ7m10m12m13.clear();
		vecQ7m10m12m13.reserve(maxLUT);
		
		auto m1213it = Q7810m1213.begin();
		run_workload(threads, 
			[this,&S,&Q7810m1213,&m1213it,&pdm12m13]
			(size_t ji, size_t jn)
			{
				vector< std::array<uint32_t,4> > tmp;
				tmp.reserve(1ULL<<20);
				fullstate_t myS = S;
				while (true)
				{
					auto it = Q7810m1213.begin();
					// grab next m1213it
					{
						lock_t lock(mut);
						it = m1213it;
						++m1213it;
						++pdm12m13;
						if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
							break;
					}
					if (it == Q7810m1213.end())
						break;

					myS.m[12] = it->second[0];
					myS.m[13] = it->second[1];
					myS.computeQtm3(13);
					myS.computeQtm3(12);
					myS.computeQtm3(11);
					for (uint32_t m10 : good_m10)
					{
						myS.m[10] = m10;

						myS.computeQtm3(10); // m10
						if (!masked_value_QtQtp1(7,myS).check(myS.Qt(7)))
							continue;

						myS.computeQtp1(21); // m10
						if (!masked_value_QtQtm1(22,myS).check(myS.Qt(22)))
							continue;

						myS.computeQtp1(22); // m15
						if (!masked_value_QtQtm1(23,myS).check(myS.Qt(23)))
							continue;

						myS.computeQtp1(23); // m4
						if (!masked_value_QtQtm1(24,myS).check(myS.Qt(24)))
							continue;

						if (!checkrotationQtQtp1(7,myS))
							continue;
						if (!checkrotationQtQtp1(21,myS))
							continue;
						if (!checkrotationQtQtp1(22,myS))
							continue;
						if (!checkrotationQtQtp1(23,myS))
							continue;
					
						tmp.push_back(std::array<uint32_t,4>({myS.Qt(7), myS.m[10], myS.m[12], myS.m[13]}));
						if (tmp.size() >= (1ULL<<20))
						{
							lock_t lock(mut);
							for (auto& v : tmp)
							{
								if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
									break;
								vecQ7m10m12m13.emplace_back(v);
							}
							tmp.clear();
						}
					}
				}
				lock_t lock(mut);
				for (auto& v : tmp)
				{
					if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
						break;
					vecQ7m10m12m13.emplace_back(v);
				}
			});
		std::cout << "Converting into look-up table..." << std::flush;

		std::sort(vecQ7m10m12m13.begin(), vecQ7m10m12m13.end(), 
			[](const std::array<uint32_t,4>& l, const std::array<uint32_t,4>& r)
			{
				return l[0] < r[0];
			});
		for (size_t i = 0; i < vecQ7m10m12m13.size(); )
		{
			Q7ptr[ vecQ7m10m12m13[i][0] ] = i;
			size_t j = i+1;
			while (j < vecQ7m10m12m13.size() && vecQ7m10m12m13[j][0] == vecQ7m10m12m13[i][0])
				++j;
			i = j;
		}
		std::cout << "done." << std::endl;
		std::cout << "vecQ7m10m12m13size: " << vecQ7m10m12m13.size() << " Q7ptrsize: " << Q7ptr.size() << std::endl;





		/* MAIN LOOP */
		S.Qt(-3) = ihv1[0]; S.Qt(-2) = ihv1[3]; S.Qt(-1) = ihv1[2]; S.Qt(0) = ihv1[1];
		
		S.computeQtp1(0); // m0 is fixed
		if (!masked_value_QtQtm1(1,S).check(S.Qt(1)))
			throw std::runtime_error("Q1-inconsistency");
		if (!checkrotationQtQtp1(0, S))
			throw std::runtime_error("Q1-inconsistency");
		S.computeQtp1(1); // m1 is fixed
		if (!masked_value_QtQtm1(2,S).check(S.Qt(2)))
			throw std::runtime_error("Q2-inconsistency");
		if (!checkrotationQtQtp1(1, S))
			throw std::runtime_error("Q2-inconsistency");
		
		counter_exponential_print m2m3cnt("m2m3cnt"), Q7attempts("Q7attempts"), Q7match("Q7match"), Q7success("Q7success"), 
					  m7ok("m7ok"), m8ok("m8ok"), m9ok("m9ok"), m12ok("m12ok"), Q24ok("Q24ok");

		std::vector<uint32_t> m2rndrange;
		for (uint32_t m2 : MA.word_range(2))
			m2rndrange.emplace_back(m2);
		randomize_vector(m2rndrange);
	
		std::atomic_size_t m2index(0);

		fullstate_t Scopy = S;

		run_workload(threads, 
		[this,&Scopy,&m2index,&m2rndrange,&Q24ok]
		(size_t ji, size_t jn)
		{
		   auto S = Scopy;
		   auto m3range = MA.word_range(3);
		   while (true)
		   {
			size_t m2i = m2index.fetch_add(1);
			if (m2i >= m2rndrange.size())
				break;

			S.m[2] = m2rndrange[m2i];

			S.computeQtp1(2);
			if (!masked_value_QtQtm1(3,S).check(S.Qt(3)))
				continue;
			if (!checkrotationQtQtp1(2,S))
				continue;

			for (uint32_t m3 : m3range)
			{
//				++m2m3cnt;

				S.m[3] = m3;
				S.computeQtp1(3);
				if (!masked_value_QtQtm1(4,S).check(S.Qt(4)))
					continue;
				if (!checkrotationQtQtp1(3,S))
					continue;

				S.computeQtp1(4); // m4 fixed
				if (!masked_value_QtQtm1(5,S).check(S.Qt(5)))
					continue;

				S.computeQtp1(5); // m5 fixed
				if (!masked_value_QtQtm1(6,S).check(S.Qt(6)))
					continue;

				S.computeQtp1(6); // m6 fixed
				if (!masked_value_QtQtm1(7,S).check(S.Qt(7)))
					continue;

				if (!checkrotationQtQtp1(4,S))
					continue;
				if (!checkrotationQtQtp1(5,S))
					continue;
				if (!checkrotationQtQtp1(6,S))
					continue;

//				++Q7attempts;
				auto it = Q7ptr.find(S.Qt(7));
				if (it == Q7ptr.end())
					continue;

				for (size_t q7i = it->second; q7i < vecQ7m10m12m13.size() && vecQ7m10m12m13[q7i][0] == S.Qt(7); ++q7i)
				{
//					++Q7match;

					S.m[10] = vecQ7m10m12m13[q7i][1];
					S.m[12] = vecQ7m10m12m13[q7i][2];
					S.m[13] = vecQ7m10m12m13[q7i][3];

					S.computeQtm3(13); // Q10
					S.computeQtm3(12); // Q9
					S.computeQtm3(11); // Q8
					S.computeWt(7);
					if (!MA.checkword(7, S.m[7]))
						continue;
					if (!checkrotationQtQtp1(7,S))
						continue;

//					++m7ok;

					// now iterate over Q9m9 tunnel
					S.computeQtp1(21); // m10 => Q22 already checked
					S.computeQtp1(22); // m15 => Q23 already checked
					S.computeQtp1(23); // m4 => Q24 already checked

					uint32_t Q9m9tunnel = ~Qvaluemask[offset+9] & ~Qprev[offset+10] & S.Qt(11) & ~S.Qt(10);
					uint32_t Q9org = S.Qt(9);
					uint32_t Q9cur = 0;
					do {
						Q9cur -= 1; Q9cur &= Q9m9tunnel;
						S.Qt(9) = Q9org ^ Q9cur;
						
						S.computeWt(9);
						if (!MA.checkword(9, S.m[9]))
							continue;
//						++m9ok;
						S.computeWt(8);
						if (!MA.checkword(8, S.m[8]))
							continue;
//						++m8ok;
						S.computeWt(12);
						if (!MA.checkword(12, S.m[12]))
							continue;
//						++m12ok;

						S.computeWt(10);
						if (!MA.checkword(10, S.m[10]))
							throw std::runtime_error("Q9m9tunnel: m10-invalid");
						S.computeWt(11);
						if (!MA.checkword(11, S.m[11]))
							throw std::runtime_error("Q9m9tunnel: m11-invalid");

/*
						// Not checking the specific conditions for Q25 and Q26 increases the success probability by about 10%,
						// as there exist other differential paths for these steps that are also possible.
						// The part of the code is only executed a number of times that is in the order of a few million times.
						// Hence, the total cost of calling the 'heavy' check_solution each time should be less than a second anyway.
						// Thus not verifying Q25 and Q26 has negligble performance impact, and a noticable success probability increase.
						int tQtok = 24;
						for (int t = 24; t < 26; ++t)
						{
							S.computeQtp1(t);
							if (!masked_value_QtQtm1(t+1,S).check(S.Qt(t+1)))
								break;
							if (!checkrotationQtQtp1(t,S))
								break;
							tQtok = t+1;
						}
						if (tQtok < 24)
							continue;
*/
						++Q24ok;
						
						check_solution(S);
						
					} while (Q9cur != 0);
				}
			}
		   }
		}); // run_workload
	
	} // m4-loop
}

void textcoll_solver_t::check_solution(const fullstate_t& sol)
{
	fullstate_t S = sol;
#if 0
	// check consistency of solution
	for (int t = 0; t <= 23; ++t)
	{
		int wt = md5_wt[t];
		S.computeWt(t);
		if (S.m[wt] != sol.m[wt])
			std::cerr << "sol:W" << t << "-inconsistency!" << std::flush;
	}
#endif

	uint32 block1[16], block2[16];
	for (int k = 0; k < 16; ++k)
	{
		block1[k] = sol.m[k];
		block2[k] = block1[k] + m_diff[k];
	}
	uint32 solihv1[4], solihv2[4];
	for (unsigned i = 0; i < 4; ++i)
	{
		solihv1[i] = solihv2[i] = ihv1[i];
	}
	md5compress(solihv1, block1);
	md5compress(solihv2, block2);
	
	if (solihv2[0]-solihv1[0] != dQ[0] + dQ[offset+61])
		return;
		
	std::cout << "\n === Partial dIHV solution ( each has prob 1/32 to be full solution ) === " << std::endl;
	std::cout << "MSG: ";
	for (unsigned t = 0; t < 16; ++t)
	{
		for (unsigned b = 0; b < 4; ++b)
			std::cout << char((block1[t]>>(8*b))&0xFF);
	}
	std::cout << std::endl;
	
	if (solihv2[1]-solihv1[1] != dQ[3] + dQ[offset+64])
		return;
	if (solihv2[2]-solihv1[2] != dQ[2] + dQ[offset+63])
		return;
	if (solihv2[3]-solihv1[3] != dQ[1] + dQ[offset+62])
		return;

	// necessary conditions for second block of the attack
	if ( (solihv1[1]>>31) != (solihv1[2]>>31) || (solihv1[2]>>31) != (solihv1[3]>>31) )
		return;
		
	std::cout << " ====== FULL SOLUTION FOUND =====" << std::endl;
        uint32 x = xrng128();
        std::string filename1 = "textcoll1_block1_" + boost::lexical_cast<string>(x) + ".txt";
        std::string filename2 = "textcoll2_block1_" + boost::lexical_cast<string>(x) + ".txt";
        std::cout << "Saving to " << filename1 << " & " << filename2 << "." << std::endl;
        ofstream of1(filename1.c_str(), ios::binary);
        ofstream of2(filename2.c_str(), ios::binary);
        if ((!of1) || (!of2)) {
                cerr << "Cannot open output files!" << endl;
        }
        save_block(of1, block1);
        save_block(of2, block2);
        of1.close();
        of2.close();
        exit(0);
}

void textcoll_solver_t::compute_good_m10(const fullstate_t& Sin)
{
	fullstate_t S = Sin;
	good_m10.clear();

	auto m10range = MA.word_range(10);
	for (uint32_t m10 : m10range)
	{
		S.m[10] = m10;

		S.computeQtp1(21); // m10
		if (!masked_value_QtQtm1(22,S).check(S.Qt(22)))
			continue;

		S.computeQtp1(22); // m15
		if (!masked_value_QtQtm1(23,S).check(S.Qt(23)))
			continue;

		S.computeQtp1(23); // m4
		if (!masked_value_QtQtm1(24,S).check(S.Qt(24)))
			continue;

		if (!checkrotationQtQtp1(21,S))
			continue;
		if (!checkrotationQtQtp1(22,S))
			continue;
		if (!checkrotationQtQtp1(23,S))
			continue;

		good_m10.emplace_back(m10);
	}
	
}
