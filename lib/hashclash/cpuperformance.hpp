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

#ifndef HASHCLASH_CPUPERFORMANCE_HPP
#define HASHCLASH_CPUPERFORMANCE_HPP


#include <vector>
#include <string>
#include <iostream>
#include <cstdio>

#ifndef __GNUC__
#include <intrin.h>
#endif

#include "types.hpp"

namespace hashclash {

	inline uint64 cpu_timestamp()
	{
#ifdef __GNUC__
		uint32 highpart, lowpart;
		asm volatile("rdtsc": "=d"(highpart), "=a"(lowpart));
		return (uint64(highpart) << 32) | uint64(lowpart);
#else
		return __rdtsc();
#endif
		
	}

	inline void start_update_counter(uint64& performancecounter)
	{
		performancecounter -= cpu_timestamp();
	}
	inline void end_update_counter(uint64& performancecounter)
	{
		performancecounter += cpu_timestamp();
	}

	class update_performance_counter {
		uint64& _counter;
	public:
		update_performance_counter(uint64& performance_counter)
			: _counter(performance_counter)
		{
			_counter -= cpu_timestamp();
		}
		~update_performance_counter()
		{
			_counter += cpu_timestamp();
		}
	};

	class performance_counter_manager {
		std::vector<uint64*> counters;
		std::vector<std::string> descriptions;
	public:
		void add_performance_counter(uint64& counter, const std::string& description)
		{
			counters.push_back(&counter);
			descriptions.push_back(description);
		}

		void show_results()
		{
			for (unsigned i = 0; i < counters.size(); ++i)
			{
				std::cout << "Counter " << i << ": \t" << descriptions[i] << std::endl;
				std::cout << "Counter " << i << ": \tValue = " << (*counters[i]) << std::endl;
			}
		}
	};
} // namespace hashclash


#endif // HASHCLASH_CPUPERFORMANCE_HPP
