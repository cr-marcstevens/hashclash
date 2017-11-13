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

#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <vector>
#include <string>

#include <boost/filesystem/operations.hpp>

#include <hashclash/types.hpp>
#include <hashclash/saveload_gz.hpp>
#include <hashclash/saveload_bz2.hpp>

using namespace hashclash;
using namespace std;

#include "birthday_types.hpp"

namespace boost {
	namespace serialization {
		template<class Archive>
		void serialize(Archive& ar, trail_type& t, const unsigned int file_version)
		{
			ar & boost::serialization::make_nvp("start0", t.start[0]);
			ar & boost::serialization::make_nvp("start1", t.start[1]);
			ar & boost::serialization::make_nvp("start2", t.start[2]);
			ar & boost::serialization::make_nvp("end0", t.end[0]);
			ar & boost::serialization::make_nvp("end1", t.end[1]);
			ar & boost::serialization::make_nvp("end2", t.end[2]);
			ar & boost::serialization::make_nvp("len", t.len);
		}
	}
}

extern std::string workdir;
extern uint32 colla1, collb1, collc1, colla2, collb2, collc2;

int dostep(birthday_parameters& parameters);
void birthday(birthday_parameters& parameters);
void determine_nrblocks_distribution(birthday_parameters& parameters);

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, uint32 block[]);

// cuda

#ifdef HASHCLASH_HAVE_CUDA
#define HAVE_CUDA
#endif

void cuda_device_query();
int get_num_cuda_devices();

#endif // MAIN_HPP
