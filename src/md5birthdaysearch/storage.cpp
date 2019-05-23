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
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include "main.hpp"
#include "storage.hpp"

boost::mutex storage_mutex;
#define LOCK_STORAGE_MUTEX	boost::mutex::scoped_lock lock(storage_mutex);

using namespace hashclash;
using namespace std;

// service for CUDA part
void boost_thread_yield()
{
	boost::this_thread::sleep(boost::posix_time::milliseconds(1));
}

void storage_type::reserve_memory(uint64 m)
{
	LOCK_STORAGE_MUTEX;
	trail_hash.resize(m >> 5);
	for (unsigned i = 0; i < trail_hash.size(); ++i)
		trail_hash[i].reserve(2 << 5);
}

void storage_type::insert_trail(const trail_type& tr)
{
	LOCK_STORAGE_MUTEX;
	uint32 hashindex = (tr.end[1] / procmodn) % trail_hash.size();
	vector<trail_type>::const_iterator 
		cit = trail_hash[hashindex].begin(),
		citend = trail_hash[hashindex].end();

	for (; cit != citend; ++cit)
	{
		if (cit->end[0] == tr.end[0]
			&& cit->end[1] == tr.end[1]
			&& cit->end[2] == tr.end[2])
		{
			if (tr != *cit) {
				++totcoll;
				collisions.push_back( make_pair(tr, *cit) );				
			}
		}
	}
	if (!memhardlimit || trail_hash[hashindex].size() < trail_hash[hashindex].capacity()) {
		trail_hash[hashindex].push_back(tr);
	} else {
		unsigned i = (++tmpptr) % trail_hash[hashindex].size();
		trail_hash[hashindex][i] = tr;
	}
}

void storage_type::insert_trails(const vector<trail_type>& trs)
{
	LOCK_STORAGE_MUTEX;
	vector<trail_type>::const_iterator cit, citend;
	vector<trail_type>::const_iterator 
		trsit = trs.begin(),
		trsend = trs.end();
	for (; trsit != trsend; ++trsit) {
		uint32 hashindex = (trsit->end[1] / procmodn) % trail_hash.size();
		cit = trail_hash[hashindex].begin();
		citend = trail_hash[hashindex].end();
		for (; cit != citend; ++cit)
		{
			if (cit->end[0] == trsit->end[0]
				&& cit->end[1] == trsit->end[1]
				&& cit->end[2] == trsit->end[2])
			{
				if (*trsit != *cit)
				{
					++totcoll;
					collisions.push_back( make_pair(*trsit, *cit) );
				}
			}
		}
		if (!memhardlimit || trail_hash[hashindex].size() < trail_hash[hashindex].capacity()) {
			trail_hash[hashindex].push_back(*trsit);
		} else {
			unsigned i = (++tmpptr) % trail_hash[hashindex].size();
			trail_hash[hashindex][i] = *trsit;
		}
	}
}

void storage_type::get_birthdaycollisions(vector< pair<trail_type,trail_type> >& coll, unsigned multiple_of)
{
	LOCK_STORAGE_MUTEX;
	coll.clear();
	if (collisions.size() >= multiple_of) {
		size_t count = (collisions.size()/multiple_of)*multiple_of;
		coll.resize(count);
		std::copy(collisions.end() - count, collisions.end(), coll.begin());
		collisions.resize(collisions.size()-count);
	}
}

unsigned storage_type::get_collqueuesize()
{
	unsigned k;
	{
		LOCK_STORAGE_MUTEX;
		k = unsigned(collisions.size());
	}
	return k;
}

unsigned storage_type::get_totcoll()
{
	unsigned result;
	{
		LOCK_STORAGE_MUTEX;
		result = totcoll;
	}
	return result;
}

storage_type main_storage;
