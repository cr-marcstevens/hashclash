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

#ifndef STORAGE_HPP
#define STORAGE_HPP

// to be extended to support hard drive storage

class storage_type {
public:
	storage_type(): totcoll(0), tmpptr(0) {}
	~storage_type() {}

	void set_parameters(const birthday_parameters& parameters)
	{
		procmodn = parameters.modn;
		memhardlimit = parameters.memhardlimit;
	}

	void reserve_memory(uint64 m);
	void insert_trail(const trail_type& tr);
	void insert_trails(const vector<trail_type>& trs);
	void get_birthdaycollisions(vector< pair<trail_type,trail_type> >& collisions, unsigned multiple_of = 1);
	unsigned get_collqueuesize();
	unsigned get_totcoll();
private:
	vector< vector<trail_type> > trail_hash;
	vector< pair<trail_type,trail_type> > collisions;
	unsigned procmodn;
	bool memhardlimit;
	uint32 tmpptr;

	unsigned totcoll;
};

extern storage_type main_storage;

#endif // STORAGE_HPP
