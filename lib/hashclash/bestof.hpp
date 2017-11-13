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

#ifndef HASHCLASH_BESTOF_HPP
#define HASHCLASH_BESTOF_HPP

#include <map>
#include <vector>
#include <utility>

#include "types.hpp"

namespace hashclash {

	template<typename _value_type, typename _count_type = uint64>
	class bestof {
	public:
		typedef _value_type value_type;
		typedef _count_type count_type;
		typedef bestof<value_type, count_type> my_type;
		typedef std::map<value_type, count_type> map_type;
		typedef std::vector< std::pair<value_type, count_type> > vector_type;
		typedef typename vector_type::const_iterator const_iterator;

		count_type& operator[](const value_type& v) { return itemtocount[v]; }
		const count_type& operator[](const value_type& v) const { return itemtocount[v]; }

		typename map_type::size_type size() const {
			return itemtocount.size();
		}
		void clear() {
			itemtocount.clear();
			counttoitem.clear();
		}
		void invert() {
			counttoitem.reserve(size());
			for (typename map_type::const_iterator cit = itemtocount.begin(); cit != itemtocount.end(); ++cit)
				counttoitem.push_back( *cit );
			sort(counttoitem.begin(), counttoitem.end(), bestof_less());
		}
		const_iterator begin() {
			if (counttoitem.size() != size()) invert();
			return counttoitem.begin();
		}
		const_iterator end() {
			if (counttoitem.size() != size()) invert();
			return counttoitem.end();
		}

		map_type itemtocount;
		vector_type counttoitem;

		struct bestof_less
			: public std::binary_function<std::pair<value_type, count_type>, std::pair<value_type, count_type>, bool>
		{
			bool operator()(const std::pair<value_type, count_type>& _Left, const std::pair<value_type, count_type>& _Right) const
			{
				return _Left.second > _Right.second;
			}
		};

	};

} // namespace hashclash

#endif //HASHCLASH_BESTOF_HPP
