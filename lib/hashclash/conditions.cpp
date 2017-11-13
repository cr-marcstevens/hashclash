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

#include <boost/utility.hpp>

#include "conditions.hpp"

namespace hashclash {

	std::ostream& operator<<(std::ostream& o, const bitcondition& v)
	{
		switch (v) {
			case bc_constant:	o << '.'; break;
			case bc_plus:		o << '+'; break;
			case bc_minus:		o << '-'; break;
			case bc_zero:		o << '0'; break;
			case bc_one:		o << '1'; break;
			case bc_prev:		o << '^'; break;
			case bc_prevn:		o << '!'; break;
			case bc_prev2:		o << 'm'; break;
			case bc_prev2n:		o << '#'; break;
			case bc_or2:		o << '?'; break;
			case bc_next:		o << 'V'; break;
			case bc_nextn:		o << 'Y'; break;
			case bc_next2:		o << 'W'; break;
			case bc_next2n:		o << 'H'; break;
			case bc_or2b:		o << 'Q'; break;
			default:
				throw;
		}
		return o;
	}

	std::istream& operator>>(std::istream& i, bitcondition& v)
	{
		char c;
		if (!(i >> c)) return i;
		switch (c) {
			case '.':	v = bc_constant; break;
			case '+':	v = bc_plus; break;
			case '-':	v = bc_minus; break;
			case '0':	v = bc_zero; break;
			case '1':	v = bc_one; break;
			case '^':	v = bc_prev; break;
			case '!':	v = bc_prevn; break;
			case 'm':	v = bc_prev2; break;
			case '#':	v = bc_prev2n; break;
			case '?':	v = bc_or2; break;
			case 'V':	v = bc_next; break;
			case 'Y':	v = bc_nextn; break;
			case 'W':	v = bc_next2; break;
			case 'H':	v = bc_next2n; break;
			case 'Q':	v = bc_or2b; break;
			default:
				i.putback(c);
				i.setstate(std::ios::failbit);
		}
		return i;
	}

	std::ostream& operator<<(std::ostream& o, const byteconditions& bc)
	{
		for (unsigned b = 8; b > 0; --b)
			o << bc[b-1];
		return o;
	}

	std::istream& operator>>(std::istream& i, byteconditions& bc)
	{
		unsigned b = 8;
		while ((i) && b > 0)
		{
			bitcondition tmp = bc_constant;
			i >> tmp;
			if (i) bc.set(b-1, tmp);
			--b;
		}
		return i;
	}

	std::ostream& operator<<(std::ostream& o, const wordconditions& w)
	{
		return o << '|' << w.bytes[3] << ' ' << w.bytes[2] << ' ' << w.bytes[1] << ' ' << w.bytes[0] << '|';
	}

	std::istream& operator>>(std::istream& i, wordconditions& w)
	{
		char c;
		if (!(i >> c)) return i;
		if (c != '|') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		if (!(i >> w.bytes[3] >> c)) return i;
		if (c != ' ') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		if (!(i >> w.bytes[2] >> c)) return i;
		if (c != ' ') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		if (!(i >> w.bytes[1] >> c)) return i;
		if (c != ' ') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		if (!(i >> w.bytes[0] >> c)) return i;
		if (c != '|') {
			i.putback(c);
			i.setstate(std::ios::failbit);
			return i;
		}
		return i;
	}

} // namespace hashclash
