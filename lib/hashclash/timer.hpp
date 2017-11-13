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

#ifndef HASHCLASH_TIMER_HPP
#define HASHCLASH_TIMER_HPP

#include "types.hpp"

namespace hashclash
{

	class timer_detail;
	class timer {
	public:
		timer(bool direct_start = false);
		~timer();
		void start();
		void stop();
		double time() const;// get time between start and stop (or now if still running) in seconds
		bool isrunning() const { return running; } // check if timer is running

	private:
		timer_detail* detail;
		bool running;
	};

} // namespace hashclash

#endif // HASHCLASH_TIMER_HPP
