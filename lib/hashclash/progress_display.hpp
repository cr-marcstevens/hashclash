// extension from boost::progress_display:
//	Copyright Beman Dawes 1994-99.  
//      Distributed under the Boost Software License, Version 1.0. 
//      (See http://www.boost.org/LICENSE_1_0.txt)

#ifndef PROGRESS_DISPLAY_HPP
#define PROGRESS_DISPLAY_HPP

#include <iostream>
#include <string>

#include <boost/cstdint.hpp>
#include <boost/utility.hpp>

namespace hashclash {

	class progress_display : private boost::noncopyable	{
	public:
		typedef boost::uint64_t uint64;

		explicit progress_display(uint64 expected_count, 
										bool show_scale = true,
										std::ostream & os = std::cout,
										const std::string & s1 = "\n",
										const std::string & s2 = "",
										const std::string & s3 = "" 
								) 
			: m_show_scale(show_scale), m_os(os), m_s1(s1), m_s2(s2), m_s3(s3) 
		{ 
			restart(expected_count); 
		}

		void restart(uint64 expected_count)
		{
			_count = _next_tic_count = _tic = 0;
			_expected_count = expected_count;
			if ( !_expected_count ) _expected_count = 1;

			if (m_show_scale)
				m_os	<< m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
						<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|\n"
						<< m_s3 << std::flush;
		}

		void redraw(bool show_scale = true)
		{
			_next_tic_count = _tic = 0;
			if (show_scale)
				m_os	<< m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
						<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|\n"
						<< m_s3 << std::flush;
			display_tic();
		}

		uint64 operator+=(uint64 increment )
		{
			if ( (_count += increment) >= _next_tic_count ) { display_tic(); }
			return _count;
		}

		uint64 operator++()           { return operator+=( 1 ); }
		uint64 count() const          { return _count; }
		uint64 expected_count() const { return _expected_count; }

	private:
		bool               m_show_scale;
		std::ostream &     m_os;
		const std::string  m_s1;
		const std::string  m_s2;
		const std::string  m_s3;

		uint64 _count, _expected_count, _next_tic_count;
		unsigned int  _tic;
		void display_tic()
		{
			unsigned int tics_needed = static_cast<unsigned int>
				( (static_cast<double>(_count)/static_cast<double>(_expected_count))*50.0 );

			for (; _tic <= tics_needed; ++_tic)
				m_os << '*' << std::flush;
			if (_count == _expected_count)
			{
				for (; _tic <= 50; ++_tic)
				m_os << '*' << std::flush;
				m_os << std::endl;
			}

			if (_tic != 49)
				_next_tic_count = static_cast<uint64>
				( (static_cast<double>(_tic+1)/50.0)*static_cast<double>(_expected_count) );
			else // to avoid a precision error at the end
				_next_tic_count = _expected_count;
		}
	};

} // namespace hashclash

#endif // PROGRESS_DISPLAY_HPP
