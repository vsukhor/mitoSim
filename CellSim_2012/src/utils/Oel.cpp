
#include "Oel.h" 

namespace Utils
{

Oel::Oel( std::ofstream& sl,
		  int precision
		)
	: sl(sl)
{
	set_formats(precision);
}

Oel::Oel( bool ou,
		  bool lo,
		  std::ofstream& sl,
		  int precision
		)
	: sl {sl}
	, ou {ou}
	, lo {lo}
{
	set_formats( precision );
}

void Oel::set_formats( int precision ) noexcept
{
	so.precision(precision);
	so.setf(std::ios::scientific);
	sl.precision(precision);
	sl.setf(std::ios::scientific);
}

void Oel::print0( const std::string& s ) noexcept
{
	print<false>(s);
}
void Oel::exit( const std::string& s ) noexcept
{	
	print<true>(s);
	::exit(EXIT_FAILURE);
}

void Oel::exit( const char *fmt, ... ) noexcept
{	
	va_list va;
	va_start(va, fmt);
	const auto n = vsprintf(buf, fmt, va);
	va_end(va);
	const auto s = std::string(buf).substr(0, std::size_t(n));
	prn(lo, sl, s, 1);
	prn(ou, so, s, 1);

	::exit(EXIT_FAILURE);
}

}
