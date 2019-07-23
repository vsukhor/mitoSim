/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann.
   All Rights Reserved.
   Developed at Helmholtz Center for Infection Research, Braunschweig, Germany.
   Please see Readme file for further information

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

============================================================================== */

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
