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

#ifndef Errors_h
#define Errors_h

#include <iostream>
#include <string>
#include <exception>
#include <algorithm>

#include "Oel.h"

#ifdef _DEBUG
	#define XASSERT(EX, msg) \
		(void)( (EX) || ( Utils::assert_fun( #EX, __FILE__, __LINE__, msg ), 0 ) )
#else
	#define XASSERT(EX, msg)
#endif

namespace Utils 
{

long long assert_fun( const char* EX,
				 const char *file, 
				 int line, 
				 const std::string& msg );

int simple_error( const std::string& s );

class MyException: public std::exception {
public:
	MyException( const std::string& s, Oel& oel ) 
	{
		oel.print<true>(s);
	}

	MyException( const std::string& s ) 
	{
		std::cerr << s << std::endl;
	}
	
	template<typename T>
	MyException( const std::string& s, const T& u, Oel& oel ) 
	{
		oel.print(s, u);
	}
	
	static int throwIt( const std::string& s, Oel& oel )
	{
		throw MyException(s, oel);
		return EXIT_FAILURE;		// pro forma: to return anything
	}
	
	template<typename T>
	static int throwIt( const std::string& s, const T& u, Oel& oel )
	{
		throw MyException(s, u, oel);
		return EXIT_FAILURE;		// pro forma: to return anything
	}
};

// template for ParOutOfRangeException xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename Q, bool isDiscrete, typename Enabler = void>
class ParOutOfRangeException: public std::exception
{};

// specialization for fundamental scalars discrete xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename Q>
class ParOutOfRangeException<Q, true,		// discrete case
	std::enable_if_t<std::is_fundamental<Q>::value>> : public std::exception {

public:

	ParOutOfRangeException(const std::string& name, const Q& p, const std::vector<Q>& r)
		:	message{generate_message(name, p, r)}
	{
		std::cerr << message;
		std::abort();
	}

	ParOutOfRangeException(const std::string& name, const Q& p, const std::vector<Q>& r, Oel& oel)
		:	message{generate_message(name, p, r)}
	{
		oel.print(message);
	}

private:

	const std::string message;
	
	std::string generate_message(const std::string& name, const Q& p, const std::vector<Q>& r)
	{
		XASSERT(r.size() == 2, "Incorrect r size in ParOutOfRangeException");

		auto print = [](const std::vector<Q>& a) {
			std::string w {"{ "};
			for(const auto o : a)
				w += std::to_string(o)+" ";
			return w+"}";
		};
		return "Error in conf specification for parameter '"+name+
			   "' = "+std::to_string(p)+" :"+"\n\tthe value provided "+
			   " is outside thw acceptable range "+
			   print(r);
	}

};

// specialization for fundamental scalars continous xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename Q>
class ParOutOfRangeException<Q, false,		// continuous case
	std::enable_if_t<std::is_fundamental<Q>::value>> : public std::exception {

public:

	ParOutOfRangeException(const std::string& name, const Q& p, const std::vector<Q>& r)
		:	message{generate_message(name, p, r)}
	{
		std::cerr << message;
		std::abort();
	}

	ParOutOfRangeException(const std::string& name, const Q& p, const std::vector<Q>& r, Oel& oel)
		:	message{generate_message(name, p, r)}
	{
		oel.print(message);
	}

private:

	const std::string message;
	
	std::string generate_message(const std::string& name, const Q& p, const std::vector<Q>& r)
	{
		XASSERT(r.size() == 2, "Incorrect r size in ParOutOfRangeException");
		return "Error in conf specification for parameter '"+name+
			   "' = "+std::to_string(p)+" :"+"\n\tthe value provided "+
			   " is outside the acceptable range "+
			   "[ "+std::to_string(r[0])+", "+std::to_string(r[1])+" ]";
	}
};

}

#endif /* Errors_h */
