#ifndef Utils_h
#define Utils_h

#include <cmath>
#include <limits>
#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <numeric>
#include <memory>
#include <stdio.h>
#include <sys/stat.h>

namespace Utils {

#define __LINUX__

extern const std::string SLASH;

#ifdef FP32
	typedef float real;
#else
	typedef double real;
#endif

#define STR(x) std::to_string(x)

typedef unsigned long ulong;
typedef unsigned int uint;
typedef std::size_t szt;

template<typename T> using vec3 = std::vector<std::vector<std::vector<T>>>;
template<typename T> using vec2 = std::vector<std::vector<T>>;
template<typename T> using vup = std::vector<std::unique_ptr<T>>;	

template<typename T> constexpr const T zero = static_cast<T>(0.L);
template<typename T> constexpr const T half = static_cast<T>(.5L);
template<typename T> constexpr const T one = static_cast<T>(1.L);
template<typename T> constexpr const T two = static_cast<T>(2.L);

template<typename T> constexpr const T pi = static_cast<T>(3.1415926535897932384626433832795L);
template<typename T> constexpr const T twopi = two<T>*pi<T>;

template<typename T, typename Enabler = void> constexpr const T huge;
template<typename T> constexpr const T huge<T,std::enable_if_t<std::is_fundamental<T>::value>> =
	std::numeric_limits<T>::has_infinity ?
   	std::numeric_limits<T>::infinity() :
   	std::numeric_limits<T>::max();

template<typename T1, typename T2>
inline vec2<T1> array_like( const vec2<T2>& as )
{
	vec2<T1> me( as.size() );
	for (szt i=0; i<as.size(); i++)
		me[i].resize( as[i].size() );
	return me;
}

// returns how many elements in b /= 0 putting their indices to j
template<typename T> 
szt find( const std::vector<T>& b, 
		  std::vector<szt>& j )	noexcept
{	j.clear();
	for (szt i=0; i<b.size(); i++) 
		if (b[i] != T(0) ) 
			j.push_back( i );
	return j.size();
}

std::string padZeros3( szt n );

bool file_exists( const std::string& name );
bool fileExists( const std::string& name );
std::string get_current_time();

}



#endif /* Utils_h */
