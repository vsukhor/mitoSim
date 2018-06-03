#ifndef Oel_hpp
#define Oel_hpp

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <stdarg.h>

namespace Utils
{


//template<typename Q, typename Enabler = void>
//class Oel
//{};

// specialization for fundamental scalars or strings xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

class Oel {

public:

	std::ostream&  so {std::cout};
	std::ofstream& sl;
	
	Oel( std::ofstream& sl,
		 int precision = 6 );

	Oel( bool ou,
		 bool lo,
		 std::ofstream& sl,
		 int precision = 6 );
	
	void set_formats( int precision ) noexcept;

	template<typename V, std::size_t N>
	void print( const std::string& name, const std::array<V,N>& v ) noexcept;

	template<typename V>
	void print( const std::string& name, const std::vector<V>& v ) noexcept;

	template<bool end=true>
	void print( const std::string& ) noexcept;
	void print0( const std::string& ) noexcept;

	template<bool end=true>
	void print( const char *fmt, ... ) noexcept;
		
	void exit( const std::string& s ) noexcept;
	void exit( const char *fmt, ... ) noexcept;
	

private:

	char buf [1024];
	
	bool ou {true};
	bool lo {true};
	
	template<typename IO, typename V>
	void prn( bool b, IO& io, const V& v, bool end ) noexcept;
	template<typename IO, typename V> 
	void prn( bool b, IO& io, const std::string& name, const V& v, bool end ) noexcept;
	template<typename IO, typename V> 
	void prn( bool b, IO& io, const std::string& name, const std::vector<V>& v, bool end ) noexcept;

};

template<typename IO, typename V> 
void Oel::prn( bool b, IO& io, const V& v, bool end ) noexcept
{ 
	if (b) {
		io << v << " "; 
		if(end) io << std::endl;
	}
}

template<typename IO, typename V> 
void Oel::prn( bool b, IO& io, const std::string& name, const V& v, bool end ) noexcept
{ 
	if (b) {
		io << name << " " << v << " "; 
		if (end) io << std::endl;
	}
}

template<typename V, std::size_t N>
void Oel::print( const std::string& name, const std::array<V,N>& v ) noexcept
{
	const std::string emp("");

	print<false>(emp+name);
	for (const auto o : v)
		print<false>(emp.c_str(), o);
	print("\n");
}

template<typename V>
void Oel::print( const std::string& name, const std::vector<V>& v ) noexcept
{
	const std::string emp("");

	print<false>(emp+name);
	for (const auto o : v)
		print<false>(emp.c_str(), o);
	print("\n");
}
	
template<bool end>
void Oel::print( const std::string& s ) noexcept
{	
	prn( lo, sl, s, end );
	prn( ou, so, s, end );
}

template<bool end>
void Oel::print( const char *fmt, ... ) noexcept
{
	va_list va;
	va_start(va, fmt);
	const auto n = vsprintf(buf, fmt, va);
	va_end(va);
	const auto s = std::string(buf).substr(0, n);
	prn(lo, sl, s, end);
	prn(ou, so, s, end);
}

}

#endif /* Oel_hpp */
