#ifndef Par_h
#define Par_h

#include <sstream>

#include "Misc.h"
#include "Errors.h"

namespace Utils {

class Oel;

template<typename Q>
class ParBase
{
	std::string name_;
	bool detect_name( std::ifstream&, 
					  std::string& ) const;
protected:
	
	bool isLoaded_ = false;

	ParBase() = default;
	ParBase( const std::string& name ) 
		: name_(name) 
	{}
	
	virtual void print( Oel& ) = 0;
	virtual void initialize( std::string ) = 0;
	std::string name() const noexcept;

public:
	void load( std::ifstream& );
	void load( const std::string& fname );
};

// if the line contains a valid parname-value combination, returns true and the value, otherwise retruns false
template<typename Q>
bool ParBase<Q>::detect_name( std::ifstream& config, 
							  std::string& value ) const
{	
	const std::string emp {" "};
	const std::string tab {"\t"};
	
	std::string line;
	getline(config, line);
	
	ulong commentpos = line.find_first_of('#');
	if (commentpos != std::string::npos)
		line.erase(commentpos);
	
	if (!line.length())
		return false;
	while (!line.substr(line.length()-1, 1).compare(emp) ||
		   !line.substr(line.length()-1, 1).compare(tab))
		line.erase(line.length()-1);
	if (!line.length())
		return false;
	
	int parnameend = -1;
	if (     line.find_first_of( emp ) == std::string::npos &&
			 line.find_first_of( tab ) != std::string::npos) parnameend = (int)line.find_first_of(tab);
	else if (line.find_first_of( emp ) != std::string::npos &&
			 line.find_first_of( tab ) == std::string::npos) parnameend = (int)line.find_first_of(emp);
	else if (line.find_first_of( emp ) != std::string::npos &&
			 line.find_first_of( tab ) != std::string::npos) parnameend = std::min( (int)line.find_first_of(emp),
																					 (int)line.find_first_of(tab) );
	const auto parname = line.substr(0, (size_t)parnameend);
	
	if (parname != name_)
		return false; 
	
	value = line.substr(line.find_last_of("=")+1);
	while( !value.substr(0, 1).compare(emp) ||
		   !value.substr(0, 1).compare(tab) )
		value.erase(value.begin());
	return true;
}

template<typename Q>
void ParBase<Q>::load( const std::string& fname ) 
{
	std::string parname, value;
	std::ifstream ifs(fname);
	if (!ifs.is_open()) {
		std::cout << "Unable to open config file at " << fname << std::endl;
		exit(0);
	}
	try {
		load(ifs);
	} catch(MyException& e) {
		return;
	}
}
template<typename Q>
void ParBase<Q>::load( std::ifstream& config )
{
	XASSERT(!isLoaded_, "Repeated load of " + name_ + "\n");
	
	config.clear();
	config.seekg(0, std::ios::beg);
	while (config.good()) {
		std::string value;
		if (!this->detect_name(config, value))
			continue;
		initialize(value);
		isLoaded_ = true;
		return;
	}
	throw MyException("Error: parrameter not loaded: " + name_);
}
template<typename Q>
std::string ParBase<Q>::name() const noexcept
{
	return name_;
} 

	
// template for Par xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename Q, bool isDiscrete, typename Enabler=void>
class Par : public ParBase<Q>
{};

// partial specialization for fundamental scalars or strings xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename T, bool isDiscrete>
class Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
												  std::is_same<T,std::string>::value>> : public ParBase<T>
{
	typedef T Q;
	using ParBase<T>::name;
	using ParBase<T>::isLoaded_;
	
	Q p_;
	
public:	
	Par( const std::string& name )
		: ParBase<T>::ParBase(name)
	{}
	Par( const std::string& name, 
		 const std::string& fname, 
		 const std::vector<T>& range )
		: ParBase<T>::ParBase(name)
	{
		this->load(fname);
		check_range(range);
//		print();
	}
	Par( const std::string& name,
		 const std::string& fname, 
		 Oel& oel,
		 const std::vector<T>& range )
		: ParBase<T>::ParBase(name)
	{
		this->load(fname);
		check_range(range, oel);
		print(oel);
	}
	
	void check_range( const std::vector<T>& r, Oel& oel )
	{
		if (!r.size()) return;		// use this case to omit string checkups

		if (isDiscrete) {
			if (std::find(r.begin(), r.end(), p_) == r.end())
				throw ParOutOfRangeException<T,isDiscrete>(name(), p_, r, oel);
		}
		else {
			XASSERT(r.size()==2, "size of r must be 2 for continuous parameters");
			if (p_<r[0] || p_>r[1])
				throw ParOutOfRangeException<T,isDiscrete>(name(), p_, r, oel);
		}
	}

	void check_range( const std::vector<T>& r )
	{
		if (!r.size()) return;		// use this case to omit string checkups

		if (isDiscrete) {
			if (std::find(r.begin(), r.end(), p_) == r.end())
				throw ParOutOfRangeException<T,isDiscrete>(name(), p_, r);
		}
		else {
			XASSERT(r.size()==2, "size of r must be 2 for continuous parameters");
			if (p_<r[0] || p_>r[1])
				throw ParOutOfRangeException<T,isDiscrete>(name(), p_, r);
		}
	}

	static auto readin(const std::string& s, const std::string& fname, Oel& oel, const std::vector<Q>& range )
	{ 
		return Par<Q,isDiscrete>({s, fname, oel, range})();
	};
	
	static auto readin(const std::string& s, const std::string& fname, const std::vector<Q>& range )
	{ 
		return Par<Q,isDiscrete>({s, fname, range})();
	};
	
	void print( Oel& oel ) final
	{
		oel.print(name()+" = "+STR(p_));
	}

	void print()
	{
		std::cout << name() << " = " << p_ << std::endl;
	}

	void set( const Q& val ) 
	{ 
		isLoaded_ = true;
		p_ = val; 
	}
	
	Q operator()() const 
	{ 
		XASSERT(isLoaded_, name());
		return p_; 
	}
	
private:

	void initialize( std::string value ) final
	{
		std::stringstream(value) >> p_;
	}

};	

// specialization for vectors of fundamental types xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename T, bool isDiscrete>
class Par<std::vector<T>, isDiscrete, std::enable_if_t<std::is_fundamental<T>::value>> : public ParBase<T>
{
	typedef std::vector<T> Q;

	using ParBase<T>::name;
	using ParBase<T>::isLoaded_;

	Q p_;
	szt expectedSize_;
	
public:	

	Par( const std::string& name, const szt& expectedSize )
		: ParBase<T>::ParBase(name)
		, expectedSize_(expectedSize)
	{}
	Par( const std::string& name,
		 const std::string& fname, 
		 Oel& oel,
		 std::vector<Q> range )
		: ParBase<Q>::ParBase(name)
	{
		this->load(fname);
		check_range(range, oel);
		print(oel);
	}
	
	void check_range( const std::vector<Q>& r, Oel& oel )
	{
		if (!r.size()) return;		// use this case to omit string checkups

		XASSERT(!isDiscrete || r.size()==2, "size of r must be 2 for continuous parameters");
		if (isDiscrete) {
			if (std::find(r.begin(), r.end(), p_) == r.end())
				throw ParOutOfRangeException<Q,isDiscrete>(name(), p_, r, oel);
		}
		else {
			if (p_<r[0] || p_>r[1])
				throw ParOutOfRangeException<Q,isDiscrete>(name(), p_, r, oel);
		}
	}

	static auto readin(const std::string& s, const std::string& fname, Oel& oel )
	{ 
		return Par<Q,isDiscrete>({s, fname, oel})();
	};
	
	void print( Oel& oel ) final
	{
		oel.print( name, p_, 1 );
	}
	
	Q operator()() const 
	{ 
		XASSERT(true, name); 
		return p_; 
	}
	
	T operator[]( const szt& i ) const 
	{
		XASSERT(isLoaded_, name()); 
		XASSERT(i<p_.size(), name()); 
		return p_[i];
	}

private:
	
	void initialize( std::string value ) final
	{
		const std::string emp {" "};
		const std::string tab {"\t"};
		const MyException improperSizeEx {"Improper Config::" + name + " initialization: Excessive data size"};
		while (value.length()) {
			ulong e {value.find(emp)};
			if (e == std::string::npos) e = value.find(tab);
			if (e == std::string::npos) e = value.length();
			const std::string val {value.substr(0, e)};
			if (val.length() < 1)
				throw MyException("Error in config file: Number of elelments in " + name() +
								  " is " + STR(p_.size()) + " which is insufficient");
			T tmp;
			std::stringstream(val) >> tmp;
			p_.push_back(tmp);
			value.erase(0, e);
			while (!value.substr(0, 1).compare(emp) ||
				   !value.substr(0, 1).compare(tab))
				value.erase(value.begin());
			if (p_.size() > expectedSize_)
				throw improperSizeEx;
		}
		if (p_.size() != expectedSize_)
			throw improperSizeEx;
	}
};

// implementation for arrays of fundamental types xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template<typename T, szt W>
class Par<std::array<T,W>, false, std::enable_if_t<std::is_fundamental<T>::value>> : public ParBase<T>
{
	typedef std::array<T,W> Q;

	using ParBase<T>::name;
	using ParBase<T>::isLoaded_;

	Q p_;
	
public:	
	Par( const std::string& name )
		: ParBase<T>::ParBase(name)
	{}
	Par( const std::string& name, 
		 const std::string& fname, 
		 Oel& oel,
		 const std::vector<Q>& range )
		: ParBase<T>::ParBase(name) 
	{
		this->load(fname);
		check_range(range, oel);
		print(oel);
	}

	void check_range( const std::vector<Q>& r, Oel& oel )
	{
//		if (isDiscrete) {
//			if (std::find(r.begin(), r.end(), p_) == r.end())
//				throw ParOutOfRangeException<Q,isDiscrete>(name(), p_, r, oel);
//		}
//		else {
	for (szt i=0; i<W; i++)
			if (p_[i]<r[0][i] || p_[i]>r[1][i])
				throw ParOutOfRangeException<T,false>(name(), p_[i], {r[0][i], r[1][i]}, oel);
//		}
	}

	void print( Oel& oel ) final
	{
		oel.print(name(), p_);
	}
	
	Q operator()() const 
	{ 
		XASSERT(true, name()); 
		return p_; 
	}
	
	T operator[]( const szt& i ) const 
	{
		XASSERT(isLoaded_, name()); 
		XASSERT(i<W, name()); 
		return p_[i];
	}

private:
	
	void initialize( std::string value ) final
	{
		const std::string emp {" "};
		const std::string tab {"\t"};
		szt i(0);
		while (value.length()) {
			if (i == W)
				throw MyException("Improper Config::" + name() + " initialization: Excessive data size");
			ulong e {value.find(emp)};
			if (e == std::string::npos) e = value.find(tab);
			if (e == std::string::npos) e = value.length();
			const std::string val {value.substr(0, e)};
			if (val.length() < 1)
				throw MyException("Error in config file: Number of elelments in " + name() +
								  " is " + STR(p_.size()) + " which is insufficient");
			std::stringstream(val) >> p_[i];
			value.erase(0, e);
			while (!value.substr(0, 1).compare(emp) ||
				   !value.substr(0, 1).compare(tab))
				value.erase(value.begin());
			i++;
		}
		if (i < W)
			throw MyException("Improper Config::" + name() + " initialization: Data size insufficient");
	}
};

}
#endif /* Par_h */
