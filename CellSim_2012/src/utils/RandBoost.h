#ifndef RandBoost_h
#define RandBoost_h

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include "Misc.h"
#include "RandCore.h"

namespace Utils {

template<typename realT>
class RandBoost : public RandCore<realT> {

public:

	RandBoost() {}
	RandBoost(uint, const std::string&, Oel&);
	RandBoost(const std::string&, szt, Oel&);

	realT r01u();
	
	constexpr int uniform0(int);
	constexpr uint uniform0(uint);
	constexpr szt uniform0(szt);
	constexpr realT uniform0(realT max);
	
	template<typename intT>
	constexpr intT uniform1(intT);
	
private:

	using RandCore<realT>::bufferSize;

	boost::random::uniform_01<float>	flt01_unifromDistr;			
	boost::random::uniform_01<double>	dbl01_unifromDistr;			

	realT								rU01[bufferSize];
	int									rU01_ind;

	boost::mt19937						g;

	void prepare_uniform_real01();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename realT>
RandBoost<realT>::
RandBoost(uint seed,
		  const std::string& runName,
		  Oel& oel)
	: RandCore<realT>(oel, seed, runName)
	, rU01_ind {-1}
{
	g.seed(this->theSeed());
	prepare_uniform_real01();
}

template<typename realT>
RandBoost<realT>::
RandBoost(const std::string& seedFname,
		  szt ii,
		  Oel& oel)
	: RandCore<realT>(oel, seedFname, ii)
	, rU01_ind {-1}
{
	
	g.seed(this->theSeed());
	prepare_uniform_real01();
}

// Generates real random numbers with uniform distribution over [0,1)
template<> inline
void RandBoost<float>::
prepare_uniform_real01()
{
	for (auto& o : rU01) 
		o = flt01_unifromDistr(g);	
}

template<> inline
void RandBoost<double>::
prepare_uniform_real01()
{
	for (auto& o : rU01) 
		o = dbl01_unifromDistr(g);	
}

// returns a random number with uniform distribution over [0,1)
template<typename realT> inline
realT RandBoost<realT>::
r01u()
{
	if (++rU01_ind == bufferSize) {
		prepare_uniform_real01();
		rU01_ind = 0;
	}
	return rU01[rU01_ind];
}

// returns int in the range [ 0, max-1 ]
template<typename realT> constexpr
int RandBoost<realT>::
uniform0(int max)
{			
	XASSERT(max > 0, "RandBoost<realT>::uniform0 requires max > 0");

	auto ir {static_cast<int>(r01u()*max)};
	
	while (ir >= max) 
		ir = static_cast<int>(r01u()*max);
	
	return ir;	
}

// returns uint in the range [ 0, max-1 ]
template<typename realT> constexpr
uint RandBoost<realT>::
uniform0(uint max)
{			
	XASSERT(max > 0, "RandBoost<realT>::uniform0 requires max > 0");

	auto ir {static_cast<uint>(r01u()*max)};
	
	while (ir >= max) 
		ir = static_cast<uint>(r01u()*max);
	
	return ir;	
}

// returns szt in the range [ 0, max-1 ]
template<typename realT> constexpr
szt RandBoost<realT>::
uniform0(szt max)
{			
	XASSERT(max > 0, "RandBoost<realT>::uniform0 requires max > 0");

	auto ir {static_cast<szt>(r01u()*max)};
	
	while (ir >= max) 
		ir = static_cast<szt>(r01u()*max);
	
	return ir;	
}

// returns outT in the range [ 1, max ]
template<typename realT>
template<typename intT> constexpr
intT RandBoost<realT>::
uniform1(intT max)
{			
	// Ensure that the template parameter is a floating type
	static_assert(std::is_integral<intT>::value,
				  "class Geometric can only be instantiated with integer types");
	
	XASSERT(max > 0, "RandBoost<realT>::uniform1 requires max > 0");
	
	return uniform0(max) + 1;
}

template<typename realT> constexpr
realT RandBoost<realT>::
uniform0(realT max)
{			
	XASSERT(max > realT(0.), "RandBoost<realT>::uniform0 requires max > 0");

	auto ir {r01u() * max};
	
	while (ir >= max) 
		ir = r01u() * max;
	
	return ir;	
}

}
	
#endif /* RandBoost_h */
