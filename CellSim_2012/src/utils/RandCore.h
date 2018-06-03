
#ifndef RandCore_h
#define RandCore_h

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "Misc.h"
#include "Errors.h"

namespace Utils {
	
template<typename realT> 
class RandCore {

	// Assure that the template parameter is a floating type
	static_assert(std::is_floating_point<realT>::value,
				  "class RandCore can only be instantiated with floating point types");

public:

	static const szt num_saved_seeds {1000001};		// # of seeds in 'seed' file
	static const int bufferSize {1000000};
	static const int mainSeed {1'234'567'890};

	static void make_seed(const std::string& seedFname);

	static int readin_seed(const std::string& seedFname, 
						   szt ii, 
						   Oel &oel);	

	uint theSeed() { return seed; }

protected:

	RandCore(Oel&, uint, const std::string&) noexcept;
	RandCore(Oel&, const std::string&, szt);

private:

	uint	seed;
	Oel&	oel;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename realT> 
RandCore<realT>::
RandCore(Oel& oel,
		 uint seed,
		 const std::string& runName) noexcept
	: seed {seed}
	, oel {oel}
{
	oel.print("RUN "+runName);
	oel.print("SEED = %d", this->seed);
}

template<typename realT> 
RandCore<realT>::
RandCore(Oel& oel,
		 const std::string& seedFname,
		 szt runInd)
	: oel{oel}
{
	if (!file_exists(seedFname)) 
			make_seed(seedFname);
	seed = readin_seed(seedFname, runInd, oel);
	oel.print("RUN = %d", runInd);
	oel.print("SEED = %d", seed);
}

template<typename realT> 
void RandCore<realT>::
make_seed( const std::string& seedFname )
{
	std::cout << "No seed file found. Creating a new seed file "+seedFname << std::endl;
	std::mt19937 g;
	g.seed(mainSeed);
	
	std::uniform_int_distribution<uint> seed_d(100000000, 2100000000);
	std::ofstream file {seedFname, std::ios::binary};
	if (!file.is_open())
		exit(simple_error("Unable to create seed file "+seedFname));

	for (szt i=0; i<num_saved_seeds; i++) {
		uint s = seed_d(g);
		file.write((char*)(&s), sizeof(uint));
	}
}

template<typename realT> 
int RandCore<realT>::
readin_seed(const std::string& seedFname,
			szt runInd,
			Oel& oel)
{
	oel.print(("Reading from file "+seedFname+" seed no: %d").c_str(), runInd);
	
	std::ifstream file {seedFname, std::ios::binary};
	if (!file.is_open())
		exit(simple_error("Unable to open file "+seedFname));
	file.seekg(static_cast<std::fstream::off_type>(runInd*sizeof(uint)), file.beg);
	uint seed;
	file.read(reinterpret_cast<char*>(&seed), sizeof(uint));

	return seed;
}

}
#endif /* Rand_h */
