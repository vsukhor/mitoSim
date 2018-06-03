#ifndef ConfigMito_h
#define ConfigMito_h

#include <array>
#include <vector>
#include <string>

#include "utils/Misc.h"
#include "utils/Par.h"
#include "utils/Errors.h"
#include "ConfigMain.h"

namespace MitoD {

class ConfigMito {

	ConfigMito( const ConfigMito& ) = delete;				// non construction-copyable: single parameter class per simulation
	ConfigMito& operator=( const ConfigMito& ) = delete;	// non copyable: single parameter class per simulation

	const std::string	cfgFname;
	Oel&				oel;
	ConfigReader		read;

	const std::vector<bool> bools {false, true};
	const std::vector<real>	allT {-huge<real>,huge<real>};
	const std::vector<real>	nonNegT {0.,huge<real>};
	const std::vector<szt>	posSzt {1, huge<szt>};

public:

	const real		timeTotal;
	const szt		logFrequency;			// (iterations): every logfreq iteration steps log is updated
	const szt		saveFrequency;			// (iterations): every savefrequency iteration steps save_mitos is executed

	const real		edgeLength;				// (micrometers)
	const szt		mtmassini;				// total chondriome length, edges
	const szt		segmassini;				// segment lengths at the beginning of simulation

	// FISSION
	const bool		use_fission;
	const real		rate_fission;	 		// basic probability of breaking up a junction

	// FUSION
	const bool		use_11_fusion;
	const real		fusion_rate_11;			// probability of a free end to bind to another free end
	const bool		use_12_fusion;
	const real		fusion_rate_12;			// probability of a free end to bind to a side
	const bool		use_1L_fusion;
	const real		fusion_rate_1L;			// probability of a free end to bind to a pure loop junction node

ConfigMito( const MitoD::ConfigMain& cfgMain,
			Oel& oel
		)
	: cfgFname {cfgMain.confname("Mito")}
	, oel {oel}
	, read {cfgFname, oel}
	, timeTotal			{read("timeTotal", nonNegT)}
	, logFrequency		{read("logFrequency", posSzt)}
	, saveFrequency		{read("saveFrequency", posSzt)}

	, edgeLength		{read("edgeLength", nonNegT)}
	, mtmassini			{read("mtmassini", posSzt)}
	, segmassini		{read("segmassini", posSzt)}
// FISSION
	, use_fission		{read("use_fission", bools)}
	, rate_fission		{read("rate_fission", nonNegT)}
// FUSION
	, use_11_fusion		{read("use_11_fusion", bools)}
	, fusion_rate_11	{read("fusion_rate_11", nonNegT)}
	, use_12_fusion		{read("use_12_fusion", bools)}
	, fusion_rate_12	{read("fusion_rate_12", nonNegT)}
	, use_1L_fusion		{read("use_1L_fusion", bools)}
	, fusion_rate_1L	{read("fusion_rate_1L", nonNegT)}
{}

};

}

#endif /* ConfigMito_h */
