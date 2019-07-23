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

#ifndef Config_h
#define Config_h

#include <array>
#include <vector>
#include <string>

#include "utils/Misc.h"
#include "utils/Par.h"
#include "utils/Errors.h"

namespace MitoD {

/**
 * The ConfigReader class.
 * Assists in reading the confuguration file.
 */
class ConfigReader {

public:

	/** Constructor.
	 * @param cfgFname configuration file name
	 * @param oel logging facility object
	 */
	ConfigReader( const std::string& cfgFname,
				  Oel& oel
				)
		: cfgFname {cfgFname}
		, oel {oel}
	{}

	bool operator()( std::string s, const std::vector<bool>& range )
	{
		return Par<bool,true>(s, cfgFname, oel, range)();
	}

	real operator()( std::string s, const std::vector<real>& range )
	{
		return Par<real,false>(s, cfgFname, oel, range)();
	}
	szt operator()( std::string s, const std::vector<szt>& range )
	{
		return Par<szt,false>(s, cfgFname, oel, range)();
	}

private:

	const std::string	cfgFname;		/**< configuration file name */
	Oel&				oel;			/**< logging facility */
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * The Config class.
 * Reads from a file and stores initial configuration parameters for the hole modeling session.
 */
class Config {

	Config(const Config&) = delete;				/**< non construction-copyable: single parameter class per simulation */
	Config& operator=(const Config&) = delete;	/**< non copyable: single parameter class per simulation */

	const std::string	configSuffix;		/** application-specific suffix of the configuration file */
	const std::string	cfgFname;			/** configuration file name */
	Oel&				oel;				/** logging facility */
	ConfigReader		read;				/** an auxiliary object */

	// The vectors below are used for chacking parameter value boundaries
	const std::vector<bool> bools {false, true};		/**< allowed booleans */
	const std::vector<real>	allT {-huge<real>,huge<real>};	/**< allowed all reals */
	const std::vector<real>	nonNegT {0.,huge<real>};		/**< allowed non-negative reals */
	const std::vector<szt>	posSzt {1, huge<szt>};			/**< allowed positive integers */

public:

	const std::string 	runName;
	const std::string	workingDirOut;

	const bool		verbose;				/**< work in verbose mode */
	const real		timeTotal;				/**< total time to simulate */
	const szt		logFrequency;			/**< (iterations): every logfreq iteration steps log is updated */
	const szt		saveFrequency;			/**< (iterations): every saveFrequency steps save_mitos is executed */

	const real		edgeLength;				/**< (micrometers) */
	const szt		mtmassini;				/**< total chondriome length, edges */
	const szt		segmassini;				/**< segment lengths at the beginning of simulation */

	// FISSION
	const bool		use_fission;			/**< a flag for the activation of the fission reaction */
	const real		rate_fission;	 		/**< basic probability of breaking up a junction */

	// FUSION
	const bool		use_11_fusion;			/**< a flag for the activation of the fussion_11 reaction */
	const real		fusion_rate_11;			/**< probability of a free end to bind to another free end */
	const bool		use_12_fusion;			/**< a flag for the activation of the fusion_12 reaction */
	const real		fusion_rate_12;			/**< probability of a free end to bind to a side */
	const bool		use_1L_fusion;			/**< a flag for the activation of the fusion_1L reaction */
	const real		fusion_rate_1L;			/**< probability of a free end to bind to a separate cycle junction node */

	/**
	 * Constructor.
	 * @par workingDirOut directory for the output
	 * @par configSuffix application-specific suffix of the configuration file
	 * @par runName run index
	 * @par verbose report in verbose mode
	 * @par oel logging facility
	 */
	explicit Config(
			const std::string& workingDirOut,
			const std::string& configSuffix,
			const std::string& runName,
			const bool verbose,
			Oel& oel
			)
		: configSuffix {configSuffix}
		, cfgFname {workingDirOut+"config"+"_"+configSuffix+".txt"}
		, oel {oel}
		, read {cfgFname, oel}
		, runName {runName}
		, workingDirOut {workingDirOut}
		, verbose {verbose}
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

#endif /* Config_h */
