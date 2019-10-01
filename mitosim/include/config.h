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

#ifndef CONFIG_H
#define CONFIG_H

#include <string>

#include "utils/common/msgr.h"
#include "config/reader.h"

namespace MitoD {

/**
 * Encapsulates and reads configuration parameters.
 * Reads from a file and stores initial configuration parameters for the hole modeling session.
 */
class Config {

	const std::string	configSuffix;	///< Application-specific suffix of the configuration file.
	const std::string	cfgFname;		///< Configuration file name.
	::Config::Reader	read;			///< Generic reader of configuraion files.

public:
/
	const std::string 	runName;		///< Current run index as a string.
	const std::string	workingDirOut;	///< Output directory.

	const real	timeTotal;				///< Total time to simulate.
	const szt	logFrequency;			///< (iterations): Every logfreq iteration steps log is updated.
	const szt	saveFrequency;			///< (iterations): Every saveFrequency steps save_mitos is executed.

	const real	edgeLength;				///< (micrometers) Edge Length.
	const szt	mtmassini;				///< Total chondriome length, edges.
	const szt	segmassini;				///< Segment lengths at the beginning of simulation.

	// FISSION
	const bool	use_fission;			///< A flag for the activation of the fission reaction.
	const real	rate_fission;	 		///< Basic probability of breaking up a junction.

	// FUSION
	const bool	use_11_fusion;			///< A flag for the activation of the fussion_11 reaction.
	const real	fusion_rate_11;			///< Probability of a free end to bind to another free end.
	const bool	use_12_fusion;			///< A flag for the activation of the fusion_12 reaction.
	const real	fusion_rate_12;			///< Probability of a free end to bind to a side.
	const bool	use_1L_fusion;			///< A flag for the activation of the fusion_1L reaction.
	const real	fusion_rate_1L;			///< Probability of a free end to bind to a separate cycle junction node.

	/**
	 * \brief Constructor.
	 * @par workingDirOut directory for the output
	 * @par configSuffix application-specific suffix of the configuration file
	 * @par runName run index
	 * @par msgr logging facility
	 */
	explicit Config(
			const std::string& workingDirOut,
			const std::string& configSuffix,
			const std::string& runName,
			Msgr& msgr
			)
		: configSuffix {configSuffix}
		, cfgFname {workingDirOut+"config"+"_"+configSuffix+".txt"}
		, read {cfgFname, msgr}
		, runName {runName}
		, workingDirOut {workingDirOut}
		, timeTotal			{read("timeTotal", zerohuge<real>)}
		, logFrequency		{read("logFrequency", onehuge<szt>)}
		, saveFrequency		{read("saveFrequency", onehuge<szt>)}

		, edgeLength		{read("edgeLength", zerohuge<real>)}
		, mtmassini			{read("mtmassini", onehuge<szt>)}
		, segmassini		{read("segmassini", onehuge<szt>)}

	// FISSION
		, use_fission		{read("use_fission", bools)}
		, rate_fission		{read("rate_fission", zerohuge<real>)}

	// FUSION
		, use_11_fusion		{read("use_11_fusion", bools)}
		, fusion_rate_11	{read("fusion_rate_11", zerohuge<real>)}
		, use_12_fusion		{read("use_12_fusion", bools)}
		, fusion_rate_12	{read("fusion_rate_12", zerohuge<real>)}
		, use_1L_fusion		{read("use_1L_fusion", bools)}
		, fusion_rate_1L	{read("fusion_rate_1L", zerohuge<real>)}
	{}
};

}	// namespace MitoD

#endif // CONFIG_H
