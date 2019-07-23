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

#include <string>
#include <vector>
#include <ctime>

#define FP32			/* comment this to switch to a double precision */
#define _DEBUG			/* toggles XASSERTs */
//#define PRINT_EDGES	/* comment this to avoid printing detailed edge information */

#include "utils/Misc.h"
#include "utils/Oel.h"
#include "utils/RandBoost.h"
using namespace Utils;
typedef RandBoost<real> RandFactory;

int time0 {};
#include "Config.h"
#include "Segment.h"
#include "Network.h"

int main( int argc, char* argv[] )
{ 
	using namespace MitoD;

	if (argc < 5)
		return simple_error("Error: not sufficient configuration data provided");

	const auto workingDir {std::string(argv[1])};	// working directory
	const auto configSuffix {std::string(argv[2])};	// application-specific suffix of the configuration file
	const auto runIni {static_cast<szt>(std::stoi(argv[3]))};	// index of the starting run
	const auto runEnd {static_cast<szt>(std::stoi(argv[4]))};	// index of the last run

	bool verbose {(argc >= 6) ? static_cast<bool>(std::stoi(argv[5]))
							  : false};

	const auto workingDirIn {workingDir};	// directory for the input
	const auto workingDirOut {workingDir};	// directory for the output

	for (szt ii=runIni; ii<=runEnd; ii++) {
		StopWatch stopwatch;
		stopwatch.start();
		std::string timestr {"Run "+STR(ii)+"started: "+stopwatch.start.str};

		const auto fname {workingDirOut+"log_m_"+STR(ii)+".txt"};
		std::ofstream logfile(fname);
		Oel oel {1, 1, logfile};
		oel.print(timestr);
		oel.print("workingDirOut = "+workingDirOut);
		oel.print("runIni = %d ", runIni);
		oel.print("runEnd = %d ", runEnd);

		Config cfg {workingDirOut, configSuffix, STR(ii), verbose, oel};

		const auto seedFileName {workingDirIn+"seeds"};
		if (!file_exists(seedFileName))
			RandFactory::make_seed(seedFileName);

		auto rnd {std::make_unique<RandFactory>(seedFileName, ii, oel)};
		const auto network {std::make_unique<Network<Segment<3>>>(cfg, ii, *rnd, oel)};

		stopwatch.stop();
		timestr = "Run "+STR(ii)+" finished: "+stopwatch.stop.str+" after "+stopwatch.duration()+" sec";
		oel.print(timestr);
		
	}

	return EXIT_SUCCESS;
}

