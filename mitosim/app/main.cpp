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

#define FP32			/* comment this to switch to a double precision */
#define _DEBUG			/* toggles XASSERTs */
//#define PRINT_EDGES	/* comment this to avoid printing detailed edge information */

#ifdef FP32
	using real = float;
#else
	using real = double;
#endif

#include "utils/common/misc.h"
#include "utils/common/msgr.h"
#include "utils/common/stop_watch.h"
#include "utils/random/with_boost.h"

namespace MitoD {
using namespace Utils::Common;
using RandFactory = Utils::Random::Boost<real>;
constexpr bool verbose {};		/**< work in verbose mode */
}

#include "config.h"
#include "segment.h"
#include "network.h"

int main( int argc, char* argv[] )
{
	using namespace MitoD;

	if (argc < 5)
		return Exceptions::simple("Error: not sufficient configuration data provided");

	const auto workingDir {std::string(argv[1])};	// working directory
	const auto configSuffix {std::string(argv[2])};	// application-specific suffix of the configuration file
	const auto runIni {static_cast<szt>(std::stoi(argv[3]))};	// index of the starting run
	const auto runEnd {static_cast<szt>(std::stoi(argv[4]))};	// index of the last run

	const auto workingDirIn {workingDir};	// directory for the input
	const auto workingDirOut {workingDir};	// directory for the output

	for (szt ii=runIni; ii<=runEnd; ii++) {
		StopWatch stopwatch;
		stopwatch.start();

		std::ofstream logfile {workingDirOut+"log_m_"+STR(ii)+".txt"};
		Msgr msgr {&std::cout, &logfile};
		msgr.print("Run "+STR(ii)+" started: "+stopwatch.start.str);
		msgr.print("workingDirOut = "+workingDirOut);
		msgr.print("runIni = %d ", runIni);
		msgr.print("runEnd = %d ", runEnd);

		MitoD::Config cfg {workingDirOut, configSuffix, STR(ii), msgr};

		const auto seedFileName {workingDirIn+"seeds"};
		if (!file_exists(seedFileName))
			RandFactory::make_seed(seedFileName, &msgr);

		auto rnd {std::make_unique<RandFactory>(seedFileName, ii, msgr)};
		const auto network {std::make_unique<Network<Segment<3>>>(cfg, ii, *rnd, msgr)};

		stopwatch.stop();
		msgr.print("Run "+STR(ii)+" finished: "+stopwatch.stop.str+"after "+stopwatch.duration()+" sec\n");
		
	}

	return EXIT_SUCCESS;
}

