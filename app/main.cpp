/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann,
   Helmholtz Center for Infection Research (Braunschweig, Germany).
   All Rights Reserved.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
============================================================================== */

#include <string>
#include <vector>

#define FP32			// comment this to switch to a double precision.
#define _DEBUG			// toggles XASSERTs.
//#define PRINT_EDGES	// comment this to avoid printing detailed edge information.

#ifdef FP32
	using real = float;
#else
	using real = double;
#endif

#include "utils/common/misc.h"
#include "utils/common/msgr.h"
#include "utils/common/stop_watch.h"
#include "utils/random/with_boost.h"

namespace MitoSim {
using namespace Utils::Common;
using RandFactory = Utils::Random::Boost<real>;
constexpr bool verbose {};		///< Work in verbose mode.
}

#include "config.h"
#include "segment.h"
#include "network.h"

int main( int argc, char* argv[] )
{
	using namespace MitoSim;

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

		MitoSim::Config cfg {workingDirOut, configSuffix, STR(ii), msgr};

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

