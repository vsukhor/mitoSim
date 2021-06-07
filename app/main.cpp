/* =============================================================================
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
================================================================================
*/

#include <filesystem>
#include <ostream>
#include <string>
#include <vector>

#include "utils/common/constants.h"
#include "utils/common/misc.h"
#include "utils/common/msgr.h"
#include "utils/common/stop_watch.h"
#include "utils/random/with_boost.h"

// #define USE_UTILS_XASSERT  // toggles XASSERTs.
// #define PRINT_EDGES  // comment this to avoid printing detailed edge info.

namespace mitosim {

using real = float;
using RandFactory = utils::random::Boost<real>;
constexpr bool verbose {};   ///< Work in verbose mode.

}   // namespace mitosim

#include "config.h"
#include "network.h"
#include "segment.h"

int main( int argc, char* argv[] )
{
    using utils::common::STR;
    using utils::common::szt;

    constexpr const int MIN_ARGC = 5;
    if (argc < MIN_ARGC)
        return utils::common::exceptions::simple(
            "Error: not sufficient configuration data provided", nullptr);

    // Working directory:
    const auto workingDir = std::filesystem::path {std::string(argv[1])};
    // Application-specific suffix of the configuration file:
    const auto configSuffix = std::string(argv[2]);
    // Index of the starting run:
    const auto runIni = static_cast<szt>(std::stoi(argv[3]));
    // Index of the last run:
    const auto runEnd = static_cast<szt>(std::stoi(argv[4]));

    const auto workingDirIn = workingDir;    // directory for the input
    const auto workingDirOut = workingDir;   // directory for the output

    for (szt ii=runIni; ii<=runEnd; ii++) {
        utils::common::StopWatch stopwatch;
        stopwatch.start();

        const auto logf = std::filesystem::directory_entry {
            workingDirOut / (STR("log_m_")+STR(ii)+".txt")
        };
        std::ofstream logfile {logf};
        constexpr const int PRINT_PRECISION {6};
        utils::common::Msgr msgr {&std::cout, &logfile, PRINT_PRECISION};
        msgr.print("Run "+STR(ii)+" started: "+stopwatch.start.str);
        msgr.print("workingDirOut = "+workingDirOut.string());
        msgr.print("runIni = "+STR(runIni));
        msgr.print("runEnd = "+STR(runEnd));

        mitosim::Config<mitosim::real> cfg {workingDirOut, configSuffix, STR(ii), msgr};

        const auto seeds =
            std::filesystem::directory_entry(cfg.workingDirOut / "seeds");
        if (!seeds.is_regular_file())
            mitosim::RandFactory::make_seed(seeds, &msgr);

        auto rnd = std::make_unique<mitosim::RandFactory>(seeds, ii, msgr);
        constexpr const int MAX_NODE_DEGREE {3};
        const auto network =
            std::make_unique<
                mitosim::Network<mitosim::Segment<MAX_NODE_DEGREE>>
                    >(cfg, *rnd, msgr);
        network->assemble()->simulate();
        stopwatch.stop();
        msgr.print("Run "+STR(ii)+
                   " finished: "+stopwatch.stop.str+
                   " after "+stopwatch.duration()+" sec\n");

    }

    return EXIT_SUCCESS;
}

