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

#include "utils/stop_watch.h"

#include "definitions.h"
#include "config.h"
#include "network.h"
#include "segment.h"

int main( int argc, char* argv[] )
{
    constexpr int minArgc = 5;
    if (argc < minArgc) {
        std::cerr << "Error: not sufficient configuration data provided\n";
        std::exit(EXIT_FAILURE);
    }

    // Working directory:
    const std::filesystem::path workingDir {std::string(argv[1])};
    const auto workingDirIn = workingDir;    // directory for the input
    const auto workingDirOut = workingDir;   // directory for the output

    // Application-specific suffix of the configuration file:
    const auto configSuffix = std::string(argv[2]);

    const auto runIni = std::stoi(argv[3]);  // index of the starting run
    const auto runEnd = std::stoi(argv[4]);  // index of the last run

    // Loop over the simulation runs:
    for (int runInd=runIni; runInd<=runEnd; runInd++) {

        utils::StopWatch stopwatch;
        stopwatch.start();

        // Set the logging:
        const auto logf {workingDirOut /
                         (std::string("log_m_")+std::to_string(runInd)+".txt")};
        std::ofstream logfile {logf};
        constexpr int printPrecision = 6;
        mitosim::Msgr msgr {&std::cout, &logfile, printPrecision};

        // Report the environment:
        msgr.print("Run ", runInd, " started: ", stopwatch.start.str);
        msgr.print("workingDirOut = ", workingDirOut);
        msgr.print("runIni = ", runIni);
        msgr.print("runEnd = ", runEnd);

        // Import the configuration settings:
        mitosim::Config<mitosim::real> cfg {
            workingDirOut, configSuffix, std::to_string(runInd), msgr
        };

        auto rnd = std::make_unique<mitosim::RandFactory>(runInd, msgr);

        // Create and simulate the network:
        constexpr int maxNodeDegree = 3;
        const auto network =
            std::make_unique<
                mitosim::Network<mitosim::Segment<maxNodeDegree>>
                    >(cfg, *rnd, msgr);
        network->assemble()->simulate();

        // Finalize:
        stopwatch.stop();
        msgr.print("Run ", runInd, " finished: "+stopwatch.stop.str,
                   " after ", stopwatch.duration(), " sec\n");
    }

    return EXIT_SUCCESS;
}

