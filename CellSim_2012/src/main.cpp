
#define __LINUX__

#include <string>
#include <vector>
#include <ctime>

#define FP32			// comment this out to switch to double precision
#define _DEBUG
//#define PRINT_EDGES

#include "utils/Misc.h"
#include "utils/Oel.h"
#include "utils/RandBoost.h"
using namespace Utils;
typedef RandBoost<real> RandFactory;

int time0 {0};
const int contentType {2};

#include "ConfigMain.h"
#include "Cell.h"

int main( int argc, char* argv[] )
{ 
	using namespace MitoD;

	if (argc < 3)
		return simple_error("Error: no config file provided");

	const auto workingDir = std::string(argv[1]);
	const auto configSuffix = std::string(argv[2]);

	ConfigMain cfg {workingDir, configSuffix};
	
	for (szt ii=cfg.runIni; ii<=cfg.runEnd; ii++) {

		auto timestr {std::string("Run started: ")+get_current_time()};

		const auto runName {STR(ii)};
		const auto fname {cfg.workingDirOut+"log_m_"+runName+".txt"};
		std::ofstream logfile(fname);
		Oel oel {1, 1, logfile};
		
		oel.print(timestr);
		cfg.print(oel);

		const auto seedFileName {cfg.workingDirIn+"seeds"};
		if (!file_exists(seedFileName))
			RandFactory::make_seed(seedFileName);

		const auto cell {std::make_unique<Cell>(cfg, seedFileName, ii, oel)};

		timestr = "Run finished: "+get_current_time();
		oel.print(timestr);
		
	}

	return EXIT_SUCCESS;
}

