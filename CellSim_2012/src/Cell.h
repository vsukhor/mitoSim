#ifndef Cell_h
#define Cell_h

#include <vector>
#include <memory>

#include "ConfigMain.h"
#include "Segment.h"
#include "Network.h"
#include "Simulation.h"

namespace MitoD {
using namespace Utils;

class Cell {

	typedef Network<Segment<contentType>> Ntw;

public:

	Cell(
			const ConfigMain& cfgMain,
			const std::string& seedfname,
			szt runIndex,
			Oel& oel
		)
		: cfgFname_{cfgMain.confname("Cell")}
		, rnd{std::make_unique<RandFactory>(seedfname, runIndex, oel)}
		, runName {STR(runIndex)}
		, chndr {cfgMain, runName, *rnd, time, it, oel}
		, sim {chndr, runName, *rnd, time, it, oel}

	{
		sim();
	}

private:

	const std::string				cfgFname_;
	std::unique_ptr<RandFactory>	rnd;
	const std::string				runName;
	Ntw								chndr;
	real							time {zero<real>};
	szt								contentTypeit {0};
	szt								it {0};
	Simulation<Ntw>					sim;
};

}
#endif /* Cell_h */
