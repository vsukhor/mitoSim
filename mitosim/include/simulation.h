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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"
#include "utils/common/gillespie.h"

#include "reactions/fission.h"
#include "reactions/fusion11.h"
#include "reactions/fusion12.h"
#include "reactions/fusion1u.h"

namespace MitoD {

/**
 * The Simulation class.
 * Handles the overall simulation process and its termination.
 * The reactions are encapsulated inside the Gillespie object, constructed here.
 * Controlls the output.
 */
template<typename Ntw>
class Simulation {

public:

	/** Constructor
	 * @param netw the network to be simulated
	 * @param rnd random number factory
	 * @param time current time
	 * @param it iteration counter
	 * @param msgr logging facility
	 */
	Simulation(
		Ntw& netw,
		RandFactory& rnd,
		double& time,
		ulong& it,
		Msgr& msgr
	);

	void operator()();	/**< Runs the simulation */

private:

	Ntw& netw;	/**< ref: simulated network */

	// Convenience references to some data fields of the network
	Msgr&			msgr;			/**< ref: logging facility */
	RandFactory&	rnd;			/**< ref: random number factory */
	double&			time;			/**< ref: current time */
	ulong&			it;				/**< ref: iteration counter */

	// Output parameters
	szt		logFrequency;	/**< frequency of short output to a log line */
	szt		saveFrequency;	/**< frequency of detailed output to flie */

	Gillespie<RandFactory,Reaction> gsp;			/**< Gillespie reactor controlling the simulation */

	void populateRc();					/**< Adds reactions to the Gillespie simulator */
	void terminate(const std::string&);	/**< Terminates upon the reactant exhaustion */

	// Logging
	void update_log();					/**< Outputs status summary to a log file */
	void update_log(std::ostream&);		/**< Outputs status summary to a log file */
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
Simulation<Ntw>::
Simulation(
		Ntw& netw,
		RandFactory& rnd,
		double& time,
		ulong& it,
		Msgr& msgr
	)
	: netw {netw}
	, msgr {msgr}
	, rnd {rnd}
	, time {time}
	, it {it}
	, logFrequency {netw.cfg.logFrequency}
	, saveFrequency {netw.cfg.saveFrequency}
	, gsp {rnd}
{
	populateRc();
	gsp.initialize();
}

template<typename Ntw>
void Simulation<Ntw>::
populateRc()
{
	szt ind {};
	if (netw.cfg.use_fission  )	gsp.add_reaction(std::make_unique<Fission <Ntw>>(msgr, ind++, netw, netw.cfg.rate_fission, it, time));
	if (netw.cfg.use_11_fusion)	gsp.add_reaction(std::make_unique<Fusion11<Ntw>>(msgr, ind++, netw, netw.cfg.fusion_rate_11, it, time));
	if (netw.cfg.use_12_fusion)	gsp.add_reaction(std::make_unique<Fusion12<Ntw>>(msgr, ind++, netw, netw.cfg.fusion_rate_12, it, time));
	if (netw.cfg.use_1L_fusion)	gsp.add_reaction(std::make_unique<Fusion1U<Ntw>>(msgr, ind++, netw, netw.cfg.fusion_rate_1L, it, time));
}

template<typename Ntw>
void Simulation<Ntw>::
operator()()
{
	netw.update_nn();
	netw.update_books();
	netw.save_mitos(true, false, 0, zero<real>);
	if (it % logFrequency == 0)
		update_log();

	// main loop
	while (time < netw.cfg.timeTotal) {
		it++;
		if (!gsp.set_asum()) {
			terminate("\nNo reaction left! Termination due to reaction *score == 0 for all reactions used.");
			break; 
		}
		XASSERT(!std::isnan(gsp.tau()), "Tau is nan\n");
		
		gsp.fire(time);

		if (it % saveFrequency == 0) {
			netw.save_mitos(false, false, it, time);		// appended
		}
		if (it % logFrequency == 0)
			update_log();
		if (!netw.mtnum) {
			terminate("No segments left! Termination due to chondriome exhaustion.");
			break;
		};
	}

	msgr.print("\nFinal state:");
	update_log();
	netw.save_mitos(true, true, it, time);				// only the last snapshot
	msgr.print("Final mtnum: %d\n", netw.mtnum);
}

template<typename Ntw>
void Simulation<Ntw>::
terminate( const std::string& s )
{
	netw.update_nn();
	update_log();
	msgr.print(s);
}

template<typename Ntw>
void Simulation<Ntw>::
update_log()
{
	update_log(*msgr.so);
	update_log(*msgr.sl);
}

template<typename Ntw>
void Simulation<Ntw>::
update_log( std::ostream &ofs )
{
	ofs << it << " t " << time;
	gsp.log_data(ofs);
	netw.print(ofs);
	gsp.printScores(ofs);
	ofs << std::endl;
}

}	// namespace MitoD

#endif // SIMULATION_H
