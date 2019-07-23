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

#ifndef Simulation_h
#define Simulation_h

#include "utils/Oel.h"
#include "utils/Gillespie.h"

#include "Reactions/Fission.h"
#include "Reactions/Fusion11.h"
#include "Reactions/Fusion12.h"
#include "Reactions/Fusion1L.h"

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
	 * @param oel logging facility
	 */
	Simulation(
		Ntw& netw,
		RandFactory& rnd,
		real& time,
		szt& it,
		Oel& oel
	);

	void operator()();	/**< Runs the simulation */

private:

	Ntw& netw;	/**< ref: simulated network */

	// Convenience references to some data fields of the network
	Oel&			oel;			/**< ref: logging facility */
	RandFactory&	rnd;			/**< ref: random number factory */
	real&			time;			/**< ref: current time */
	szt&			it;				/**< ref: iteration counter */

	// Output parameters
	szt		logFrequency;	/**< frequency of short output to a log line */
	szt		saveFrequency;	/**< frequency of detailed output to flie */
	bool	verbose;		/**< verbosity of the short logs */

	Gillespie<Reaction>	gsp;			/**< Gillespie reactor controlling the simulation */

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
		real& time,
		szt& it,
		Oel& oel
	)
	: netw {netw}
	, oel {oel}
	, rnd {rnd}
	, time {time}
	, it {it}
	, logFrequency {netw.cfg.logFrequency}
	, saveFrequency {netw.cfg.saveFrequency}
	, verbose {netw.verbose}
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
	if (netw.cfg.use_fission  )	gsp.add_reaction(std::make_unique<Fission <Ntw>>(oel, ind++, netw, netw.cfg.rate_fission, it, time, netw.cfg.verbose));
	if (netw.cfg.use_11_fusion)	gsp.add_reaction(std::make_unique<Fusion11<Ntw>>(oel, ind++, netw, netw.cfg.fusion_rate_11, it, time, netw.cfg.verbose, 1, 1));
	if (netw.cfg.use_12_fusion)	gsp.add_reaction(std::make_unique<Fusion12<Ntw>>(oel, ind++, netw, netw.cfg.fusion_rate_12, it, time, netw.cfg.verbose, 1, 2));
	if (netw.cfg.use_1L_fusion)	gsp.add_reaction(std::make_unique<Fusion1L<Ntw>>(oel, ind++, netw, netw.cfg.fusion_rate_1L, it, time, netw.cfg.verbose, 1, 2));
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
		
		time += gsp.tau();
		gsp.make();

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

	oel.print("\n Final state:.");
	update_log();
	netw.save_mitos(true, true, it, time);				// only the last snapshot
	oel.print("Final mtnum: %d", netw.mtnum);
	const int curtime = int(clock()/CLOCKS_PER_SEC);
	oel.print("Clock: %d sec",curtime-time0);
	time0 = curtime;
}

template<typename Ntw>
void Simulation<Ntw>::
terminate( const std::string& s )
{
	netw.update_nn();
	update_log();
	oel.print(s);
}

template<typename Ntw>
void Simulation<Ntw>::
update_log()
{
	update_log(oel.so);
	update_log(oel.sl);
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

}

#endif /* Simulation_h */
