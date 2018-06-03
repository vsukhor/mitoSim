#ifndef Simulation_h
#define Simulation_h

#include "utils/Oel.h"
#include "utils/Gillespie.h"

#include "Reactions/Fission.h"
#include "Reactions/Fusion11.h"
#include "Reactions/Fusion12.h"
#include "Reactions/Fusion1L.h"

namespace MitoD {

template<typename Ntw>
class Simulation {

public:

	Simulation(
		Ntw& netw,
		const std::string& runName,
		RandFactory& rnd,
		real& time,
		szt& it,
		Oel& oel
	);

	void operator()();

private:

	Ntw&					netw;
	Oel&					oel;
	const std::string&		runName;
	RandFactory&			rnd;
	szt						logFrequency;
	szt						saveFrequency;
	real&					time;
	szt&					it;
	bool					verbose;
	bool					isEquilibrated {false};
	Gillespie<Reaction>		gsp;

	void populateRc();
	void terminate(const std::string&);

	void log();
	void log(std::ostream&);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
Simulation<Ntw>::
Simulation(
		Ntw& netw,
		const std::string& runName,
		RandFactory& rnd,
		real& time,
		szt& it,
		Oel& oel
	)
	: netw {netw}
	, oel {oel}
	, runName {runName}
	, rnd {rnd}
	, logFrequency {netw.logFrequency}
	, saveFrequency {netw.saveFrequency}
	, time {time}
	, it {it}
	, verbose {netw.verbose}
	, gsp {rnd}
{
	populateRc();
	gsp.initialize();
}

template<typename Ntw>
void Simulation<Ntw>::
populateRc() {
	szt ind {0};
	if (netw.use_fission  )	gsp.add_reaction(std::make_unique<Fission <Ntw>>(oel, ind++, netw, netw.rate_fission, it, time, netw.verbose));
	if (netw.use_11_fusion)	gsp.add_reaction(std::make_unique<Fusion11<Ntw>>(oel, ind++, netw, netw.fusion_rate_11, it, time, netw.verbose, 1, 1));
	if (netw.use_12_fusion)	gsp.add_reaction(std::make_unique<Fusion12<Ntw>>(oel, ind++, netw, netw.fusion_rate_12, it, time, netw.verbose, 1, 2));
	if (netw.use_1L_fusion)	gsp.add_reaction(std::make_unique<Fusion1L<Ntw>>(oel, ind++, netw, netw.fusion_rate_1L, it, time, netw.verbose, 1, 2));
}

template<typename Ntw>
void Simulation<Ntw>::
operator()()
{
	netw.update_nn();
	netw.update_books();
//	netw.print_mitos( "IT " );
	netw.io.save_mitos(true, false, 0, zero<real>);
	if (it % logFrequency == 0)
		log();
	while (time < netw.timeTotal) {
		it++;
		if (!gsp.set_asum()) {
			terminate("\nNo reaction left! Termination due to reaction *score == 0 for all reactions used.");
			break; 
		}
		XASSERT(!std::isnan(gsp.tau()), "Tau is nan");
		
		time += gsp.tau();
		gsp.make();

		if (it % saveFrequency == 0) {
//			netw.save_mitos( false, false, it, t );		// appended
//			netw.save_mitos( true, true, it, t );		// only the last snapshot
//			netw.save_edges( false, "removed_segments", false, t, netw.segr );		// appended
//			netw.save_edges( true, "removed_segments", true, t, netw.segr );		// only the last snapshot
		}
		if (it % logFrequency == 0)
			log();
		if (!netw.mtnum) {
			terminate("No mitos left! Termination due to chondriome exhaustion.");
			break;
		};
	}
	oel.print("\n Final state:.");
	log();
	netw.io.save_mitos(true, true, it, time);				// only the last snapshot
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
	log();
	oel.print(s);
}

template<typename Ntw>
void Simulation<Ntw>::
log()
{
	log(oel.so);
	log(oel.sl);
}
template<typename Ntw>
void Simulation<Ntw>::
log( std::ostream &ofs )
{
	ofs << it << " t " << time;
	gsp.log_data(ofs);
	netw.print(ofs);
	gsp.printScores(ofs);
	ofs << std::endl;
}

}

#endif /* Simulation_h */
