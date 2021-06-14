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

/**
* @file simulation.h
* @brief The simulation control class.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_SIMULATION_H
#define MITOSIM_SIMULATION_H

#include "utils/common/misc.h"
#include "utils/msgr.h"
#include "utils/stochastic/gillespie.h"

#include "reactions/fission.h"
#include "reactions/fusion11.h"
#include "reactions/fusion12.h"
#include "reactions/fusion1u.h"

namespace mitosim {

/**
 * @brief The Simulation class.
 * @details Handles the overall simulation process and its termination.
 * The reactions are encapsulated inside the Gillespie object, constructed here.
 * Controlls the output.
 */
template<typename Ntw>
class Simulation {

public:

    /**
     * @brief Constructor
     * @param netw the network to be simulated
     * @param rnd random number factory
     * @param time current time
     * @param it iteration counter
     * @param msgr Output message processor.
     */
    explicit Simulation(
        Ntw& netw,
        RandFactory& rnd,
        double& time,
        ulong& it,
        utils::Msgr& msgr
    );

    /// Make everything ready for start.
    auto initialize() -> Simulation<Ntw>&;

    void operator()();  ///< Runs the simulation.

private:

    Ntw& netw;  ///< ref: Simulated network.

    // Convenience references to some data fields of the network.
    utils::Msgr& msgr;   ///< ref: Output message processor.
    RandFactory& rnd;    ///< ref: random number factory.
    double&      time;   ///< ref: current time.
    ulong&       it;     ///< ref: iteration counter.

    // Output parameters
    szt logFrequency;   ///< Frequency of short output to a log line.
    szt saveFrequency;  ///< Frequency of detailed output to flie.

    /// Gillespie reactor controlling the simulation.
    utils::stochastic::Gillespie<RandFactory,
                                 utils::stochastic::Reaction<RandFactory>> gsp;

    void populateRc();  ///< Add reactions to the Gillespie simulator.

    /**
     * @brief Terminate the simulation early due to the reactant exhaustion.
     * @param s Message to pring on termination.
     */
    void terminate(const std::string& s);

    // Logging
    void update_log();               ///< Output status summary to a log file.
    void update_log(std::ostream&);  ///< Output status summary to a log file.
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
Simulation<Ntw>::
Simulation(
        Ntw& netw,
        RandFactory& rnd,
        double& time,
        ulong& it,
        utils::Msgr& msgr
    )
    : netw {netw}
    , msgr {msgr}
    , rnd {rnd}
    , time {time}
    , it {it}
    , logFrequency {netw.cfg.logFrequency}
    , saveFrequency {netw.cfg.saveFrequency}
    , gsp {rnd}
{}

template<typename Ntw>
auto Simulation<Ntw>::
initialize() -> Simulation<Ntw>&
{
    populateRc();
    gsp.initialize();

    return *this;
}


template<typename Ntw>
void Simulation<Ntw>::
populateRc()
{
    szt ind {};
    if (netw.cfg.use_fission)
        gsp.add_reaction(std::make_unique<Fission <Ntw>>(
            msgr, ind++, netw, netw.cfg.rate_fission));

    if (netw.cfg.use_11_fusion)
        gsp.add_reaction(std::make_unique<Fusion11<Ntw>>(
            msgr, ind++, netw, netw.cfg.fusion_rate_11));

    if (netw.cfg.use_12_fusion)
        gsp.add_reaction(std::make_unique<Fusion12<Ntw>>(
            msgr, ind++, netw, netw.cfg.fusion_rate_12));

    if (netw.cfg.use_1L_fusion)
        gsp.add_reaction(std::make_unique<Fusion1U<Ntw>>(
            msgr, ind++, netw, netw.cfg.fusion_rate_1L));
}


template<typename Ntw>
void Simulation<Ntw>::
operator()()
{
    netw.update_node_numbers();
    netw.update_books();
    netw.save_mitos(true, false, 0, utils::common::zero<real>);
    if (it % logFrequency == 0)
        update_log();

    // main loop
    while (time < netw.cfg.timeTotal) {
        it++;
        if (!gsp.set_asum()) {
            terminate(std::string("\nNo reaction left! ") +
                      "Termination due to reaction *score == 0 "+
                      "for all reactions used.");
            break; 
        }
        XASSERT(!std::isnan(gsp.tau()), "Tau is nan\n");
        
        gsp.fire(time);

        if (it % saveFrequency == 0)
            netw.save_mitos(false, false, it, time);  // appended

        if (it % logFrequency == 0)
            update_log();

        if (!netw.mtnum) {
            terminate("No segments left! Termination due to chondriome exhaustion.");
            break;
        };
    }

    msgr.print("\nFinal state:");
    update_log();
    netw.save_mitos(true, true, it, time);   // only the last snapshot
    msgr.print("Final mtnum: ", netw.mtnum, "\n");
}


template<typename Ntw>
void Simulation<Ntw>::
terminate( const std::string& s )
{
    netw.update_node_numbers();
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

}  // namespace mitosim

#endif  // MITOSIM_SIMULATION_H
