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
 * @file network.h
 * @brief High-level network components and functionality.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_NETWORK_H
#define MITOSIM_NETWORK_H

#include "ability_for_fusion.h"
#include "config.h"
#include "definitions.h"
#include "ntw_fission.h"
#include "ntw_fusion11.h"
#include "ntw_fusion12.h"
#include "ntw_fusion1u.h"
#include "simulation.h"

namespace mitosim {

/**
 * @brief The Network class template.
 * @details Represents a fully dynamic network, capable for
 * both fusion and division.
 * @tparam SegmentT type of the segment used by the network.
 */
template<typename SegmentT>
class Network
    :  public AbilityForFusion<SegmentT> {

public:

    using thisT = Network<SegmentT>;
    using ST = SegmentT;

    using Structure<SegmentT>::add_disconnected_segment;
    using Structure<SegmentT>::clnum;
    using Structure<SegmentT>::glm;
    using Structure<SegmentT>::msgr;
    using Structure<SegmentT>::mt;
    using Structure<SegmentT>::mtmass;
    using Structure<SegmentT>::mtnum;
    using Structure<SegmentT>::nn;

    friend Fission<thisT>;
    friend Fusion<1,0,thisT>;
    friend Fusion<1,1,thisT>;
    friend Fusion<1,2,thisT>;
    friend Fusion11<thisT>;
    friend Fusion12<thisT>;
    friend Fusion1U<thisT>;
    friend NtwFission<thisT>;
    friend NtwFusion11<thisT>;
    friend NtwFusion12<thisT>;
    friend NtwFusion1U<thisT>;
    friend Simulation<thisT>;

    RandFactory&  rnd;   ///< Random number factory.
    double        time;  ///< Current time.
    unsigned long it;    ///< Iteration counter.
    const Config<real>& cfg;   ///< Configuration.

    // Reaction slots:
    NtwFission<thisT>  fis;   ///< Slot for fission reaction.
    NtwFusion11<thisT> fu11;  ///< Slot for fusion raction of nodes degr. 1+1.
    NtwFusion12<thisT> fu12;  ///< Slot for fusion reaction of nodes degr. 1+2.
    NtwFusion1U<thisT> fu1L;  ///< Slot for fusion raction of nodes degr. 1 and a cycle.

    /**
     * @brief Constructor.
     * @param cfg Configuration object.
     * @param rnd Random number factory.
     * @param msgr Output message processor.
     */
    explicit Network(
            const Config<real>& cfg,
            RandFactory& rnd,
            Msgr& msgr
    );

    /// Produce everything necessary for the simulation to start.
    auto assemble() -> thisT*;

    /// Run simulation on this network.
    void simulate();

private:

    /// Generate the network components
    void generate_components();

    /// Update the network state variables
    void update_books() noexcept;

    /**
     * @brief Write network to a file.
     * @param startnew Start a new file vs. adding new data records.
     * @param last Is the final writeout.
     * @param itr Current simulation iteration.
     * @param t Current simulation time.
     */
    void save_mitos(
        bool startnew,
        bool last,
        szt itr,
        real t
    ) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename SegmentT>
Network<SegmentT>::
Network(
        const Config<real>& cfg,
        RandFactory& rnd,
        Msgr& msgr
    )
    : AbilityForFusion<SegmentT> {msgr}
    , rnd {rnd}
    , time {zero<double>}
    , it {}
    , cfg {cfg}
    , fis {*this}
    , fu11 {*this}
    , fu12 {*this}
    , fu1L {*this}
{}


template<typename SegmentT>
auto Network<SegmentT>::
assemble() -> thisT*
{
    generate_components();
    update_books();
    return this;
}


template<typename SegmentT>
void Network<SegmentT>::
simulate()
{
    Simulation sim {*this, rnd, time, it, msgr};
    sim.initialize()();
}


template<typename SegmentT>
void Network<SegmentT>::
generate_components()
{
    // Desired initial number of sedments.
    const szt num = cfg.mtmassini / cfg.segmassini;
    if (num < 1)
        msgr.exit("The system should have at least one segment initially");

    szt m {mtnum};      // initial number of segments
    while (mtnum - m <= num-1)
        add_disconnected_segment(cfg.segmassini);

    msgr.print("Generated mtnum ", num, " of mtmass: ", mtmass);
}


template<typename SegmentT>
void Network<SegmentT>::
update_books() noexcept
{
    this->update_structure();
}


template<typename SegmentT>
void Network<SegmentT>::
save_mitos(
    const bool startnew,
    const bool last,
    const szt itr,
    const real t
) const
{
    const auto file {last
        ? cfg.workingDirOut / (std::string("mitos_last_")+cfg.runName)
        : cfg.workingDirOut / (std::string("mitos_")     +cfg.runName)};
    const auto flags = startnew ? std::ios::binary | std::ios::trunc
                                : std::ios::binary | std::ios::app;
    std::ofstream ofs {file, flags};
    if (ofs.fail())
        msgr.print("Cannot open file: ", file);

    ofs.write(reinterpret_cast<const char*>(&t), sizeof(t));
    ofs.write(reinterpret_cast<const char*>(&mtnum), sizeof(szt));

    static szt mtnummax;
    static szt nn1max;
    static szt nn2max;
    if (!last) {
        if (startnew) {
            mtnummax = 0;
            nn1max = 0;
            nn2max = 0;
        }
        if (mtnum > mtnummax)
            mtnummax = mtnum;
    }
    for (szt q=1; q<=mtnum; q++) {
        mt[q].write(ofs);
        if (!last) {
            if (mt[q].nn[1] > nn1max) nn1max = mt[q].nn[1];
            if (mt[q].nn[2] > nn2max) nn2max = mt[q].nn[2];
        }
    }
    ofs.write(reinterpret_cast<const char*>(&mtnummax), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&nn1max), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&nn2max), sizeof(szt));

    szt nst2save = last
                 ? szt{}
                 : static_cast<szt>(itr / cfg.saveFrequency);
    ofs.write(reinterpret_cast<const char*>(&nst2save), sizeof(szt));
}

}  // namespace mitosim

#endif  // MITOSIM_NETWORK_H
