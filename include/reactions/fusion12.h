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
* @file fusion12.h
* @brief Contains class encapsulating tip-to-side fusion reaction.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_FUSION12_H
#define MITOSIM_FUSION12_H

#include "utils/common/misc.h"
#include "utils/msgr.h"
#include "utils/stochastic/gillespie.h"
#include "utils/stochastic/reaction.h"

#include "fusion.h"

namespace mitosim {
    
/**
 * Reaction slot for fusion of a degree 1 node with a degree 2 node.
 * @tparam Ntw The network class.
 */
template<typename Ntw>
class Fusion12
    : public Fusion<1,2,Ntw> {
    
    using Reaction = utils::stochastic::Reaction<RandFactory>;

    friend utils::stochastic::Gillespie<Reaction,RandFactory>;

public:

    using Msgr = utils::Msgr;

    /**
     * Constructor.
     * @param msgr Output message processor.
     * @param ind Reaction id.
     * @param netw The network object.
     * @param rate Rate constant.
     */
    Fusion12(
            Msgr& msgr,
            const szt ind,
            Ntw& netw,
            const real rate
        )
        : Fusion<1,2,Ntw>{msgr, ind, netw, rate, name}
    {}

    /// Sets the Gillespie score for this reaction.
    void set_score() noexcept override;

    /**
     * Gillespie score for this reaction.
     * @result Total weight of this reaction in the simulation set.
     */
    real get_score() const noexcept override { return *score; };

    /**
     * Updates propensity for a pair of network components.
     * @param c1 Index of the 1st component to update.
     * @param c2 Index of the 2nd component to update.
     */
    void update_prop(szt c1, szt c2) noexcept override;

    /// Execute the raction event.
    void fire() noexcept override;

     /// The number of times this reaction was fired.
    szt event_count() const noexcept override
    {
        return eventCount;
    }

    /**
     * Prints the reaction parameters.
     * @param le True if new line after the output.
     */
    void print(bool le) const override;

private:

    using Fusion<1,2,Ntw>::cc;
    using Fusion<1,2,Ntw>::netw;
    using Fusion<1,2,Ntw>::print;
    using Fusion<1,2,Ntw>::update_netw_stats;
    using Reaction::msgr;
    using Reaction::rate;

    static const std::string name;  ///< Reaction name constant.

    real* score {};       ///< Current rate as seen by the Gillespie reactor.
    szt   eventCount {};  ///< Number of times this reaction was fired.

    /// Total propensity for this reaction over all network components.
    szt propTotal {};

    /// Sets this reaction propensity for the whole network.
    void set_prop() noexcept override;

    /**
     * Attaches this score to the Gillespie mechanism.
     * @param a Placeholder in the Gillespie object responsible for this
     * reaction score.
     */
    void attach_score_pointer(real* a) noexcept override
    {
        score = a;
    };
};
    
// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> const std::string Fusion12<Ntw>::name {"fu12"};


template<typename Ntw>
void Fusion12<Ntw>::
set_score() noexcept
{
    *score = rate * propTotal;
}


template<typename Ntw>
void Fusion12<Ntw>::
set_prop() noexcept
{
    propTotal = netw.fu12.set_prop();
}


template<typename Ntw>
void Fusion12<Ntw>::
update_prop( szt  /*unused*/, szt  /*unused*/ ) noexcept
{
    set_prop();
}


template<typename Ntw>
void Fusion12<Ntw>::
fire() noexcept
{
    if constexpr (verbose) print(true);

    eventCount++;

    cc = netw.fu12.fire();

    update_netw_stats();
}


template<typename Ntw>
void Fusion12<Ntw>::
print( const bool le ) const
{
    Reaction::print(false);
    msgr.template print<false>(" score ", *score);
    msgr.template print<false>(" eventCount ", eventCount);
    if (le) msgr.print("\n");
}


}  // namespace mitosim

#endif  // MITOSIM_FUSION12_H
