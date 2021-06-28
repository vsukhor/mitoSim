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

/// @file fusion1u.h
/// @brief Contains class encapsulating tip-to-cycle fusion reaction.
/// @author Valerii Sukhorukov

#ifndef MITOSIM_FUSION1U_H
#define MITOSIM_FUSION1U_H

#include "utils/stochastic/gillespie.h"
#include "utils/stochastic/reaction.h"

#include "definitions.h"
#include "fusion.h"

namespace mitosim {

/// Reaction slot for fusion of a degree 1 node with a cycle segment.
/// @tparam Ntw The network class.
template<typename Ntw>
class Fusion1U
    : public Fusion<1,0,Ntw> {

    using Reaction = utils::stochastic::Reaction<RandFactory>;

    friend utils::stochastic::Gillespie<Reaction,RandFactory>;

public:

    /// Constructor.
    /// @param msgr Output message processor.
    /// @param ind Reaction id.
    /// @param netw The network object.
    /// @param rate Rate constant.
    explicit Fusion1U(
            Msgr& msgr,
            const szt ind,
            Ntw& netw,
            const real rate
        )
        : Fusion<1,0,Ntw> {msgr, ind, netw, rate, name}
    {}

    /// Sets the Gillespie score for this reaction.
    void set_score() noexcept override;

    /// Updates propensity for a pair of network components.
    /// @param c1 Index of the 1st component to update.
    /// @param c2 Index of the 2nd component to update.
    void update_prop(szt c1, szt c2) noexcept override;

    /// Executes the raction event.
    void fire() noexcept override;

    /// Prints the reaction parameters.
    /// @param le True if new line after the output.
    void print(bool le) const override;

private:

    using Fusion<1,0,Ntw>::cc;
    using Fusion<1,0,Ntw>::netw;
    using Fusion<1,0,Ntw>::print;
    using Fusion<1,0,Ntw>::update_netw_stats;
    using Reaction::eventCount;
    using Reaction::msgr;
    using Reaction::rate;
    using Reaction::score;

    static const std::string name;  ///< Reaction name constant.
    
    /// Total propensity for this reaction over all network components.
    szt propTotal {};

    /// Sets this reaction propensity for the whole network.
    void set_prop() noexcept override;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> const std::string Fusion1U<Ntw>::name {"fu1L"};


template<typename Ntw>
void Fusion1U<Ntw>::
set_score() noexcept
{
    *score = rate * propTotal;
}


template<typename Ntw>
void Fusion1U<Ntw>::
set_prop() noexcept
{
    propTotal = netw.fu1L.set_prop();
}


template<typename Ntw>
void Fusion1U<Ntw>::
update_prop(const szt  /*unused*/,
            const szt  /*unused*/ ) noexcept
{
    set_prop();
}


template<typename Ntw>
void Fusion1U<Ntw>::
fire() noexcept
{
    if constexpr (verbose) print(true);

    eventCount++;
    
    cc = netw.fu1L.fire();

    update_netw_stats();
}


template<typename Ntw>
void Fusion1U<Ntw>::
print( const bool le ) const
{
    Fusion<1,0,Ntw>::print(false);
    if (le) msgr.print("\n");
}


}  // namespace mitosim

#endif  // MITOSIM_FUSION1U_H
