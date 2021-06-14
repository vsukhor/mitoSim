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
 * @file fission.h
 * @brief Contains class encapsulating the fission reaction.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_FISSION_H
#define MITOSIM_FISSION_H

#include <string>
#include <vector>

#include "utils/common/misc.h"
#include "utils/msgr.h"
#include "utils/stochastic/gillespie.h"
#include "utils/stochastic/reaction.h"

namespace mitosim {

using utils::common::szt;
using utils::common::uint;
using utils::common::ulong;
using utils::common::vup;

template<uint, uint, typename>
class Fusion;

/**
 * Fission reaction class.
 * @tparam Ntw The network class.
 */
template<typename Ntw>
class Fission
: public utils::stochastic::Reaction<RandFactory> {

    friend utils::stochastic::Gillespie<Reaction,RandFactory>;

public:

    /**
     * Constructor.
     * @param msgr Output message processor.
     * @param ind reaction id.
     * @param netw The network object.
     * @param rate Reaction rate constant.
     */
    Fission( utils::Msgr& msgr,
             const szt ind,
             Ntw& netw,
             const real rate
       )
        : Reaction {msgr, ind, rate, name, "fission"}
        , netw {netw}
    {}


    /// Sets the Gillespie score for this reaction.
    void set_score() noexcept override;

    /**
     * Gillespie score for this reaction.
     * @result Total weight of this reaction in the simulation set.
     */
    real get_score() const noexcept override { return *score; };


    /// Updates propensity for two network components indexed in the parameters.
    void update_prop(szt c1, szt c2) noexcept override;


    /// Executes the raction event.
    void fire() noexcept override;


    /**
     * Activity status of the reaction.
     * @return True if the reaction is used in the current simulation session.
     * @param r Pointer to the reaction.
     */
    static constexpr auto is_active(const std::unique_ptr<Reaction>& r) noexcept;


    /**
     * Populates the vector of reactions that need a score update.
     * The update is performed after *this has fired
     * and initializes the propensities and effective rate.
     * @param rc Vector of unique pointers to all reactions taking part in
     * the simulation.
     */
    void initialize_dependencies(const vup<Reaction>& rc) noexcept override;


    /// Returns the number of times this reaction was fired.
    szt event_count() const noexcept override { return eventCount; }


    /**
     * Prints the reaction parameters.
     * @param le True if new line after the output.
     */
    void print(bool le) const override;


private:

    using Reaction::msgr;
    using Reaction::rate;

    Ntw& netw;  ///< ref: The network.

    std::array<szt,2> cc;

    /// Current rate as seen by the Gillespie reactor.
    real* score {};

    /// Number of times this reaction was fired.
    szt eventCount {};

    /// Reactions that need a score update after *this has fired.
    std::vector<Reaction*> dependents;

    /// Reaction name constant.
    static const std::string name;

    /// Set this reaction propensity for the whole network.
    void set_prop() noexcept override;

    /// Attach this score to the Gillespie mechanism.
    void attach_score_pointer(real*a) noexcept override { score = a; };

    /// Updates necessary after a reaction event was executed.
    void update_netw_stats() override;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> const std::string Fission<Ntw>::name {"fiss"};


template<typename Ntw> constexpr
auto Fission<Ntw>::
is_active( const std::unique_ptr<Reaction>& r ) noexcept
{
    return r->shortName == name;
};


template<typename Ntw>
void Fission<Ntw>::
initialize_dependencies( const vup<Reaction>& rc ) noexcept 
{
    for (const auto& o : rc) {
        if (
            Fusion<1,1,Ntw>::is_active(o) ||
            Fusion<1,2,Ntw>::is_active(o) ||
            Fission<Ntw>::is_active(o)
            )
            dependents.push_back(o.get());
    }
    set_prop();
    set_score();
}


template<typename Ntw> inline
void Fission<Ntw>::
set_score() noexcept
{
    *score = rate * netw.fis.get_prTotal();
}


template<typename Ntw> inline
void Fission<Ntw>::
set_prop() noexcept
{
    netw.fis.set_prop();
}


template<typename Ntw> inline
void Fission<Ntw>::
update_prop( szt c0, szt c1 ) noexcept
{
    netw.fis.update_prop(c0);
    if (c1 != c0 && c1 != huge<szt>)
        netw.fis.update_prop(c1);
}


template<typename Ntw>
void Fission<Ntw>::
update_netw_stats()
{
    for(auto& o : dependents) {
        o->update_prop(cc[0], cc[1]);
        o->set_score();
    }
}


template<typename Ntw>
void Fission<Ntw>::
fire() noexcept
{
    if constexpr (verbose) print(true);

    eventCount++;

    cc = netw.fis.fire();

    update_netw_stats();
}


template<typename Ntw>
void Fission<Ntw>::
print( const bool le ) const
{
    Reaction::print(false);
    msgr.print<false>(" score ", *score);
    msgr.print<false>(" eventCount ", eventCount);
    if (le) msgr.print("\n");
}

}  // namespace mitosim

#endif  // MITOSIM_FISSION_H
