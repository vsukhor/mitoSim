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
* @file fusion.h
* @brief Contains base class for fusions reactions.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_FUSION_h
#define MITOSIM_FUSION_h

#include <string>
#include <vector>

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "reaction.h"

namespace mitosim {

using utils::common::uint;

/**
 * @brief Base class for fusion reaction classes.
 * @tparam D1 node degree of the 1st fusion participant
 * @tparam D2 node degree of the 2nd fusion participant
 * @tparam Ntw type of the network
 */
template<uint D1,
         uint D2,
         typename Ntw>
class Fusion
    : public Reaction {

    friend utils::common::Gillespie<Reaction,RandFactory>;

public:

    /**
     * @brief Constructor.
     * @param msgr Output message processor.
     * @param ind reaction id
     * @param netw The network object.
     * @param rate rate constant
     * @param rnd Random number factory.
     * @param it iteration counter
     * @param time current time
     * @param srt reaction name literal
     */
    Fusion( utils::common::Msgr& msgr,
            const szt ind,
            Ntw& netw,
            const real rate,
            RandFactory& rnd,
            const ulong& it,        // const ref
            const real& time,        // const ref
            const std::string& srt
        )
        : Reaction {msgr, ind, rate, it, time, "fusion", srt}
        , netw {netw}
        , rnd {rnd}
    {}

    /**
     * @brief Activity status of the reaction.
     * @return True if the reaction is used in the current simulation session.
     * @param r Pointer to the reaction.
     */
    static constexpr auto is_active(const std::unique_ptr<Reaction>& r);

    /**
     * @brief Populate the vector of reactions that need a score update.
     * @details The update is performed after *this has fired
     *          and initializes the propensities and effective rate.
     * @param rc Vector of unique pointers to all reactions taking part in
     *        the simulation.
     */
    void initialize_dependencies(const vup<Reaction>& rc) noexcept override;

    /**
    * @brief Print the parameters.
    * @param le True if new line after the output.
    */
    void print(bool le) const override;

protected:

    using Reaction::msgr;

    // Convenience references.
    Ntw&         netw;        ///< ref: The host network for this reaction.
    RandFactory& rnd;        ///< ref: Random number factory.

    std::array<szt,2> cc;   ///< Indices of the clusters.

    void update_netw_stats() override;

private:

    using Reaction::it;
    using Reaction::set_score;

    /// Reactions that need a score update after *this has fired.
    std::vector<Reaction*>  dependents;

    /// Sets this reaction propensity for the whole network.
    virtual void set_prop() noexcept = 0;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<uint D1, uint D2, typename Ntw> constexpr
auto Fusion<D1,D2,Ntw>::
is_active( const std::unique_ptr<Reaction>& r )
{
    return  r->srt == "fu1L" ||
            r->srt == "fu11" ||
            r->srt == "fu12";
}


template<uint D1, uint D2, typename Ntw>
void Fusion<D1,D2,Ntw>::
initialize_dependencies( const vup<Reaction>& rc ) noexcept
{
    for (const auto& o : rc)
        if (
            Fusion<D1,D2,Ntw>::is_active(o) ||
            Fission<Ntw>::is_active(o)
            )
            dependents.push_back(o.get());

    set_prop();
    set_score();
}


template<uint D1, uint D2, typename Ntw>
void Fusion<D1,D2,Ntw>::
update_netw_stats()
{
    netw.update_books();
    
    for (auto& o : dependents) {
        o->update_prop(cc[0], cc[1]);
        o->set_score();
    }
}


template<uint D1, uint D2, typename Ntw>
void Fusion<D1,D2,Ntw>::
print( const bool le ) const
{
    Reaction::print(false);
    msgr.print<false>(" deg1 ", D1);
    msgr.print<false>(" deg2 ", D2);
    if (le) msgr.print("");
}

}  // namespace mitosim

#endif  // MITOSIM_FUSION_h
