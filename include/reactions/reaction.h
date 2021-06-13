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
* @file reaction.h
* @brief Base class for the simulated reactions.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_REACTION_H
#define MITOSIM_REACTION_H

#include <string>

#include "utils/common/misc.h"
#include "utils/msgr.h"
#include "utils/stochastic/gillespie.h"

namespace mitosim {

/// An abstract base class for all the reactions
class Reaction {

public:    // Only constant parameters are public.

    using szt = utils::common::szt;
    using ulong = utils::common::ulong;

    /// Index in Simulation::rc, i.e. index among all used and not used reactions.
    const szt ind {};
    /// Reaction rate constant.
    const real rate {};

    // Convenience references
    const ulong& it;    ///< ref: Internal network iteration counter.
    const real&  time;  ///< ref: Internal network time.

    const std::string shortName;  ///< Reaction name.
    const std::string srt;        ///< Reaction name abbreviation.

    /**
     * @brief Constructor.
     * @param msgr Output message processor.
     * @param ind reaction id.
     * @param rate Reaction rate constant.
     * @param it Iteration counter.
     * @param time Current time.
     * @param shortName Reaction name.
     * @param srt Reaction name abbreviated.
     */
    Reaction( utils::Msgr& msgr,
              const szt ind,
              const real rate,
              const ulong& it,              // const ref
              const real& time,             // const ref
              const std::string shortName,  // value + move
              const std::string srt         // value + move
        )
        : ind {ind}
        , rate {rate}
        , it {it}
        , time {time}
        , shortName {std::move(shortName)}
        , srt {std::move(srt)}
        , msgr {msgr}
    {}


    /// Virtual destructor.
    virtual ~Reaction() = default;


    /// Set the Gillespie score for this reaction.
    virtual void set_score() noexcept = 0;


    /// Return the Gillespie score for this reaction.
    virtual real get_score() const noexcept = 0;

    /**
     * @brief Update propensity for a pair of network components.
     * @param c1 Index of the 1st component to update.
     * @param c2 Index of the 2nd component to update.
     */
    virtual void update_prop(szt c1, szt c2) noexcept = 0;


    /// Execute the raction event.
    virtual void fire() noexcept = 0;


    /**
     * @brief Attach this score to the Gillespie mechanism.
     * @param a Placeholder in the Gillespie object responsible for this
     *        reaction score.
     */
    virtual void attach_score_pointer(real* a) noexcept = 0;


    /**
     * @brief Populate the vector of reactions that need a score update.
     * @details The update is performed after *this has fired
     *          and initializes the propensities and effective rate.
     * @param rc Vector of unique pointers to all reactions taking part in
     *        the simulation.
     */
    virtual void initialize_dependencies(
        const utils::common::vup<Reaction>& rc
    ) noexcept = 0;


    /**
     * @brief The number of times this reaction was fired.
     * @result The number of times this reaction was fired.
     */
    virtual szt event_count() const noexcept = 0;

    /**
     * @brief Print the parameters common to all reactions.
     * @param le True if new line after the output.
     */
    virtual void print( const bool le ) const
    {
        msgr.print<false>(" it ", it);
        msgr.print<false>(" srt ", srt);
        msgr.print<false>(" rate ", rate);
        if (le) msgr.print("\n");
    }

protected:

    utils::Msgr& msgr;  ///< ref: Output message processor.

    /** All necessary updates after the given reaction event was executed.
     * @details Pure virtual function: Network and reaction updates necessary
     * after the given reaction event was executed.
     */
    virtual void update_netw_stats() = 0;

};

}  // namespace mitosim

#endif  // MITOSIM_REACTION_H
