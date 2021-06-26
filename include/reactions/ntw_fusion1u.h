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
 * @file ntw_fusion1u.h
 * @brief Contains class encapsulating slot on the graph which enables tip-to-loop fusion.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_NTW_FUSION1U_H
#define MITOSIM_NTW_FUSION1U_H

#include <array>
#include <vector>

#include "../fusion_candidates.h"
#include "definitions.h"

namespace mitosim {

template<typename> class Fusion1U;

/**
 * Reaction slot for fusion of a degree 1 node with a looped segment.
 * @details Network-specific reaction slot for fusion of a degree 1 node with
 * a looped segment.
 * @tparam Ntw Type of the network.
 */
template<typename Ntw>
class NtwFusion1L {

    friend Fusion1U<Ntw>;

public:

    explicit NtwFusion1L(Ntw&);  ///< The only constructor.

    /// Sets this reaction propensity for the whole network.
    auto set_prop() noexcept -> szt;

private:

    Ntw& host;  ///< ref: the host network for this reaction.

    // Convenience references to some of the host members.
    RandFactory&                          rnd;
    const std::vector<szt>&               mt11;
    const std::vector<std::array<szt,2>>& mt13;
    const std::vector<szt>&               mt22;

    FusionCandidatesXU cnd; ///< Node pairs suitable for this type of fusion.

    /// Populates the vector of node pairs suitable for this type of fusion.
    void populate() noexcept;

    /// Executes the reaction event.
    auto fire() noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion1L<Ntw>::
NtwFusion1L( Ntw& host )
    : host {host}
    , rnd {host.rnd}
    , mt11 {host.mt11}
    , mt13 {host.mt13}
    , mt22 {host.mt22}
{}

template<typename Ntw>
auto NtwFusion1L<Ntw>::
set_prop() noexcept -> szt
{
     populate();
     return cnd.size();
}

template<typename Ntw>
void NtwFusion1L<Ntw>::
populate() noexcept
{
    cnd.clear();
    for (const auto w2 : mt22) {
        for (const auto w1 : mt11)
            for (const auto e1 : {szt(1),szt(2)})
                cnd.add({w1,e1}, w2);        // e2 is 1 by convention

        for (const auto& we1 : mt13)
            cnd.add(we1, w2);
    }
}

template<typename Ntw>
auto NtwFusion1L<Ntw>::
fire() noexcept
{
    const auto r = rnd.uniform0(cnd.size());

    return host.fuse1L(cnd.u[r][0], cnd.u[r][1],
                       cnd.v[r]);
}


}  // namespace mitosim

#endif  // MITOSIM_NTW_FUSION1U_H
