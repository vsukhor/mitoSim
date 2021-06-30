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
 * @file ntw_fusion12.h
 * @brief Class encapsulating slot on the graph which enables tip-to-side fusion.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_NTW_FUSION12_H
#define MITOSIM_NTW_FUSION12_H

#include <array>
#include <vector>

#include "../fusion_candidates.h"
#include "definitions.h"

namespace mitosim {

template<typename> class Fusion12;

/**
 * Reaction slot for fusion of a degree 1 node with a degree 2 node,
 * Network-specific reaction slot for fusion of a degree 1 node
 * with a degree 2 node.
 * @tparam Ntw Type of the network.
 */
template<typename Ntw>
class NtwFusion12 {

public:

    friend Fusion12<Ntw>;

    explicit NtwFusion12(Ntw&); ///< Constructor.

    /// Sets this reaction propensity for the whole network.
    auto set_prop() noexcept -> szt;

    const FusionCandidatesXX& get_cnd() { return cnd; }

private:

    Ntw& host;  ///< ref: the host network for this reaction.
    
    // Convenience references to some of the host members.
    RandFactory& rnd;

    const typename Ntw::Reticulum&        mt;
    const std::vector<szt>&               mt11;
    const std::vector<std::array<szt,2>>& mt13;
    const std::vector<szt>&               mt22;
    const std::vector<szt>&               mt33;

    FusionCandidatesXX cnd;  ///< Node pairs suitable for this type of fusion.

    /// Populates the vector of node pairs suitable for this type of fusion.
    void populate() noexcept;

    /// Executes the raction event.
    auto fire() noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion12<Ntw>::
NtwFusion12( Ntw& host )
    : host {host}
    , rnd {host.rnd}
    , mt {host.mt}
    , mt11 {host.mt11}
    , mt13 {host.mt13}
    , mt22 {host.mt22}
    , mt33 {host.mt33}
{}

template<typename Ntw>
auto NtwFusion12<Ntw>::
set_prop() noexcept -> szt
{
     populate();
     return cnd.size();
}

template<typename Ntw>
void NtwFusion12<Ntw>::
populate() noexcept
{
    constexpr auto minLL = Structure<typename Ntw::ST>::minLoopLength;
    cnd.clear();
    for (const auto w1 : mt11)                           // 11 ends to ...
        for (const auto e1 : {szt(1),szt(2)}) {
            const std::array<szt,2> we1 {w1,e1};
            for (const auto w2 : mt11)                   // ... 11 bulk
                for (szt i=1; i<mt[w2].g.size(); i++) {
                    const auto skip = w1 == w2 && (
                                    (e1 == 1 && i < minLL) ||
                                    (e1 == 2 && mt[w2].g.size()-i < minLL));
                    if (!skip) {
                        cnd.add(we1, {w2,i});
                        cnd.add(we1, {w2,i});
                    }
                }
            for (const auto& we2 : mt13)                 // ... 13 bulk
                for (szt i=1; i<mt[we2[0]].g.size(); i++)
                    cnd.add(we1, {we2[0],i});

            for (const auto w2 : mt33)                    // ... 33 bulk
                for (szt i=1; i<mt[w2].g.size(); i++)
                    cnd.add(we1, {w2,i});

            for (const auto w2 : mt22)                    // ... 22 bulk
                for (szt i=1; i<mt[w2].g.size(); i++)
                    cnd.add(we1, {w2,i});
        }

    for (const auto& we1 : mt13) {                        // a free end of 13 to ...
        for (const auto w2 : mt11)                        // ... 11 bulk
            for (szt i=1; i<mt[w2].g.size(); i++)
                cnd.add(we1, {w2,i});

        for (const auto& we2 : mt13) {                    // ... 13 bulk
            for (szt i=1; i<mt[we2[0]].g.size(); i++) {
                const auto skip = we1[0] == we2[0] &&
                                 ((we1[1] == 1 && i < minLL) ||
                                  (we1[1] == 2 && mt[we2[0]].g.size()-i < minLL));
                if (!skip)
                    cnd.add(we1, {we2[0],i});
            }
        }
        for (const auto w2 : mt33)                       // ... 33 bulk
            for (szt i=1; i<mt[w2].g.size(); i++)
                cnd.add(we1, {w2,i});

        for (const auto w2 : mt22)                       // ... 22 bulk
            for (szt i=1; i<mt[w2].g.size(); i++)
                cnd.add(we1, {w2,i});
    }
}

template<typename Ntw>
auto NtwFusion12<Ntw>::
fire() noexcept
{
    const auto r = rnd.uniform0(cnd.size());

    return host.fuse12(cnd.u[r][0], cnd.u[r][1],
                       cnd.v[r][0], cnd.v[r][1]);
}

}  // namespace mitosim

#endif  // MITOSIM_NTW_FUSION12_H
