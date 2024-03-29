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
 * @file ntw_fission.h
 * @brief Contains class encapsulating slot on the graph which enables fission.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_NTW_FISSION_H
#define MITOSIM_NTW_FISSION_H

#include <array>
#include <vector>

#include "definitions.h"

namespace mitosim {

template<typename> class Fission;

/**
 * Base class for network-specific fission reaction slots.
 * @tparam Ntw Type of the network.
 */
template<typename Ntw>
class NtwFission {

public:

    using Prop = typename Ntw::ST::EdgeT::FinT;

    friend Fission<Ntw>;

    explicit NtwFission(Ntw& host);  ///< Constructor.

    /// Sets this reaction propensity for the whole network.
    Prop set_prop()  noexcept;

    /**
     * Update this reaction propensity for the whole network.
     * This is done after updating it for the cluster indexed.
     * @param c Cluster index that triggers the update.
     */
    void update_prop(szt c) noexcept;

    /// prTotal getter.
    constexpr auto get_prTotal() const noexcept { return prTotal; }

private:

    Ntw& host;  ///< ref: The host network for this reaction.

    // Convenience references to some of the host members
    typename Ntw::Reticulum& mt;     ///< ref: The segments.
    const szt&               clnum;  ///< ref: Current number fo clusters.

    // Propensities:
    std::vector<Prop> pr;   ///< Propensities per cluster.
    Prop prTotal {};  ///< Total propensity.

    /**
     * Set this reaction propensity for the indexed cluster.
     * @note Does not update the whole network propensity.
     * @param c Index of the cluster that is updateed
     */
    void set_prop(szt c) noexcept;

    /// Executes the raction event.
    auto fire() noexcept;

    /**
     * @brief Finds sa random node from those suitable for this reaction.
     * @param w Index of random segment.
     * @param a Random position inside the segment.
     */
    bool find_random_node(szt& w, szt& a) const noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFission<Ntw>::
NtwFission( Ntw& host )
    : host {host}
    , mt {host.mt}
    , clnum {host.clnum}
{}

template<typename Ntw>
auto NtwFission<Ntw>::
set_prop() noexcept -> Prop
{
    pr.resize(clnum);
    for (szt ic=0; ic<clnum; ic++)
        set_prop(ic);

    return std::accumulate(pr.begin(), pr.end(), zero<Prop>);
}

template<typename Ntw>
void NtwFission<Ntw>::
set_prop( const szt ic ) noexcept
{
    pr[ic] = zero<Prop>;
    for (const auto w : host.clmt[ic]) {
        pr[ic] += mt[w].template set_end_fin<1>() +
                  mt[w].template set_end_fin<2>();
        for (szt a=0; a<mt[w].g.size()-1; a++)
            pr[ic] += 2UL * mt[w].set_bulk_fin(a);
    }
}

template<typename Ntw>
void NtwFission<Ntw>::
update_prop( const szt c ) noexcept
{
    // We assume ncremental clnum changes.
    if (pr.size() > clnum) {
        pr.resize(clnum);
    }
    else if (pr.size() < clnum)
        pr.resize(clnum);

    if (c < clnum)
        set_prop(c);

    prTotal = std::accumulate(pr.begin(), pr.end(), zero<Prop>);
}

template<typename Ntw>
auto NtwFission<Ntw>::
fire() noexcept
{
    auto w = undefined<szt>;
    auto a = undefined<szt>;

    find_random_node(w, a);

    return host.fiss(w, a);
}

template<typename Ntw>
bool NtwFission<Ntw>::
find_random_node( szt& w, szt& a ) const noexcept
{
    auto k = host.rnd.uniform0(prTotal);
    typename Ntw::ST::EdgeT::FinT ksum {};
    for (w=1; w<=host.mtnum; w++) {
        const auto& g = mt[w].g;
        a = 0;
        ksum += g[a].get_fin(0);
        if (k <= ksum)
            return true;
        for (; a<g.size()-1;) {
            ksum += g[a].get_fin(1);
            a++;
            ksum += g[a].get_fin(0);
            if (k <= ksum)
                return true;
        }
        ksum += g[a].get_fin(1);
        if (k <= ksum) {
            a++;
            return true;
        }
    }
    return false;
}

}  // namespace mitosim

#endif  // MITOSIM_NTW_FISSION_H
