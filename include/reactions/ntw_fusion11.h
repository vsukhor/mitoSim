/* ==============================================================================
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

============================================================================== */

/**
* @file ntw_fusion11.h
* @brief Contains class encapsulating slot on the graph which enables tip-to-tip fusion.
* @author Valerii Sukhorukov
*/

#ifndef NTW_FUSION11_H
#define NTW_FUSION11_H

#include "utils/common/misc.h"
#include "fusion_candidates.h"

namespace MitoSim {

template<typename> class Fusion11;

/**
 * @brief Network-specific reaction slot for fusion of two nodes of degree 1.
 * @tparam Ntw Type of the network.
 */
template<typename Ntw>
class NtwFusion11 {

public:

    friend Fusion11<Ntw>;

    explicit NtwFusion11(Ntw&);	    ///< Constructor.

    /// Set this reaction propensity for the whole network.
    szt set_prop() noexcept;

private:

    Ntw& host;    ///< ref: the host network for this reaction.

    // Convenience references to some of the host members.
    RandFactory&    	    	    	    rnd;
    const std::vector<szt>&    	    	    mt11;
    const std::vector<std::array<szt,2>>&    mt13;

    FusionCandidatesXX    cnd; ///< node pairs suitable for this type of fusion.

    /// Populate the vector of node pairs suitable for this type of fusion.
    void populate() noexcept;

    /// Execute the raction event.
    auto fire() noexcept;
};
// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion11<Ntw>::
NtwFusion11( Ntw& host )
    : host {host}
    , rnd {host.rnd}
    , mt11 {host.mt11}
    , mt13 {host.mt13}
{}

template<typename Ntw>
szt NtwFusion11<Ntw>::
set_prop() noexcept
{
     populate();
     return cnd.size();
}

template<typename Ntw>
void NtwFusion11<Ntw>::
populate() noexcept
{
    constexpr auto minLL = Structure<typename Ntw::ST>::minLoopLength;

    cnd.clear();
    const auto mtn11 = mt11.size();
    for (szt i1=0; i1<mtn11; i1++) {    	    	    // 11 ends to ...
	    const auto w1 = mt11[i1];
	    if (host.mt[w1].g.size() >= minLL)	    	    // ... opposite end in the same segment
    	    cnd.add({w1,1}, {w1,2});

	    for (const auto e1 : {szt(1),szt(2)}) {
    	    for (szt i2=i1+1; i2<mtn11; i2++)    	    // ... other 11 segments (both ends to both ens)
	    	    for (const auto e2 : {szt(1),szt(2)})
    	    	    cnd.add({w1,e1}, {mt11[i2],e2});

    	    for (const auto& we2 : mt13)	    	    // ... free ends of 13
	    	    cnd.add({w1, e1}, we2);

    	    for (const auto& we2 : mt13)	    	    // ... free ends of 14
	    	    cnd.add({w1, e1}, we2);
	    }
    }
    const auto mtn13 = mt13.size();
    for (szt i1=0; i1<mtn13; i1++) 	    	    	    // free ends of 13 to ...
	    for (szt i2=i1+1; i2<mtn13; i2++)	    	    // ... free ends of other 13
    	    cnd.add(mt13[i1], mt13[i2]);
}

template<typename Ntw>
auto NtwFusion11<Ntw>::
fire() noexcept
{
    const auto r = rnd.uniform0(cnd.size());

    return host.fuse11(cnd.u[r][0], cnd.u[r][1], cnd.v[r][0], cnd.v[r][1]);
}

}    // namespace MitoSim

#endif // NTW_FUSION11_H
