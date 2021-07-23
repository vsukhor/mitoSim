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
 * @file ability_for_fission.h
 * @brief Contains class responsible for the graph fillsin capacity.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_ABILITY_FOR_FISSION_H
#define MITOSIM_ABILITY_FOR_FISSION_H

#include <algorithm>
#include <vector>

#include "core_transformer.h"
#include "definitions.h"

namespace mitosim {

template<typename> class Fission;

/**
 * @brief The AbilityForFission class template.
 * @details Adds node type-specific fission capability
 * and updates the network for it.
 * @tparam Mt Type of the Edge forming the network.
 */
template<typename Mt>
class AbilityForFission
    : public CoreTransformer<Mt> {

    using thisT = AbilityForFission<Mt>;
    using ThisFission = Fission<thisT>;

protected:

    using Structure<Mt>::msgr;
    using Structure<Mt>::mt;
    using Structure<Mt>::mtnum;
    using Structure<Mt>::clnum;
    using Structure<Mt>::glm;
    using CoreTransformer<Mt>::copy_neigs;
    using CoreTransformer<Mt>::update_neigs;
    using CoreTransformer<Mt>::fuse_antiparallel;
    using CoreTransformer<Mt>::fuse_parallel;
    using CoreTransformer<Mt>::update_gIndcl;

    using Structure<Mt>::update_structure;

public:

    /**
     * @brief Constructor.
     * @param msgr Output message processor.
     */
    explicit AbilityForFission(Msgr& msgr);

    /**
     * @brief Perform fission of a segment.
     * @param w Segment index.
     * @param a Division positon inside the segment.
     */
    auto fiss(szt w, szt a) -> std::array<szt,2>;

    /**
     * @brief Divide the segment at a node of degree 2.
     * @param w Global segment index.
     * @param a The node position inside the segment.
     */
    auto fiss2(szt w, szt a) -> std::array<szt,2>;

    /**
     * @brief Divide the segment at a node of degree 3.
     * @param w Global segment index.
     * @param end The the segment end index at which the cut is done.
     */
    auto fiss3(szt w, szt end) -> std::array<szt,2>;

private:

    /// Auxiliary, indicating visited status of segments during the graph search.
    std::vector<szt>  vis {};

    /**
     * @brief Update network component for changes resulting from its division.
     * @return A flag indicating if the division produces a pair
     * of disconnected components.
     * @param w Segment index.
     * @param e Segment end.
     */
    bool update_cl_fiss(szt w, szt e) noexcept;

    /**
     * @brief Depth-first search of the network graph.
     * @param w1 Initial segment index.
     * @param e1 Initial segment end.
     * @param w2 Final segment index.
     * @param e2 Final segment end.
     * @return True if there is a connection between {w1,e1} and {w2,e2}.
     */
    bool dfs(szt w1, szt e1, szt w2, szt e2);
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFission<Mt>::
AbilityForFission(Msgr& msgr)
    : CoreTransformer<Mt> {msgr}
{}


template<typename Mt>
bool AbilityForFission<Mt>::
update_cl_fiss( const szt w,
                const szt e ) noexcept
{
    vis.resize(mtnum + 2);
    // Nothing is visited in the beginning of the search.
    std::fill(vis.begin() + 1,
              vis.begin() + static_cast<long>(mtnum) + 1,
              0);

    const auto oe = e == 1 ? 2UL : 1UL;

    const bool is_cycle = dfs(w, e, w, oe);
    if (!is_cycle) {
        clnum++;
        szt clind {};
        for (szt i=1; i<=mtnum; i++)
            if (vis[i])
                clind = mt[i].setCl(clnum - 1, clind);
    }
    return is_cycle;
}


template<typename Mt>
bool AbilityForFission<Mt>::
dfs( const szt w1, const szt e1,
     const szt w2, const szt e2 )
{
    for (szt i=1; i<=mt[w1].nn[e1]; i++) {
        const auto cn = mt[w1].neig[e1][i];
        const auto ce = mt[w1].neen[e1][i];
        if (cn == w2) {
            if (ce == e2)
                return true;
        }
        else {
            if (!vis[cn]) {
                vis[cn] = 1;
                const auto ne = ce == 1 ? 2UL : 1UL;
                if (dfs(cn, ne, w2, e2))
                    return true;
            }
        }    
    }
    return false;
}


// 'a' is counted from 1 and is the last element to remain in the old segment.
template<typename Mt>
auto AbilityForFission<Mt>::
fiss( const szt w,
      const szt a ) -> std::array<szt,2>
{
    // node 2 -> nodes 1+1; cuts between g[a-1] and g[a]
    if (a && a < mt[w].g.size())
        return fiss2(w, a);

    // node 3 -> nodes 2+1 and node 2 -> nodes 1+1 in pure loops
    if (!a && mt[w].nn[1] <= 2)
        return fiss3(w, 1);
    if (a == mt[w].g.size() && mt[w].nn[2] <= 2)
        return fiss3(w, 2);

    msgr.exit("ERROR: Attempt of an unpropriate fission!");

    return {huge<szt>, huge<szt>};
}


// divides at a node of degree 2
// 'w' is a global segment index
// 'a' is the node position inside the segment 'w'
// 'a' is counted from 1 and is the last element to remain in the original segment
template<typename Mt>
auto AbilityForFission<Mt>::
fiss2( const szt w,
       const szt a ) -> std::array<szt,2>
{
    if constexpr (verbose)
        mt[w].print(w, "fission2:  ", a);     // cuts between g[a-1] and g[a]

    XASSERT(a && a < mt[w].g.size(), "Error: fiss2 at the segment border.");

    [[maybe_unused]] const auto clini = mt[w].get_cl();

    const auto ind1 = mt[w].g[a-1].get_ind();
    const auto ind2 = mt[w].g[a].get_ind();

    bool inCycle {};
    mt[w].nn[2] ? (inCycle = update_cl_fiss(w, 2))
                : clnum++;

    mt.emplace_back(msgr);
    ++mtnum;

    std::move(mt[w].g.begin() + static_cast<long>(a),
              mt[w].g.end(),
              std::back_inserter(mt[mtnum].g));
    mt[w].g.erase(mt[w].g.begin() + static_cast<long>(a),
                  mt[w].g.end());

    mt[mtnum].nn[1] = 0;

    this->copy_neigs(w, 2, mtnum, 2);

    mt[mtnum].set_cl(inCycle
                     ? mt[w].get_cl()
                     : clnum - 1);

    if (!inCycle) {
        // Renumber Edge::indcl of the remaining part of the original cluster.
        update_gIndcl(mt[w].get_cl());
        // Renumber Edge::indcl of the newly formed cluster.
        update_gIndcl(clnum-1);
    }

    mt[w].nn[2] = 0;

    [[maybe_unused]] bool isSelfLooped {};
    if (1 == mt[w].nn[1] &&
        1 == mt[mtnum].nn[2] && 
        w == mt[mtnum].neig[2][mt[mtnum].single_neig_index(2)] && 
        mtnum == mt[w].neig[1][mt[w].single_neig_index(1)]) {
        // Before fission 'w' was looped into itself.
        isSelfLooped = true;
        if constexpr (verbose) {
            mt[w].print(w, "transiently producing ");
            mt[mtnum].print(mtnum, "      and ");
            std::cout << std::endl;
        }
        update_neigs(w, 1, 1, 1);
        fuse_parallel(w, mtnum);
    }

    update_structure();
    const auto w1 = glm[ind1];
    const auto w2 = glm[ind2];

    XASSERT(mt[w1].get_cl() == clini ||
            mt[w2].get_cl() == clini,
            "Error in fiss2: mt[w1].cl != clini && mt[w2].cl != clini\n");

    if constexpr (verbose) {
        mt[w1].print(w1, "producing ");
        if (isSelfLooped)
            msgr.print("from a segment looped into itself");
        else mt[w2].print(w2, "      and ");
        std::cout << std::endl;
    }

    return {mt[w].get_cl(), mt[mtnum].get_cl()};
}

// divides at a node of degree 3
// 'w' is a global segment index
// 'at' is the node position inside the segment 'w'
// 'at' is counted from 1 and is the last element to remain in the original segment
template<typename Mt>
auto AbilityForFission<Mt>::
fiss3( const szt w,
       const szt end ) -> std::array<szt,2>
{
    if constexpr (verbose)
        mt[w].print(w, "fission3:  ", end);

    const auto clini = mt[w].get_cl();
    bool f {};
    std::array<szt,2> n;
    std::array<szt,2> e;
    auto ind1 = huge<szt>;
    auto ind2 = huge<szt>;

    if (end == 1) {
        ind1 = mt[w].g.front().get_ind();
        ind2 = mt[mt[w].neig[1][1]]
                .g[mt[mt[w].neig[1][1]].end2a(mt[w].neen[1][1])]
                .get_ind();

        if (mt[w].nn[1] == 2) {
            const auto ninds = mt[w].double_neig_indexes(1);
            f = true;
            for (szt j=0; j<2; j++) {
                n[j] = mt[w].neig[1][ninds[j]];
                e[j] = mt[w].neen[1][ninds[j]];
            }
        }
        else if (mt[w].nn[1] == 1)
            n[0] = mt[w].neig[1][1];

        // If not a cycle, this increments clnum and forms a new cluster
        // from w's end 1 neigs and beyond (excluding w itself).
        const auto inCycle = update_cl_fiss(w, 1);
        if (!inCycle)
            // Renumber Edge::indcl of the remaining part of the original cluster.
            update_gIndcl(clini);

        update_neigs(w, 1, 1, mt[w].nn[1]);

        if (f && n[0] != n[1]) {
            if (mt[n[0]].nn[e[0]] == 1 && 
                mt[n[1]].nn[e[1]] == 1) {

                const auto indneig0 = mt[n[0]].single_neig_index(e[0]);
                const auto indneig1 = mt[n[1]].single_neig_index(e[1]);

                if (mt[n[0]].neig[e[0]][indneig0] == n[1] && 
                    mt[n[0]].neen[e[0]][indneig0] == e[1] && 
                    mt[n[1]].neig[e[1]][indneig1] == n[0] && 
                    mt[n[1]].neen[e[1]][indneig1] == e[0]) {

                    update_neigs(n[0], e[0], 1, 1);

                    if (e[0] == e[1])
                        fuse_antiparallel(e[0], n[0], n[1]);
                    else
                        (e[0] == 1 && e[1] == 2) ? fuse_parallel(n[0], n[1])
                                                 : fuse_parallel(n[1], n[0]);
                }
            }
            else if (mt[n[0]].nn[e[0]] == 0 && 
                     mt[n[1]].nn[e[1]] == 1) {
                mt[w].print(w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1");
            }
            else if (mt[n[0]].nn[e[0]] == 1 && 
                     mt[n[1]].nn[e[1]] == 0) {
                mt[w].print(w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0");
            }
            else if (mt[n[0]].nn[e[0]] == 0 && 
                     mt[n[1]].nn[e[1]] == 0) {
                mt[w].print(w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0");
            }
        }
    }
    else if (end == 2) {
        ind1 = mt[w].g.back().get_ind();
        ind2 = mt[mt[w].neig[2][1]]
                .g[mt[mt[w].neig[2][1]].end2a(mt[w].neen[2][1])]
                .get_ind();
        if (mt[w].nn[2] == 2) {
            const auto ninds = mt[w].double_neig_indexes(2);
            f = true;
            for (szt j=0; j<2; j++) {
                n[j] = mt[w].neig[2][ninds[j]];
                e[j] = mt[w].neen[2][ninds[j]];
            }
        }
        else if (mt[w].nn[2] == 1)
            n[0] = mt[w].neig[2][1];

        // If not a cycle, this increments clnum and forms a new cluster
        // from w's end 2 neigs and beyond (excluding w itself).
        const auto inCycle = update_cl_fiss(w, 2);
        if (!inCycle)
            // Renumber Edge::indcl of the remaining part of the original cluster.
            update_gIndcl(clini);

        update_neigs(w, 2, 1, mt[w].nn[2]);

        if (f && n[0] != n[1]) {
            if (mt[n[0]].nn[e[0]] == 1 && 
                mt[n[1]].nn[e[1]] == 1 ) {
                const auto indneig0 = mt[n[0]].single_neig_index(e[0]);
                const auto indneig1 = mt[n[1]].single_neig_index(e[1]);
                if (mt[n[0]].neig[e[0]][indneig0] == n[1] && 
                    mt[n[0]].neen[e[0]][indneig0] == e[1] && 
                    mt[n[1]].neig[e[1]][indneig1] == n[0] && 
                    mt[n[1]].neen[e[1]][indneig1] == e[0]) {
                    update_neigs(n[0], e[0], 1, 1);
                    if (e[0] == e[1])
                        fuse_antiparallel(e[0], n[0], n[1]);
                    else
                        (e[0] == 1 && e[1] == 2) ? fuse_parallel(n[0], n[1])
                                                 : fuse_parallel(n[1], n[0]);
                }
            }
            else if (mt[n[0]].nn[e[0]] == 0 && 
                     mt[n[1]].nn[e[1]] == 1) {
                mt[w].print(w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1");
            }
            else if (mt[n[0]].nn[e[0]] == 1 && 
                     mt[n[1]].nn[e[1]] == 0) {
                mt[w].print(w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0");
            }
            else if (mt[n[0]].nn[e[0]] == 0 && 
                     mt[n[1]].nn[e[1]] == 0) {
                mt[w].print( w, "w ");
                mt[n[0]].print(n[0], "n[0] ");
                mt[n[1]].print(n[1], "n[1] ");
                msgr.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0");
            }
        }
    }

    update_structure();

    const auto w1 = glm[ind1];
    const auto w2 = glm[ind2];
    XASSERT(mt[w1].get_cl() == clini ||
            mt[w2].get_cl() == clini,
            "Error in fiss3: mt[w1].cl != clini && mt[w2].cl != clini\n");

    if constexpr (verbose) {
        mt[w1].print(w1, "producing ");
        if (w1 != w2)
            mt[w2].print(w2, "      and ");
        std::cout << std::endl;
    }

    return {mt[w1].get_cl(),
            mt[w2].get_cl()};
}

}  // namespace mitosim

#endif  // MITOSIM_ABILITY_FOR_FISSION_H
