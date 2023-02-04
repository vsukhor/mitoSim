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
 * @file structure.h
 * @brief Contains low-level network transformations implemented in class CoreTransformer.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_STRUCTURE_H
#define MITOSIM_STRUCTURE_H

#include <array>
#include <vector>

#include "definitions.h"

namespace mitosim {

/**
 * @brief The Structure class.
 * @details Encapsulates the major structure-related proterties, such a
 * collection of the Segments, but is unaware of any dynamics.
 * Forms base for clases adding the network reconfiguration dynamics.
 * @tparam Mt Type of the Edge forming the network.
*/
template<typename Mt>
class Structure {

public:

    using Reticulum = std::vector<Mt>;

    /// Mapping of the edge indexes to segment indexes.
    std::vector<szt> glm;
    /// Mapping of the edge indexes to element index inside segments.
    std::vector<szt> gla;

    /// The segments.
    Reticulum mt;

    /// Total number of nodes by node degree.
    std::array<szt,Mt::maxDegree> nn {};

    /// Actual number of segments.
    szt mtnum {};
    ///< Actual number of clusters.
    szt clnum {};
    ///<Current number of edges.
    szt mtmass {};

    /// Segment indices segregated into clusters: clmt - total.
    vec2<szt> clmt;
    /// Cluster sizes measured in edges.
    std::vector<szt> cls;

    // Segment indices of the corresponding end degrees.

    /// Indexes of disconnected segments not looped onto itself.
    std::vector<szt> mt11;
    /// Indexes of disconnected segments not looped onto itself.
    std::vector<szt> mtc11;

    /// Indexes of disconnected looped segments.
    std::vector<szt> mt22;
    /// Indexes of disconnected looped segments.
    std::vector<szt> mtc22;

    /// Indexes of segments between nodes of degree 3 and 3: all together.
    std::vector<szt> mt33;
    /// Indexes of segments between nodes of degree 3 and 3: sorted into clusters.
    vec2<szt> mtc33;

    /// {index,end} Pairs for segments between nodes of degs. 1 and 3 together.
    std::vector<std::array<szt,2>> mt13;
    /// {index,end} Pairs for segments between nodes of degs. 1 and 3: sorted into clusters.
    vec2<std::array<szt,2>> mtc13;

    /// Output message processor.
    Msgr& msgr;

    /// Minimal length of a segment that can bend into a cycle.
    static constexpr szt minLoopLength {2};

    /**
    * @brief Constructor.
    * @param msgr Output message processor.
    */
    explicit Structure(Msgr& msgr);

    /// Appends a disconnected segment to the reticulum.
    void add_disconnected_segment(szt segmass);

    /// Updates internal vectors.
    void update_structure() noexcept;

    /// Initializes or updates glm and gla vectors.
    void make_indma() noexcept;

    /// Populates 'mt??', 'mtc??', 'nn' and 'clmt' vectors
    void populate_cluster_vectors() noexcept;

    /**
     * Updates 'nn' for the specific node degree.
     * @tparam I Node degree to consider.
     */
    template<int I>
    void update_nn() noexcept;

    /// Updates 'nn' for all node degrese.
    void update_node_numbers() noexcept;

    /**
     * Print the network components using prefix specified.
     * @param tag Prefix.
     */
    void print_mitos(const std::string& tag) const;

    /**
     * Print the network components to a text-formatted stream.
     * @param ofs Output stream.
     */
    void print(std::ostream& ofs) const;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
Structure<Mt>::
Structure(
        Msgr& msgr
    )
    : msgr {msgr}
{}


template<typename Mt>
void Structure<Mt>::
add_disconnected_segment( const szt segmass )
{
    if (mt.empty())
        mt.emplace_back(msgr);  // Mock segment necessary for 1-based counting.

    mt.emplace_back(segmass, clnum, mtmass, msgr);
    mtnum++;
    clnum++;
    mtmass += segmass;
}


template<typename Mt> inline
void Structure<Mt>::
update_structure() noexcept
{
    make_indma();
    populate_cluster_vectors();
}

template<typename Mt> inline
void Structure<Mt>::
make_indma() noexcept
{
    cls.resize(clnum);
    std::fill(cls.begin(), cls.end(), 0);    // cluster size, # of edges
    for (szt j=1; j<=mtnum; j++)
        cls[mt[j].get_cl()] += mt[j].g.size();

    glm.resize(mtmass);
    gla.resize(mtmass);
    for (szt j=1; j<=mtnum; j++)
        for (szt k=0; k<mt[j].g.size(); k++) {
            const auto& g = mt[j].g[k];
            glm[g.get_ind()] = j;
            gla[g.get_ind()] = k;
        }
}

template<typename Mt>
void Structure<Mt>::
populate_cluster_vectors() noexcept
{
    mt11.clear();
    mtc11.resize(clnum);
    std::fill(mtc11.begin(), mtc11.end(), undefined<szt>);

    mt22.clear();
    mtc22.resize(clnum);
    std::fill(mtc22.begin(), mtc22.end(), undefined<szt>);

    mt33.clear();
    mtc33.resize(clnum);
    for (auto& o : mtc33) o.clear();

    mt13.clear();
    mtc13.resize(clnum);
    for (auto& o : mtc13) o.clear();

    nn = {{0}};
    clmt.resize(clnum);
    for (auto& o : clmt) o.clear();    // # of segments
    
    for (szt j=1; j<=mtnum; j++) {
        const auto& m = mt[j];
        clmt[m.get_cl()].push_back(j);    // mitochondria indexes clusterwise
        nn[1] += m.num_nodes(2);

        const auto e = m.has_one_free_end();
        if (e) {
            const szt oe {e == 1 ? 2UL : 1UL};
            nn[0]++;
            if (m.nn[oe] == 2) {
                const std::array<szt,2> je {j, e};
                mtc13[m.get_cl()].emplace_back(je);   // segment index, free end index
                mt13.emplace_back(je);
                nn[2]++;
            }
        }
        else if (m.nn[1] == 0 && m.nn[2] == 0) {
            mtc11[m.get_cl()] = j;  // having both ends free, it is a separate segment
            mt11.push_back(j);
            nn[0] += 2;
        }
        else if (m.is_cycle()) {
            mtc22[m.get_cl()] = j;
            mt22.push_back(j);  // it is a separate segment since it has two free ends
        }
        else if (m.nn[1] == 2 && m.nn[2] == 2) {
            mtc33[m.get_cl()].push_back(j);
            mt33.push_back(j);
            nn[2] += 2;
        }
        else {
            ; XASSERT(false,
                      "Error in populate_cluster_vectors: failed classification for "
                      +std::to_string(j)+"\n");
        }
    }
    nn[2] /= 3;
}

template<typename Mt>
template<int I>
void Structure<Mt>::
update_nn() noexcept
{
    auto count_nodes = [&](const szt deg) noexcept {
        szt k {};
        for (szt i=1; i<=mtnum; i++)
            k += mt[i].num_nodes(deg);
        return (deg == 3) ? k/3 : k;
    };

    nn[I-1] = count_nodes(I);
}

template<typename Mt>
void Structure<Mt>::
update_node_numbers() noexcept
{
    update_nn<1>();
    update_nn<2>();
    update_nn<3>();
}

template<typename Mt>
void Structure<Mt>::
print_mitos( const std::string& tag ) const
{
    for (szt j=1; j<=mtnum; j++)
        mt[j].print(j, tag);
    msgr.print("");
}
template<typename Mt>
void Structure<Mt>::
print( std::ostream& ofs ) const
{
    ofs << " X ";
    for (const auto o : nn)
        ofs << o << " ";
    ofs <<  "m11 " << mt11.size()
        << " m22 " << mt22.size()
        << " m33 " << mt33.size()
        << " m13 " << mt13.size()
        << " mtm " << mtmass
        << " mtn " << mtnum
        << " cln " << clnum;
}

}  // namespace mitosim

#endif  // MITOSIM_STRUCTURE_H
