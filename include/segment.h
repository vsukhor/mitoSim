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
 * @file segment.h
 * @brief Contains Segment class template and its specialization for graphs.
 * @details Only graphs of max degree 3 are considered.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_SEGMENT_H
#define MITOSIM_SEGMENT_H

#include <algorithm>
#include <array>
#include <vector>

#include "definitions.h"
#include "edge.h"

namespace mitosim {

/**
 * @brief Class template for the Network Segments.
 * @details Segment is a sequence of edges linked linearly (without branches).
 * Segment ends may form branching sites, where it is connected to other segments.
 * A segment not connected to other segments or
 * a complete collection of segments connected to each other form
 * a disconnected network component (aka 'cluster').
 * The class handles the tasks and properties specific to a single segment
 * and its relation to other network components.
 * @tparam _ Max node degree that the graph is able to handle.
 */
template<unsigned _>
class Segment {};

/**
 * @brief Segment class specification for max node degree equal to 3.
 * @details Segment is a sequence of edges linked linearly (without branches).
 * Segment ends may form branching sites, where it is connected to other segments.
 * A segment not connected to other segments or
 * a complete collection of segments connected to each other form
 * a disconnected network component (aka 'cluster').
 * The class handles the tasks and properties specific to a single segment
 * and its relation to other network components.
 */
template<>
class Segment<3> {

public:

    static constexpr szt numEnds {2};    ///< A segment has two ends.
    static constexpr szt maxDegree {3};  ///< Maximal node degree allowed.

    using EdgeT = Edge<maxDegree>;
    using thisT = Segment<maxDegree>;

    std::vector<EdgeT> g;  ///< The edges.

    /// Number of neighbours (for each of the two ends, counting from 1).
    std::array<szt,numEnds+1> nn {{}};

    std::array<std::vector<szt>,maxDegree> neig;  ///< Neighbour indexes.
    std::array<std::vector<szt>,maxDegree> neen;  ///< Neighbour ends.

    /**
     * @brief Constructor
     * @param msgr Output message processor.
     */
    explicit Segment(Msgr& msgr);

    /**
     * @brief Constructor.
     * @param cl Index of subnetwork to which the sebment belongs.
     * @param msgr Output message processor.
     */
    explicit Segment(Msgr& msgr,
                     szt cl);

    /**
     * @brief Constructor.
     * @param segmass Segment mass.
     * @param cl Index of subnetwork to which the sebment belongs.
     * @param ei Index of the last edge in this segment.
     * @param msgr Output message processor.
     */
    explicit Segment(
        szt segmass,
        szt cl,
        szt ei,
        Msgr& msgr );

    constexpr auto get_cl() const noexcept { return cl; }
    void set_cl( szt newcl ) noexcept { cl = newcl; }

    /// Reflect the vector containing the segment edges.
    void reflect_g();

    /**
     * @brief Change cluster index keeping the segment index unoltered.
     * @details Change disconnected network component-related indexes of the
     * segment edges, keeping the segment index unoltered.
     * @param newcl New disconnected component index.
     * @param initind Starting edge index in the current disconnected component.
     */
    auto set_gCl(szt newcl, szt initind) noexcept -> szt;

    /**
     * @brief Changecluster index.
     * @details Change disconnected network component-related indexes of the
     * segment edges, and the the segment itself.
     * @param newcl New disconnected component index.
     * @param initind Starting edge index in the current disconnected component.
     * @return The last edge index in the current disconnected network component.
     */
    auto setCl(szt newcl, szt initind) noexcept -> szt;

    /**
     * @brief Convert segment index to internal position.
     * @details Convert the segment end index of the boundary edge to internal
     * position in the segment.
     * @param e Segment end index.
     * @return Internal position.
     * @return The last edge index in the current disconnected component.
     */
    constexpr auto end2a(szt e) const noexcept -> decltype(g.size());

    /**
     * @brief Determine if the segment has one free end.
     * @return The end index if true.
     */
    constexpr auto has_one_free_end() const noexcept -> szt;

    /**
     * @brief Neigbour indexes at a segment end.
     * @param e Segment end.
     * @return the Neighbour index.
     */
    constexpr auto single_neig_index(szt e) const noexcept -> szt;

    /**
     * @brief Neigbour indexes at a segment end
     * @param e Segment end.
     * @return the Neighbour indexes.
     */
    auto double_neig_indexes(szt e) const -> std::vector<szt>;


    /// Report if the segment is looped onto itself.
    constexpr auto is_cycle() const noexcept -> bool;

    /**
     * @brief Report the number of nodes of a given degree.
     * @param deg Node degree.
     * @return the Number of nodes.
     */
    constexpr auto num_nodes(szt deg) const noexcept -> szt;

    /**Segment<3>::
     * @brief Report the segment length measured in edges.
     * @return the Segment length measured in edges.
     */
    auto length() const noexcept -> szt { return g.size(); }

    /**
     * @brief Set fission-specific factor for an end node.
     * @tparam E In-segment node posiiton.
     */
    template <unsigned E>
    constexpr auto set_end_fin() noexcept -> real;

    /**
     * @brief Set fission-specific factor for a bulk node.
     * @param a In-segment node posiiton.
     */
     auto set_bulk_fin(szt a) -> real;


    /// Print segment parameters.
    void print(
        szt w,
        const std::string& tag,
        szt at=huge<szt>
    ) const;


    /// Print segment parameters.
    void print(
        std::ostream& os,
        szt w,
        const std::string& tag,
        szt at=huge<szt>
    ) const;


    /**
     * @brief Write the segment to a binary file.
     * @param ofs std::ofstream to write to.
     */
    void write(std::ofstream& ofs) const;

protected:

    szt cl {};   ///< cluster index.
    
    Msgr& msgr;  ///< Output message processor.

    /**
     * @brief Insert an edge imediately after g[a] making it g[a+1].
     * @param a Position of the edge preceding the newly inserted
     * one relative to the segment end 1.
     * @param p Edge to be inserted.
     * @return The pointer to the newly inserted edge.
     */
    auto increment_length(long a, EdgeT p) -> EdgeT*;


    /// Initialise the neigbour vectors at both ends.
    void init_ends();
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

inline
Segment<3>::
Segment(Msgr& msgr)
    : msgr {msgr}
{
    init_ends();
}

inline
Segment<3>::
Segment(
    Msgr& msgr,
    const szt cl
)
    : cl {cl}
    , msgr {msgr}
{
    init_ends();
}

inline
Segment<3>::
Segment(
      const szt segmass,
      const szt cl,
      szt ei,
      Msgr& msgr     // var ref
)
    : Segment {msgr, cl}
{
    for (szt a=0; a<segmass; a++)
        increment_length(static_cast<long>(a-1),
                         EdgeT{ei++, a, cl});
}


inline
void Segment<3>::
init_ends()
{
    neig[1].resize(maxDegree);
    neig[2].resize(maxDegree);

    neen[1].resize(maxDegree);
    neen[2].resize(maxDegree);
}


// Inserts a particle imediately after g[a] making it g[a+1].
inline
auto Segment<3>::
increment_length(
    const long a,
    Segment<3>::EdgeT p
) -> EdgeT*
{
    g.insert(g.begin()+(a+1), std::move(p));
    return &g[static_cast<szt>(a) + 1];
}


inline
void Segment<3>::
reflect_g()
{
    std::reverse(g.begin(), g.end());
    for (auto& o : g)
        o.reflect();
}


inline
auto Segment<3>::
set_gCl(
    const szt newcl,
    const szt initind
) noexcept -> szt
{
    for (szt i=0; i<g.size(); i++) {
        g[i].cl = newcl;
        g[i].indcl = initind + i;
    }
    
    return initind + static_cast<szt>(g.size());
}


inline
auto Segment<3>::
setCl(
    const szt newcl,
    const szt initind
) noexcept -> szt
{
    cl = newcl;
    return set_gCl(newcl, initind);
}


constexpr
auto Segment<3>::
end2a( const szt e ) const noexcept -> decltype(g.size())
{
    XASSERT(e == 1 || e == 2, "Incorrect end index.");

    return (e == 1) ? 0 : g.size() - 1;
}


constexpr
auto Segment<3>::
has_one_free_end() const noexcept -> szt  // return the end index if true
{
    if (!nn[1] &&  nn[2]) return 1;
    if ( nn[1] && !nn[2]) return 2;
    /* else */            return 0;
}


constexpr
auto Segment<3>::
single_neig_index( const szt e ) const noexcept -> szt
{
    for (szt i=1; i<=nn[e]; i++) 
        if (neig[e][i])
            return i;
            
    return huge<szt>;
}


inline
auto Segment<3>::
double_neig_indexes( const szt e ) const -> std::vector<szt>
{
    XASSERT(nn[e] == 2,
            "Error in Segment::double_neig_indexes: nn[e] != 2 in cluster "+
            std::to_string(cl)+"\n");

    std::vector<szt> neigInds(nn[e]);
    for (szt j{}, i{1}; i<=nn[e]; i++)
        if (neig[e][i])
            neigInds[j++] = i;
    
    return neigInds;
}


constexpr
auto Segment<3>::
is_cycle() const noexcept -> bool
{
    return nn[1] == 1 && 
           nn[2] == 1 && 
           neig[1][single_neig_index(1)] == neig[2][single_neig_index(2)];
}


constexpr
auto Segment<3>::
num_nodes( const szt deg ) const noexcept -> szt // deg = 1, 2, 3
{                                                        
    if (deg == 1) {    // count nodes of degree 1
        if ( nn[1] &&  nn[2]) return 0;
        if (!nn[1] && !nn[2]) return 2;
        /* else */            return 1;
    }

    if (deg == 2)      // count nodes of degree 2
        return nn[1] && nn[2] && is_cycle()
               ? g.size()
               : g.size() - 1;

    if (deg == 3) {     // count nodes of degree 3
        if (nn[1] == 2 && nn[2] == 2) return 2;
        if (nn[1] == 2 || nn[2] == 2) return 1;
        if (nn[1] != 2 && nn[2] != 2) return 0;
    }

    msgr.exit("Error in Segment::num_nodes(). Not implemented for degree ", deg);
    return huge<szt>;
}


template <unsigned E> constexpr
auto Segment<3>::
set_end_fin() noexcept -> EdgeT::FinT
{
    static_assert(E == 1 || E == 2, "Incorrectt segment end index");

    auto& f = g[end2a(E)].fin;
    f[E-1] = nn[E] ? one<EdgeT::FinT>
                   : zero<EdgeT::FinT>;

    return f[E-1];
}


inline
auto Segment<3>::
set_bulk_fin( const szt a ) -> EdgeT::FinT
{
    XASSERT(a >= 0 || a <= g.size() - 1,
            std::string("Incorrectt segment edge index: ") + std::to_string(a));

    g[a].fin[1] = g[a+1].fin[0] = one<EdgeT::FinT>;

    return g[a].fin[1];
}


inline
void Segment<3>::
print( const szt w,
       const std::string& tag, 
       const szt at ) const
{
    if (msgr.so) print(*msgr.so, w, tag, at);
    if (msgr.sl) print(*msgr.sl, w, tag, at);
}


inline
void Segment<3>::
print( std::ostream& os, 
       const szt w,
       const std::string& tag, 
       const szt at) const
{
    os << "        " << tag << w;
    if (at == huge<szt>)
        os << "(of ";
    else
        os << "(at " << at << " of ";
    os << g.size() << ") [ ";    
    for (szt i=1; i<=nn[1]; i++) os << neig[1][i] << " ";
    os << "] { ";        
    for (szt i=1; i<=nn[1]; i++) os << neen[1][i] << " ";    
    os << "} [ ";
    for (szt i=1; i<=nn[2]; i++) os << neig[2][i] << " ";
    os << "] { ";    
    for (szt i=1; i<=nn[2]; i++) os << neen[2][i] << " ";    
    os << "} " << cl;    

    if constexpr (print_edges) {
        os << std::endl;
        for (szt i=0; i<g.size(); i++)
            g[i].print(os, i);
    }
    else
        os << " len " << g.size();

    os << std::endl;
}


inline
void Segment<3>::
write( std::ofstream& ofs ) const
{
    const auto len = static_cast<szt>(g.size());
    ofs.write(reinterpret_cast<const char*>(&len), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&cl), sizeof(szt));
    
    ofs.write(reinterpret_cast<const char*>(&nn[1]), sizeof(szt));
    
    for (szt j=1; j<=nn[1]; j++) {
        ofs.write(reinterpret_cast<const char*>(&neig[1][j]), sizeof(szt));
        ofs.write(reinterpret_cast<const char*>(&neen[1][j]), sizeof(szt));
    }
    
    ofs.write(reinterpret_cast<const char*>(&nn[2]), sizeof(int));
    
    for (szt j=1; j<=nn[2]; j++) {
        ofs.write(reinterpret_cast<const char*>(&neig[2][j]), sizeof(szt));
        ofs.write(reinterpret_cast<const char*>(&neen[2][j]), sizeof(szt));
    }
    for (const auto& a : g)
        a.write(ofs);
}

}  // namespace mitosim

#endif  // MITOSIM_SEGMENT_H
