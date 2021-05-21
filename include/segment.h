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
* @file segment.h
* @brief Contains Segment class template and its specialization for graphs.
* @details Only graphs of max degree 3 are considered.
* @author Valerii Sukhorukov
*/

#ifndef SEGMENT_H
#define SEGMENT_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "edge.h"

namespace MitoSim {

/**
 * @brief Class template for the Network Segments.
 * @details Segment is a sequence of edges linked linearly (without branches).
 * Segment ends may form branching sites, where it is connected to other segments.
 * A segment not connected to other segments or
 * a complete collection of segments connected to each other form
 * a disconnected network component (aka 'cluster').
 * The class handles the tasks and properties specific to a single segment
 * and its relation to other network components.
 * @tparam MAXDEG Max node degree that the graph is able to handle.
 */
template<szt MAXDEG>
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

    static constexpr szt maxDeg {3};    ///< Maximal node degree allowed.

    using EdgeT = Edge<maxDeg>;
    using thisT = Segment<maxDeg>;

    Msgr&	    	    	    	    msgr;	    ///< Output message processor.
    std::vector<EdgeT>    	    	    g;    	    ///< the edges.
    szt    	    	    	    	    cl {};	    ///< cluster index.
    std::array<szt,maxDeg>	    	    nn {{}};    ///< number of neighbours.
    std::array<std::vector<szt>,maxDeg> neig;	    ///< neighbours.
    std::array<std::vector<szt>,maxDeg> neen;	    ///< neighbour ends.

    /**
     * @brief Constructor
     * @param msgr Output message processor.
     */
    explicit Segment(Msgr& msgr);

    /**
     * @brief Constructor.
     * @param segmass Segment mass.
     * @param cl Index of subnetwork to which the sebment belongs.
     * @param mtmass Total mass of the network.
     * @param ei Index of the last edge in this segment.
     * @param msgr Output message processor.
     */
    explicit Segment(
	      const szt segmass,
	      const szt cl,
	      szt& mtmass, 
	      szt& ei, 
	      Msgr& msgr );

    /// Reflect the vector containing the segment edges.
    void reflect_g();

    /**
     * @brief Change cluster index keeping the segment index unoltered.
     * @details Change disconnected network component-related indexes of the
                segment edges, keeping the segment index unoltered.
     * @param newcl New disconnected component index.
     * @param initind Starting edge index in the current disconnected component.
     */
    szt set_gCl(
        const szt newcl,
        const szt initind
    );

    /**
     * @brief Changecluster index.
     * @details Change disconnected network component-related indexes of the
     *          segment edges, and the the segment itself.
     * @param newcl New disconnected component index.
     * @param initind Starting edge index in the current disconnected component.
     * @return The last edge index in the current disconnected network component.
     */
    szt setCl(
        const szt newcl,
        const szt initind
    );

    /**
     * @brief Convert segment index to internal position.
     * @details Convert the segment end index of the boundary edge to internal
     *          position in the segment.
     * @param e Segment end index.
     * @return Internal position.
     * @return The last edge index in the current disconnected component.
     */
    constexpr szt end2a(const szt e) const;

    /**
     * @brief Determine if the segment has one free end.
      * @return The end index if true.
     */
    constexpr szt has_one_free_end() const;

    /**
     * @brief Neigbour indexes at a segment end.
     * @param e Segment end.
      * @return the Neighbour index.
     */
    szt single_neig_index(const szt e) const;

    /**
     * @brief Neigbour indexes at a segment end
     * @param e Segment end.
      * @return the Neighbour indexes.
     */
    std::vector<szt> double_neig_indexes(const szt e) const;


    /// Report if the segment is looped onto itself.
    constexpr bool is_cycle() const;

    /**
     * @brief Report the number of nodes of a given degree.
     * @param deg Node degree.
     * @return the Number of nodes.
     */
    szt num_nodes(const szt deg) const;

    /**
     * @brief Report the segment length measured in edges.
     * @return the Segment length measured in edges.
     */
    szt length() const noexcept { return g.size(); }

    /**
     * @brief Set fission-specific factor for an end node.
     * @param a In-segment node posiiton.
     */
    ulong set_end_fin(const szt a);

    /**
     * @brief Set fission-specific factor for a bulk node.
     * @param a In-segment node posiiton.
     */
    ulong set_bulk_fin(const szt a);


    /// Print segment parameters.
    void print(
        const szt w,
        const std::string& tag,
        const szt at=huge<szt>
    ) const;


    /// Print segment parameters.
    void print(
        std::ostream& os,
        const szt w,
        const std::string& tag,
        const szt at=huge<szt>
    ) const;


    /**
     * @brief Write the segment to a binary file.
     * @param ofs std::ofstream to write to.
     */
    void write(std::ofstream& ofs) const;

private:

    /**
     * @brief Insert an edge imediately after g[a] making it g[a+1].
     * @param a Position of the edge preceding the newly inserted
     *        one relative to the segment end 1.
     * @param p Edge to be inserted.
     * @return The pointer to the newly inserted edge.
     */
    EdgeT* increment_length( const long a, EdgeT p );


    /// Initialise the neigbour vectors at both ends.
    void init_ends();
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Segment<3>::
Segment( Msgr& msgr )
    : msgr {msgr}
{
    init_ends();
}


Segment<3>::
Segment(
      const szt segmass,
      const szt cl,
      szt& mtmass, 	    // var ref
      szt& ei,     	    // var ref
      Msgr& msgr	    // var ref
    )
    : msgr {msgr}
    , cl {cl}
{
    init_ends();

    for (szt a=0; a<segmass; a++)
	    increment_length(long(a-1), EdgeT{ei++, a, cl});

    mtmass += segmass;
}


inline
void Segment<3>::
init_ends()
{
    neig[1].resize(maxDeg);
    neig[2].resize(maxDeg);

    neen[1].resize(maxDeg);
    neen[2].resize(maxDeg);
}


// Inserts a particle imediately after g[a] making it g[a+1].
inline
typename Segment<3>::EdgeT* Segment<3>::
increment_length(
    const long a,
    Segment<3>::EdgeT p
)
{
    g.insert(g.begin()+(a+1), std::move(p));
    return &g[a+1];
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
szt Segment<3>::
set_gCl(
    const szt newcl,
    const szt initind
)
{
    for (szt i=0; i<g.size(); i++) {
	    g[i].cl = newcl;
	    g[i].indcl = initind + i;
    }
    
    return initind + (szt)g.size();
}


inline
szt Segment<3>::
setCl(
    const szt newcl,
    const szt initind
)
{
    cl = newcl;
    return set_gCl(newcl, initind);
}


constexpr
szt Segment<3>::
end2a( const szt e ) const
{
    return (e == 1) ? 0 : (szt)g.size()-1;
}


constexpr
szt Segment<3>::
has_one_free_end() const    // return the end index if true
{
    if (    !nn[1] &&  nn[2]) return 1;
    else if (nn[1] && !nn[2]) return 2;
    else	    	    	  return 0;
}


inline
szt Segment<3>::
single_neig_index( const szt e ) const
{
    for (szt i=1; i<=nn[e]; i++) 
	    if (neig[e][i])
    	    return i;
    	    
    return huge<szt>;
}


inline
std::vector<szt> Segment<3>::
double_neig_indexes( const szt e ) const
{
    XASSERT(nn[e] == 2,
            "Error in Mito::double_neig_indexes: nn[e] != 2 in cluster "+
            STR(cl)+"\n");

    std::vector<szt> neigInds(nn[e]);
    szt j {};
    for (szt i=1; i<=nn[e]; i++)
	    if (neig[e][i])
    	    neigInds[j++] = i;
    
    return neigInds;
}


constexpr
bool Segment<3>::
is_cycle() const
{
    return nn[1] == 1 && 
	       nn[2] == 1 && 
	       neig[1][single_neig_index(1)] == neig[2][single_neig_index(2)];
}


szt Segment<3>::
num_nodes( const szt deg ) const	    	    // deg = 1, 2, 3
{	    	    	    	    	    	    	    
    if (deg == 1) {	    	    	// count nodes of degree 1
	    if (      nn[1] &&  nn[2]) return 0;
	    else if (!nn[1] && !nn[2]) return 2;
	    else    	    	       return 1;
    }
    else if (deg == 2)      // count nodes of degree 2
	    return nn[1] && nn[2] && is_cycle()
               ? g.size()
               : g.size() - 1;
    else if (deg == 3) {	 // count nodes of degree 3
	    if (     nn[1] == 2 && nn[2] == 2) return 2;
	    else if (nn[1] == 2 || nn[2] == 2) return 1;
	    else if (nn[1] != 2 && nn[2] != 2) return 0;
    }
    else msgr.exit("Error in Mito::num_nodes(). Not implemented deg", deg);
    return huge<szt>;
}


inline
ulong Segment<3>::
set_end_fin( const szt e )
{
    auto& f = g[end2a(e)].fin;
    f[e-1] = nn[e] ? 1 : 0;

    return f[e-1];
}


inline
ulong Segment<3>::
set_bulk_fin( const szt a )
{
    g[a].fin[1] = g[a+1].fin[0] = 1;

    return g[a].fin[1];
}


void Segment<3>::
print( const szt w,
       const std::string& tag, 
       const szt at ) const
{
    if (msgr.so) print(*msgr.so, w, tag, at);
    if (msgr.sl) print(*msgr.sl, w, tag, at);
}


void Segment<3>::
print( std::ostream& os, 
       const szt w,
       const std::string& tag, 
       const szt at ) const
{
    os << "\t" << tag << w;
    if (at == huge<szt>) os << "(of ";
    else os << "(at " << at << " of ";
    os << g.size() << ") [ ";    
    for (szt i=1; i<=nn[1]; i++) os << neig[1][i] << " ";
    os << "] { ";	    
    for (szt i=1; i<=nn[1]; i++) os << neen[1][i] << " ";    
    os << "} [ ";
    for (szt i=1; i<=nn[2]; i++) os << neig[2][i] << " ";
    os << "] { ";    
    for (szt i=1; i<=nn[2]; i++) os << neen[2][i] << " ";    
    os << "} " << cl;    
#ifdef PRINT_EDGES
    os << std::endl;
    for (szt i=0; i<g.size(); i++) 
	    g[i].print(os, i, 1);
#else
    os << " len " << g.size();
#endif
    os << std::endl;
}


void Segment<3>::
write( std::ofstream& ofs ) const
{
    szt len = (szt)g.size();
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
    for (szt j=0; j<len; j++)
	    g[j].write(ofs);
}

}    // namespace MitoSim

#endif // SEGMENT_H
