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
* @file edge.h
* @brief The graph Edge class.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_EDGE_H
#define MITOSIM_EDGE_H

#include <fstream>
#include <ostream>

#include "utils/common/constants.h"

namespace mitosim {

template <unsigned> class Segment;

/**
 * @brief The Network Edge class.
 * @details Edge is a minimal structural unit of the network.
 * The class handles the tasks and properties specific to a single edge
 * and its relation to other network components.
 * @tparam ContentT Slot for specifying strucutre of the internal content
 * the Edge can hold; currently not used.
 */
template<int ContentT>
class Edge {

    friend Segment<3>;

public:

    using szt = utils::common::szt;
    using ulong = utils::common::ulong;

    static constexpr auto hugeszt = utils::common::huge<szt>;

private:

    szt ind {hugeszt};    ///< Index network-wide: starts from 0.
    szt indcl {hugeszt};  ///< Index cluster-wide: starts from 0.
    szt cl {hugeszt};     ///< Current cluster index.

    /// Contribution to fission propensity at each end.
    std::array<ulong,2> fin {{}};

public:

    /**
     * @brief Constructor.
     * @param ind Index network-wide.
     * @param indcl Index cluster-wide.
     * @param cl Current cluster index.
     */
    explicit Edge(
        szt ind,
        szt indcl,
        szt cl
    );

    constexpr auto get_ind() const noexcept { return ind; }
    constexpr auto get_indcl() const noexcept { return indcl; }
    constexpr auto get_cl() const noexcept { return cl; }
    constexpr auto get_fin(const int i) const noexcept { return fin[i]; }
    void set_fin(const int i,
                 const ulong f) noexcept {
        fin[i] = f;
    }

    /// Swap the edge ends.
    void reflect();

    /**
     * @brief Write the edge to a file.
     * @param ofs Output file stream.
     */
    void write(std::ofstream& ofs) const;

    /**
     * @brief Print the edge to a stream.
     * @tparam ENDL Flag to end current line.
     * @param os Output stream.
     * @param a Position inside segment.
    */
    template<bool ENDL=true>
    void print(std::ostream& os, szt a) const;
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int ContentT>
Edge<ContentT>::
Edge(
    const szt ind,
    const szt indcl,
    const szt cl
    )
    : ind {ind}
    , indcl {indcl}
    , cl {cl}
{}


template<int ContentT>
void Edge<ContentT>::
reflect()
{
    std::swap(fin[0], fin[1]);
}


template<int ContentT>
void Edge<ContentT>::
write( std::ofstream &ofs ) const
{
    ofs.write(reinterpret_cast<const char*>(&ind), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&indcl), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&cl), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&fin[0]), sizeof(ulong));
    ofs.write(reinterpret_cast<const char*>(&fin[1]), sizeof(ulong));
}


template<int ContentT>
template<bool ENDL>
void Edge<ContentT>::
print( std::ostream& os,
       const szt a ) const
{
    os << "[" << a << "] ";
    os << " ind " << ind; 
    os << " indcl " << indcl;
    os << " fin " << fin[0] << " " << fin[1];
    if constexpr (ENDL) os << "\n";
}

}  // namespace mitosim

#endif  // MITOSIM_EDGE_H
