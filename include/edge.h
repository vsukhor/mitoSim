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

#include "definitions.h"

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

    using FinT = real;

private:

    szt ind {huge<szt>};    ///< Index network-wide: starts from 0.
    szt indcl {huge<szt>};  ///< Index cluster-wide: starts from 0.
    szt cl {huge<szt>};     ///< Current cluster index.

    /// Contribution to fission propensity at each end.
    std::array<FinT,2> fin {};

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


    /**
     * @brief Constructor.
     * @param ifs Input file stream supplying edge attributes.
     */
    explicit Edge(std::ifstream& ifs)
    {
        read(ifs);
    };

    constexpr auto get_ind() const noexcept { return ind; }
    void set_ind(const szt i) noexcept { ind = i; }

    constexpr auto get_indcl() const noexcept { return indcl; }
    void set_indcl(const szt i) noexcept { indcl = i; }

    constexpr auto get_cl() const noexcept { return cl; }
    void set_cl(const szt c) noexcept { cl = c; }

    constexpr auto get_fin(const szt i) const noexcept { return fin[i]; }
    void set_fin(const szt i,
                 const FinT f) noexcept {
        fin[i] = f;
    }

    /// Swap the edge ends.
    void reflect();

    /**
     * @brief Read the edge from a binary file.
     * @param ofs Output file stream.
     */
    void read(std::ifstream& ofs);

    /**
     * @brief Write the edge to a binary file.
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
read( std::ifstream &ifs )
{
    ifs.read(reinterpret_cast<char*>(&ind), sizeof(szt));
    ifs.read(reinterpret_cast<char*>(&indcl), sizeof(szt));
    ifs.read(reinterpret_cast<char*>(&cl), sizeof(szt));
    ifs.read(reinterpret_cast<char*>(&fin[0]), sizeof(FinT));
    ifs.read(reinterpret_cast<char*>(&fin[1]), sizeof(FinT));
}

template<int ContentT>
void Edge<ContentT>::
write( std::ofstream &ofs ) const
{
    ofs.write(reinterpret_cast<const char*>(&ind), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&indcl), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&cl), sizeof(szt));
    ofs.write(reinterpret_cast<const char*>(&fin[0]), sizeof(FinT));
    ofs.write(reinterpret_cast<const char*>(&fin[1]), sizeof(FinT));
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
