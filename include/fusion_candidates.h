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
 * @file fusion_candidates.h
 * @brief Contains classes for storing nodes available for fusion.
 * @author Valerii Sukhorukov
 */

#ifndef MITOSIM_FUSION_CANDIDATES_H
#define MITOSIM_FUSION_CANDIDATES_H

namespace mitosim {

/**
 * @brief Container for fusion candidate nodes.
 * @note Should be used only for the reactions not involving fusion to a
 * disconnected cycle segment.
 */
struct alignas(8) FusionCandidatesXX {

    static constexpr int MIN_ALIGNMENT = 8;

    /// Segment and end indexes of the 1st participant.
    alignas(MIN_ALIGNMENT) std::vector<std::array<szt,2>> u;
    ///< Segment and end indexes of the 2nd participant.
    alignas(MIN_ALIGNMENT) std::vector<std::array<szt,2>> v;

    /// Empty the container.
    void clear() noexcept
    {
        u.clear();
        v.clear();
    }

    /**
     * @brief Add a node pair.
     * @param uc Segment and end indexes of the 1st participant.
     * @param vc Segment and end indexes of the 2nd participant.
     */
    void add(const std::array<szt,2>& uc,
             const std::array<szt,2>& vc)
    {
        u.emplace_back(uc);
        v.emplace_back(vc);
    }

    /**
     * @brief Report the number of elements.
     * @result Current number of candidate fusion pairs.
     */
    szt size() const noexcept { return u.size(); }


};

////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Container for fusion candidate nodes.
 * @note Should be used only for fusion reactions between a cycle and a
 * non-cycle segments.
 */
struct alignas(8) FusionCandidatesXU {

    static constexpr int MIN_ALIGNMENT = 8;

    /// Segment and end indexes of the non-looped participant.
    alignas(MIN_ALIGNMENT) std::vector<std::array<szt,2>> u;
    /// Segment indexes of the cycle participant.
    alignas(MIN_ALIGNMENT) std::vector<szt> v;

    /// Empty the container.
    void clear() noexcept
    {
        u.clear();
        v.clear();
    }

    /**
     * @brief Add a node pair.
     * @param uc Segment and end indexes of the non-cycle participant.
     * @param vc Segment indexes of the cycle participant.
    */
    void add(std::array<szt,2> uc, szt vc)
    {
        u.emplace_back(uc);
        v.emplace_back(vc);
    }

    /**
     * @brief Report the number of elements.
     * @return Current number of elements.
     */
    szt size() const noexcept { return u.size(); }
};

}  // namespace mitosim

#endif  // MITOSIM_FUSION_CANDIDATES_H
