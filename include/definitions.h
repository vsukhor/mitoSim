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
 * @file definitions.h
 * @brief mitosim namespace-scope gefinitions.
 * @author Valerii Sukhorukov.
 */

#ifndef MITOSIM_DEFINITIONS_H
#define MITOSIM_DEFINITIONS_H

#include <filesystem>
#include <string>

#include "utils/constants.h"
#include "utils/msgr.h"
#include "utils/random/with_boost.h"

// #define USE_UTILS_XASSERT  // toggles XASSERTs.

namespace mitosim {

using namespace utils::constants;

using Msgr = utils::Msgr;


using real = float;
using RandFactory = utils::random::Boost<real>;

constexpr bool verbose {};   ///< Work in verbose mode.
constexpr bool print_edges {};   ///< Print detailed data on edges.

}  // namespace mitosim

#endif  // MITOSIM_DEFINITIONS_H
