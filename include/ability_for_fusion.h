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
* @file ability_for_fusion.h
* @brief Contains class responsible for the graph fusion capacity.
* @author Valerii Sukhorukov
*/

#ifndef ABILITY_FOR_FUSION_H
#define ABILITY_FOR_FUSION_H

#include <vector>

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "structure.h"
#include "core_transformer.h"
#include "ability_for_fission.h"

namespace MitoSim {

/**
 * @brief The AbilityForFusion class.
 * @details Adds node type-specific fusion capability and updates the network for it.
 * Forms base for clases adding more specific tapes of dynamics.
 * @tparam Mt type of the Edge forming the network
 */
template<typename Mt>
class AbilityForFusion
	: public AbilityForFission<Mt> {

public:

	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::msgr;
	using CoreTransformer<Mt>::update_cl_fuse;
	using CoreTransformer<Mt>::fuse_antiparallel;
	using CoreTransformer<Mt>::fuse_parallel;
	using AbilityForFission<Mt>::fiss2;

	/**@brief Constructor.
	 * @param msgr Output message processor.
	 */
	explicit AbilityForFusion(Msgr& msgr);

	/**@brief Fuse two nodes of degree 1.
	 * @param w1 Segment index of the 1st fusion partner.
	 * @param e1 Segment end of the 1st fusion partner.
	 * @param w2 Segment index of the 2nd fusion partner.
	 * @param e2 Segment end of the 2nd fusion partner.
	 */
	std::array<szt,2> fuse11(const szt w1, const szt e1, const szt w2, const szt e2) noexcept;

	/**@brief Fuse a node of degree 1 to a node of degree 2.
	 * @param w1 Segment index of the fusion partner containing the node of degree 1.
	 * @param end Segment end of the fusion partner containing the node of degree 1.
	 * @param w2 Segment index of the fusion partner containing the node of degree 2.
	 * @param a2 Position of the node of degree 2 relative to the segment starting edge.
	 */
	std::array<szt,2> fuse12(const szt w1, const szt end, const szt w2, const szt a2) noexcept;

	/**@brief Fuse a node of degree 1 to the end node in a disconnected loop.
	 * @param w1 Segment index of the fusion partner containing the node of degree 1.
	 * @param e1 Segment end of the fusion partner containing the node of degree 1.
	 * @param w2 Segment index of the loop segment.
	 */
	std::array<szt,2> fuse1L(const szt w1, const szt e1, const szt w2) noexcept;

	/**@brief Fuse end nodes a disconnected segment having free ends to form a loop.
	 * @param w segment index
	 */
	std::array<szt,2> fuse_to_loop(const szt w) noexcept;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFusion<Mt>::
AbilityForFusion(
		Msgr& msgr
	)
	: AbilityForFission<Mt> {msgr}
{}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse11( const szt w1, const szt e1, const szt w2, const szt e2 ) noexcept
{
	if (     w2 == w1)	return fuse_to_loop(w1);
	else if (e1 == e2)	return fuse_antiparallel(e1, w1, w2);
	else if (e1 == 1 )  return fuse_parallel(w1, w2);
	else				return fuse_parallel(w2, w1);
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse12( const szt w1, const szt end, const szt w2, const szt a2 ) noexcept
{
	if constexpr (verbose) {
		msgr.print("Fusion12: %d(of %d e %d) with %d(of %d at %d)\n", w1, mt[w1].g.size(), end, w2, mt[w2].g.size(), a2);
		mt[w1].print(w1, "     before s: ");
		mt[w2].print(w2, "     before s: ");
	}
	const auto cl1 {mt[w1].cl};
	const auto cl2 {mt[w2].cl};

	auto mi {mt[w2].is_cycle() ? w2 : mtnum+1};	// mt[mi] is to be  produced by the fission

	fiss2(w2, a2);

	if (w1 == w2) {
		// then this is not a cycle segment because the cycle requires neighbs at both ends
		// while w1 is allowed to have a neig at only one end
		if (end == 1) {
			mt[w1].nn[1] = 2;
			mt[w1].neig[1][1] = w1;		mt[w1].neen[1][1] = 2;
			mt[w1].neig[1][2] = mi;		mt[w1].neen[1][2] = 1;

			mt[w1].nn[2] = 2;
			mt[w1].neig[2][1] = w1;		mt[w1].neen[2][1] = 1;
			mt[w1].neig[2][2] = mi;		mt[w1].neen[2][2] = 1;

			mt[mi].nn[1] = 2;
			mt[mi].neig[1][1] = w1;		mt[mi].neen[1][1] = 1;
			mt[mi].neig[1][2] = w1;		mt[mi].neen[1][2] = 2;
		}
		else {
			mt[w1].nn[2] = 2;
			mt[w1].neig[2][1] = mi;		mt[w1].neen[2][1] = 1;
			mt[w1].neig[2][2] = mi;		mt[w1].neen[2][2] = 2;

			mt[mi].nn[1] = 2;
			mt[mi].neig[1][1] = w1;		mt[mi].neen[1][1] = 2;
			mt[mi].neig[1][2] = mi;		mt[mi].neen[1][2] = 2;

			mt[mi].nn[2] = 2;
			mt[mi].neig[2][1] = w1;		mt[mi].neen[2][1] = 2;
			mt[mi].neig[2][2] = mi;		mt[mi].neen[2][2] = 1;
		}
	}
	else {
		mt[w1].nn[end] = 2;
		mt[w1].neig[end][1] = w2;		mt[w1].neen[end][1] = 2;
		mt[w1].neig[end][2] = mi;		mt[w1].neen[end][2] = 1;

		mt[w2].nn[2] = 2;
		mt[w2].neig[2][1] = w1;			mt[w2].neen[2][1] = end;
		mt[w2].neig[2][2] = mi;			mt[w2].neen[2][2] = 1;

		mt[mi].nn[1] = 2;
		mt[mi].neig[1][1] = w1;			mt[mi].neen[1][1] = end;
		mt[mi].neig[1][2] = w2;			mt[mi].neen[1][2] = 2;
	}

	if(mt[w2].cl != mt[mi].cl) update_cl_fuse(mt[w2].cl, mt[mi].cl);
	if(mt[w2].cl != mt[w1].cl) update_cl_fuse(mt[w1].cl, mt[w2].cl);

	if constexpr (verbose) {
		mt[w1].print(w1, "       producing ");
		if (w2 != w1)
			mt[w2].print(w2, "                 ");
		if (!mt[w2].is_cycle())
			mt[mi].print( mi, "             and " );
		if(msgr.so) *msgr.so << std::endl;
		if(msgr.sl) *msgr.sl << std::endl;
	}
	return {cl1, cl2};
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse1L( const szt w1, const szt e1, const szt w2 ) noexcept
{
	if constexpr (verbose) {
		msgr.print("Fusion1U: %d(of %d e %d) with a CYCLE %d(of %d)\n", w1, mt[w1].g.size(), e1, w2, mt[w2].g.size());
		mt[w1].print( w1, "     before s: " );
		mt[w2].print( w2, "     before s: " );
	}
	const auto cl1 = mt[w1].cl;
	const auto cl2 = mt[w2].cl;

	// update w1 at e1
	mt[w1].nn[e1] = 2;
	mt[w1].neig[e1][1] = w2;	mt[w1].neen[e1][1] = 1;
	mt[w1].neig[e1][2] = w2;	mt[w1].neen[e1][2] = 2;

	// update w2 at end 1
	mt[w2].nn[1] = 2;
	mt[w2].neig[1][1] = w2;		mt[w2].neen[1][1] = 2;
	mt[w2].neig[1][2] = w1;		mt[w2].neen[1][2] = e1;

	// update w2 at end 2
	mt[w2].nn[2] = 2;
	mt[w2].neig[2][1] = w2;		mt[w2].neen[2][1] = 1;
	mt[w2].neig[2][2] = w1;		mt[w2].neen[2][2] = e1;

	if (mt[w1].cl != mt[w2].cl )
		update_cl_fuse(mt[w1].cl, mt[w2].cl);

	if constexpr (verbose) {
		mt[w1].print(w1, "       producing ");
		mt[w2].print(w2, "             and ");
		if (msgr.so) *msgr.so << std::endl;
		if (msgr.sl) *msgr.sl << std::endl;
	}
	return {cl1, cl2};
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse_to_loop( const szt w ) noexcept
{
	XASSERT(!mt[w].is_cycle(), "Error: attempt to fuse_to_loop a separate cycle.\n");
	XASSERT(!mt[w].nn[1] && !mt[w].nn[2], "Error: attempt to fuse_toLoop a not separate segment.\n");

	if constexpr (verbose) {
		msgr.print("Fused to cycle %d of length %d", w, mt[w].g.size());
		mt[w].print(w, "Before ", 0);
	}

	mt[w].nn[1] = mt[w].nn[2] = 1;
	mt[w].neig[1][1] = mt[w].neig[2][1] = w;
	mt[w].neen[1][1] = 2; mt[w].neen[2][1] = 1; 

	if constexpr (verbose) {
		msgr.print("Producing ");
		mt[w].print(w, "After ", 0);
	}
	return {mt[w].cl, mt[w].cl};
}

}	// namespace MitoSim

#endif // ABILITY_FOR_FUSION_H
