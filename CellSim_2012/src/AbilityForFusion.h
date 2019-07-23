/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann.
   All Rights Reserved.
   Developed at Helmholtz Center for Infection Research, Braunschweig, Germany.
   Please see Readme file for further information

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

============================================================================== */

#ifndef AbilityForFusion_h
#define AbilityForFusion_h

#include <vector>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "Structure.h"
#include "CoreTransformer.h"
#include "AbilityForFission.h"

namespace MitoD {

/**
 * The AbilityForFusion class.
 * Adds node type-specific fusion capability and updates the network for it.
 * Forms base for clases adding more specific tapes of dynamics.
 */
template<typename Mt>
class AbilityForFusion
	: public AbilityForFission<Mt> {

public:

	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::verbose;
	using Structure<Mt>::oel;
	using CoreTransformer<Mt>::update_cl_fuse;
	using CoreTransformer<Mt>::fuse_antiparallel;
	using CoreTransformer<Mt>::fuse_parallel;
	using AbilityForFission<Mt>::fiss2;

	/** Constructor.
	 * @param cfg configuration object
	 * @param oel logging facility object
	 */
	explicit AbilityForFusion(const Config& cfg, Oel& oel);

	/** Fuse two nodes of degree 1.
	 * @param w1 segment index of the 1st fusion partner
	 * @param e1 segment end of the 1st fusion partner
	 * @param w2 segment index of the 2nd fusion partner
	 * @param e2 segment end of the 2nd fusion partner
	 * @param verb verbasity
	 */
	std::array<szt,2> fuse11(const szt w1, const szt e1, const szt w2, const szt e2, const bool verb) noexcept;

	/** Fuse a node of degree 1 to a node of degree 2.
	 * @param w1 segment index of the fusion partner containing the node of degree 1
	 * @param end segment end of the fusion partner containing the node of degree 1
	 * @param w2 segment index of the fusion partner containing the node of degree 2
	 * @param a2 position of the node of degree 2 relative to the segment starting edge
	 * @param verb verbasity
	 */
	std::array<szt,2> fuse12(const szt w1, const szt end, const szt w2, const szt a2, const bool verb) noexcept;

	/** Fuse a node of degree 1 to the end node in a disconnected loop.
	 * @param w1 segment index of the fusion partner containing the node of degree 1
	 * @param e1 segment end of the fusion partner containing the node of degree 1
	 * @param w2 segment index of the loop segment
	 * @param verb verbasity
	 */
	std::array<szt,2> fuse1L(const szt w1, const szt e1, const szt w2, const bool verb) noexcept;

	/** Fuse end nodes a disconnected segment having free ends to form a loop.
	 * @param w segment index
	 * @param verb verbasity
	 */
	std::array<szt,2> fuse_to_loop(const szt w, const bool verb) noexcept;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFusion<Mt>::
AbilityForFusion(
		const Config& cfgMain,
		Oel& oel
	)
	: AbilityForFission<Mt> {cfgMain, oel}
{}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse11( const szt w1, const szt e1, const szt w2, const szt e2, const bool verbose ) noexcept
{
	if (     w2 == w1)	return fuse_to_loop(w1, verbose);
	else if (e1 == e2)	return fuse_antiparallel(e1, w1, w2, verbose);
	else if (e1 == 1 )  return fuse_parallel(w1, w2, verbose);
	else				return fuse_parallel(w2, w1, verbose);
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse12( const szt w1, const szt end, const szt w2, const szt a2, const bool verbose ) noexcept
{
	if (verbose) {
		oel.print("Fusion12: %d(of %d e %d) with %d(of %d at %d)\n", w1, mt[w1].g.size(), end, w2, mt[w2].g.size(), a2);
		mt[w1].print(w1, "     before s: ");
		mt[w2].print(w2, "     before s: ");
	}
	const auto cl1 {mt[w1].cl};
	const auto cl2 {mt[w2].cl};

	auto mi {mt[w2].isCycle() ? w2 : mtnum+1};	// mt[mi] is to be  produced by the fission

	fiss2(w2, a2, verbose);

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

	if (verbose) {
		mt[w1].print(w1, "       producing ");
		if (w2 != w1)
			mt[w2].print(w2, "                 ");
		if (!mt[w2].isCycle())
			mt[mi].print( mi, "             and " );
		oel.so << std::endl;
		oel.sl << std::endl;
	}
	return {cl1, cl2};
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse1L( const szt w1, const szt e1, const szt w2, const bool verbose ) noexcept
{
	if (verbose) {
		oel.print("Fusion1L: %d(of %d e %d) with a CYCLE %d(of %d)\n", w1, mt[w1].g.size(), e1, w2, mt[w2].g.size());
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

	if (verbose) {
		mt[w1].print(w1, "       producing ");
		mt[w2].print(w2, "             and ");
		oel.so << std::endl;
		oel.sl << std::endl;
	}
	return {cl1, cl2};
}

template<typename Mt>
std::array<szt,2> AbilityForFusion<Mt>::
fuse_to_loop( const szt w, const bool verbose ) noexcept
{
	XASSERT(!mt[w].isCycle(), "Error: attempt to fuse_to_loop a separate cycle.\n");
	XASSERT(!mt[w].nn[1] && !mt[w].nn[2], "Error: attempt to fuse_toLoop a not separate segment.\n");

	if (verbose) {
		oel.print("Fused to cycle %d of length %d", w, mt[w].g.size());
		mt[w].print(w, "Before ", 0);
	}

	mt[w].nn[1] = mt[w].nn[2] = 1;
	mt[w].neig[1][1] = mt[w].neig[2][1] = w;
	mt[w].neen[1][1] = 2; mt[w].neen[2][1] = 1; 

	if (verbose) {
		oel.print("Producing ");
		mt[w].print(w, "After ", 0);
	}
	return {mt[w].cl, mt[w].cl};
}

}
#endif /* AbilityForFusion_h */
