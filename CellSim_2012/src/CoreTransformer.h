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

#ifndef CoreTransformer_h
#define CoreTransformer_h

#include <vector>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "Structure.h"

namespace MitoD {

/**
 * The CoreTransformer class.
 * Encapsulates a elementary dynamics and updates the network for it.
 * Forms base for clases adding more specific tapes of dynamics.
 */
template<typename Mt>
class CoreTransformer
	: public Structure<Mt> {

public:

	using Structure<Mt>::oel;
	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::clnum;

	/** Constructor,
	 * @param cfg configuration object
	 * @param oel logging facility object
	 */
	explicit CoreTransformer(const Config& cfg, Oel& oel);

	/** Updates the structure of diaconnected components resulting from a fusion event.
	 * @param w1 index of the 1st participant segment
	 * @param w2 index of the 2nd participant segment
	 */
	constexpr void update_mtcl_fuse(const szt w1, const szt w2) noexcept;

	/** Updates the structure of diaconnected components resulting from a fusion event.
	 * @param c1 index of the 1st participant component
	 * @param c2 index of the 2nd participant component
	 */
	constexpr void update_cl_fuse(const szt c1, const szt c2) noexcept;

	/** Updates the diaconnected component indexes.
	 * @param cf initial index
	 * @param ct final index
	 */
	constexpr void update_cl(const szt cf, const szt ct) noexcept;

	/** Updates the diaconnected component indexes.
	 * @param c initial index
	 */
	constexpr void update_gIndcl(const szt c) noexcept;

protected:

	/** Performs fusion of two antiparallel oriented segments.
	 * Antiparallel orientation is defined as
	 * fusion of end 1 of the 1st segment to to end 1 of the 2nd segment, or
	 * fusion of end 2 of the 1st segment to to end 2 of the 2nd segment.
	 * @param end end index at the fusion position
	 * @param w1 index of the 1st segment
	 * @param w2 index of the 2nd segment
	 * @param verb verbasity
	 */
	std::array<szt,2> fuse_antiparallel(const szt end, const szt w1, const szt w2, const bool verb) noexcept;

	/** Performs fusion of two parallel oriented segments.
	 * Parallel orientation is defined as
	 * fusion of end 1 of the 1st segment to to end 2 of the 2nd segment, or
	 * fusion of end 2 of the 1st segment to to end 1 of the 2nd segment.
	 * @param w1 index of the 1st segment
	 * @param w2 index of the 2nd segment
	 * @param verb verbasity
	 */
	 std::array<szt,2> fuse_parallel(const szt w1, const szt w2, const bool verb) noexcept;

	/** Update the network segment indexes such that segment indexed as source will become target.
	 * @param f source segment index
	 * @param t target segment index
	*/
	void rename_mito(const szt f, const szt t);

protected:

	/** Copy connection partners to a new segment.
	 * @param f source segment index
	 * @param ef source segment end
	 * @param t target segment index
	 * @param et target segment end
	*/
	void copy_neigs(const szt f, const szt ef, const szt t, const szt et) noexcept;

	/** Update segment neighbours.
	 * @param oldn old neighbour segment index
	 * @param oend old neighbour segment end
	 * @param n1 neighbour initial index
	 * @param n2 neighbour final index
	 * @param newn new neighbour segment index
	 * @param nend new neighbour segment end
	 * @param removefromneigs flag if the neighbouring segments should be removed
	*/
	void update_neigs(	const szt oldn, const szt oend,
						const szt n1, const szt n2,
						const szt newn, const szt nend,
						const bool removefromneigs ) noexcept;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
CoreTransformer<Mt>::
CoreTransformer(
		const MitoD::Config& cfg,
		Oel& oel
	)
	: Structure<Mt> {cfg, oel}
{}

template<typename Mt>
void CoreTransformer<Mt>::
rename_mito( const szt f, const szt t )
{
	copy_neigs(f, 1, t, 1);
	copy_neigs(f, 2, t, 2);
	mt[t].g = std::move(mt[f].g);
	mt[t].cl = mt[f].cl;
}

template<typename Mt>
void CoreTransformer<Mt>::
copy_neigs( const szt f, const szt ef, const szt t, const szt et ) noexcept
{
	for (szt j=1; j<=mt[f].nn[ef]; j++) {
		mt[t].neig[et][j] = mt[f].neig[ef][j];
		mt[t].neen[et][j] = mt[f].neen[ef][j];
	}
	mt[t].nn[et] = mt[f].nn[ef];

	update_neigs(f, ef, 1, mt[f].nn[ef], t, et, false);		// substitute f in f's neig's neigs for t
}

template<typename Mt>
void CoreTransformer<Mt>::
update_neigs( const szt oldn, const szt oend,
			  const szt n1,	  const szt n2,
			  const szt newn, const szt nend,
			  const bool removefromneigs ) noexcept
{
	for (szt j=n1; j<=n2; j++) {
		const auto cn = mt[oldn].neig[oend][j];			// our neig currently processed
		const auto ce = mt[oldn].neen[oend][j];			// our neig currently processed

		szt i1 {};
		while (1) {
			XASSERT(++i1 <= mt[cn].nn[ce], "A neig was not found: "+STR(oldn)+" "+STR(cn)+"\n");
			if (mt[cn].neig[ce][i1] == oldn &&
				mt[cn].neen[ce][i1] == oend)
				break;
		}
		if (removefromneigs) {
			// remove oldn from the neig list of its j-th neig
			mt[cn].neig[ce][i1] = mt[cn].neig[ce][mt[cn].nn[ce]];
			mt[cn].neen[ce][i1] = mt[cn].neen[ce][mt[cn].nn[ce]--];
			// remove the j-th neig from the oldn's list of neigs
			mt[oldn].neig[oend][j] = mt[oldn].neig[oend][mt[oldn].nn[oend]];
			mt[oldn].neen[oend][j] = mt[oldn].neen[oend][mt[oldn].nn[oend]--];
		}
		else {
			mt[cn].neig[ce][i1] = newn;
			mt[cn].neen[ce][i1] = nend;
		}
	}
}
template<typename Mt>
std::array<szt,2> CoreTransformer<Mt>::
fuse_antiparallel( const szt end, const szt w1, const szt w2, const bool verbose ) noexcept
{
	const auto len1 {mt[w1].g.size()};
	const auto len2 {mt[w2].g.size()};
	const auto cl1 {mt[w1].cl};
	const auto cl2 {mt[w2].cl};

	if (verbose) {
		oel.print("Fusion11a: %d(of %d) with %d(of %d) at end %d", w1, len1, w2, len2, end);
		mt[w1].print(w1, "     before a: ");
		mt[w2].print(w2, "     before a: ");
	}
	XASSERT(w1 != w2, "Error during antiparallel fusion: w1 == w2: fuse_toLoop should be used instead.\n");
	XASSERT(!mt[w1].nn[end], "Error during antiparallel fusion: end of w1 is not free.\n");
	XASSERT(!mt[w2].nn[end], "Error during antiparallel fusion: end of w2 is not free.\n");

	const szt opend = (end==2) ? 1 : 2;
	if (end == 1)
		copy_neigs(w1, 2, w1, 1);		// copy w1's 1-end neigs to its 0-end
	copy_neigs(w2, opend, w1, 2);		// copy w2's 1-end neigs to w1's 1-end

	if (mt[w2].cl != mt[w1].cl)
		update_mtcl_fuse(w1, w2);

	if (end == 1) mt[w1].reflect_g();	//  Edge::a takes values [1:g.size()]; for w1 reflect positions if 1-ends are joined;
	else 		  mt[w2].reflect_g();	//  Edge::a takes values [1:g.size()]; for w2 reflect positions if 2-ends are joined;

	std::move(mt[w2].g.begin(), mt[w2].g.end(), std::back_inserter(mt[w1].g));
	mt[w2].g.clear();

	if (w2 != mtnum)
		rename_mito(mtnum, w2);
	mt.pop_back();
	mtnum--;

	update_gIndcl(cl1);
	if (cl1 != cl2)
		update_gIndcl(cl2);

	if (verbose) {
		if (w1 == mtnum+1)	{ mt[w2].print( w2, "       producing " ); oel.so << std::endl; oel.sl << std::endl; }
		else 				{ mt[w1].print( w1, "       producing " ); oel.so << std::endl; oel.sl << std::endl; }
	}
	return {cl1, cl2};
}

template<typename Mt>
std::array<szt,2> CoreTransformer<Mt>::
fuse_parallel( const szt w1, const szt w2, const bool verbose ) noexcept
{
	const auto len1 {mt[w1].g.size()};
	const auto len2 {mt[w2].g.size()};
	const auto cl1 {mt[w1].cl};
	const auto cl2 {mt[w2].cl};

	if(verbose) {
		oel.print("Fusion11p: %d(of %d) with %d(of %d)", w1, len1, w2, len2);
		mt[w1].print(w1, "     before p: ");
		mt[w2].print(w2, "     before p: ");
	}
	XASSERT(w1 != w2, "Error during parallel fusion: w1 == w2: fuse_toLoop should be used instead.\n");
	XASSERT(!mt[w1].nn[1], "Error during parallel fusion: end 1 of w1 is not free.\n");
	XASSERT(!mt[w2].nn[2], "Error during parallel fusion: end 2 of w2 is not free.\n");

	copy_neigs(w2, 1, w1, 1);
	if(mt[w2].cl != mt[w1].cl)
		update_mtcl_fuse(w1, w2);

	std::move(mt[w1].g.begin(), mt[w1].g.end(), std::back_inserter(mt[w2].g));
	mt[w1].g = std::move(mt[w2].g);

	if (w2 != mtnum)
		rename_mito(mtnum, w2);
	mt.pop_back();
	mtnum--;

	update_gIndcl(cl1);
	if (cl1 != cl2)
		update_gIndcl(cl2);

	if (verbose) {
		if (w1 == mtnum+1)	{ mt[w2].print( w2, "       producing " ); oel.so << std::endl; oel.sl << std::endl; }
		else 				{ mt[w1].print( w1, "       producing " ); oel.so << std::endl; oel.sl << std::endl; }
	}

	return {cl1, cl2};
}

template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_mtcl_fuse( const szt w1, const szt w2 ) noexcept
{
	const auto w1cl = mt[w1].cl;
	const auto w2cl = mt[w2].cl;

	for (szt i=1; i<=mtnum; i++)
		if (mt[i].cl == w2cl)
			mt[i].cl = w1cl;

	if (w2cl != clnum-1)
		for(szt i=1; i<=mtnum; i++)
			if(mt[i].cl == clnum-1)
				mt[i].cl = w2cl;
	clnum--;
}

template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_cl_fuse( const szt c1, const szt c2 ) noexcept	// args by value
{
	update_cl(c2, c1);
	if (c2 != clnum - 1)
		update_cl(clnum-1, c2);
	clnum--;
}

template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_cl( const szt cf, const szt ct ) noexcept			// args by value
{
	for (szt i=1; i<=mtnum; i++)
		if (mt[i].cl == cf)
			mt[i].cl = ct;

	update_gIndcl(ct);
}

template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_gIndcl( const szt cl ) noexcept
{
	szt indcl {};
	for (szt j=1; j<=mtnum; j++)
		if (mt[j].cl == cl)
			indcl = mt[j].set_gCl(cl, indcl);
}

}

#endif /* CoreTransformer_h */
