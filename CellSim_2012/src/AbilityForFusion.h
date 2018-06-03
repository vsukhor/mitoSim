#ifndef AbilityForFusion_h
#define AbilityForFusion_h

#include <vector>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "Structure.h"
#include "CoreTransformer.h"
#include "AbilityForFission.h"

namespace MitoD {

template<typename Mt>
class AbilityForFusion : public AbilityForFission<Mt> {

public:

	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::verbose;
	using Structure<Mt>::oel;
	using CoreTransformer<Mt>::update_cl_fuse;
	using CoreTransformer<Mt>::fuse_antiparallel;
	using CoreTransformer<Mt>::fuse_parallel;
	using AbilityForFission<Mt>::fiss2;

	explicit AbilityForFusion(const MitoD::ConfigMain&, Oel&);

	std::array<szt,2> fuse11( const szt, const szt, const szt, const szt, const bool ) noexcept;
	std::array<szt,2> fuse12( const szt w1, const szt end, const szt w2, const szt a2, const bool verbose ) noexcept;
	std::array<szt,2> fuse1L( const szt w1, const szt e1, const szt w2, const bool verbose ) noexcept;
	std::array<szt,2> fuse_to_loop( const szt w, const bool verbose ) noexcept;

private:

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFusion<Mt>::
AbilityForFusion(
		const MitoD::ConfigMain& cfgMain,
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
	const auto cl1 = mt[w1].cl;
	const auto cl2 = mt[w2].cl;

	auto mi = mt[w2].isPureLoop() ? w2 : mtnum+1;	// mt[mi] is to be  produced by the fission

	fiss2(w2, a2, verbose);

	if (w1 == w2) {			// then this is not a loop mito because the loop requires neibs at both ends while w1 is allowed to have a neib at only one end
		if (end == 1) {
			mt[w1].nn[1] = 2;
			mt[w1].neib[1][1] = w1;		mt[w1].neen[1][1] = 2;
			mt[w1].neib[1][2] = mi;		mt[w1].neen[1][2] = 1;

			mt[w1].nn[2] = 2;
			mt[w1].neib[2][1] = w1;		mt[w1].neen[2][1] = 1;
			mt[w1].neib[2][2] = mi;		mt[w1].neen[2][2] = 1;

			mt[mi].nn[1] = 2;
			mt[mi].neib[1][1] = w1;		mt[mi].neen[1][1] = 1;
			mt[mi].neib[1][2] = w1;		mt[mi].neen[1][2] = 2;
		}
		else {
			mt[w1].nn[2] = 2;
			mt[w1].neib[2][1] = mi;		mt[w1].neen[2][1] = 1;
			mt[w1].neib[2][2] = mi;		mt[w1].neen[2][2] = 2;

			mt[mi].nn[1] = 2;
			mt[mi].neib[1][1] = w1;		mt[mi].neen[1][1] = 2;
			mt[mi].neib[1][2] = mi;		mt[mi].neen[1][2] = 2;

			mt[mi].nn[2] = 2;
			mt[mi].neib[2][1] = w1;		mt[mi].neen[2][1] = 2;
			mt[mi].neib[2][2] = mi;		mt[mi].neen[2][2] = 1;
		}
	}
	else {
		mt[w1].nn[end] = 2;
		mt[w1].neib[end][1] = w2;		mt[w1].neen[end][1] = 2;
		mt[w1].neib[end][2] = mi;		mt[w1].neen[end][2] = 1;

		mt[w2].nn[2] = 2;
		mt[w2].neib[2][1] = w1;			mt[w2].neen[2][1] = end;
		mt[w2].neib[2][2] = mi;			mt[w2].neen[2][2] = 1;

		mt[mi].nn[1] = 2;
		mt[mi].neib[1][1] = w1;			mt[mi].neen[1][1] = end;
		mt[mi].neib[1][2] = w2;			mt[mi].neen[1][2] = 2;
	}

	if(mt[w2].cl != mt[mi].cl) update_cl_fuse(mt[w2].cl, mt[mi].cl);
	if(mt[w2].cl != mt[w1].cl) update_cl_fuse(mt[w1].cl, mt[w2].cl);

	if (verbose) {
		mt[w1].print(w1, "       producing ");
		if (w2 != w1)
			mt[w2].print(w2, "                 ");
		if (!mt[w2].isPureLoop())
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
	mt[w1].neib[e1][1] = w2;		mt[w1].neen[e1][1] = 1;
	mt[w1].neib[e1][2] = w2;		mt[w1].neen[e1][2] = 2;

	// update w2 at end 1
	mt[w2].nn[1] = 2;
	mt[w2].neib[1][1] = w2;			mt[w2].neen[1][1] = 2;
	mt[w2].neib[1][2] = w1;			mt[w2].neen[1][2] = e1;

	// update w2 at end 2
	mt[w2].nn[2] = 2;
	mt[w2].neib[2][1] = w2;			mt[w2].neen[2][1] = 1;
	mt[w2].neib[2][2] = w1;			mt[w2].neen[2][2] = e1;

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
	XASSERT(!mt[w].isPureLoop(), "Error: attempt to fuse_to_loop a pure loop.");
	XASSERT(!mt[w].nn[1] && !mt[w].nn[2], "Error: attempt to fuse_toLoop a not separate mito.");

	if (verbose) {
		oel.print("Fused to loop %d of length %d", w, mt[w].g.size());
		mt[w].print(w, "Before ", 0);
	}

	mt[w].nn[1] = mt[w].nn[2] = 1;
	mt[w].neib[1][1] = mt[w].neib[2][1] = w;
	mt[w].neen[1][1] = 2; mt[w].neen[2][1] = 1; 

	if (verbose) {
		oel.print("Producing ");
		mt[w].print(w, "After ", 0);
	}
	return {mt[w].cl, mt[w].cl};
}

}
#endif /* AbilityForFusion_h */
