#ifndef CoreTransformer_h
#define CoreTransformer_h

#include <vector>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "Structure.h"

namespace MitoD {

template<typename Mt>
class CoreTransformer
	: public Structure<Mt> {

public:


	explicit CoreTransformer( const MitoD::ConfigMain&, Oel& );

	using Structure<Mt>::oel;
	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::clnum;

	constexpr void update_mtcl_fuse( const szt w1, const szt w2 ) noexcept;
	constexpr void update_cl_fuse( const szt c1, const szt c2 ) noexcept;
	constexpr void changeCl( const szt cf, const szt ct ) noexcept;
	constexpr void update_gIndcl( const szt cl ) noexcept;

protected:

	std::array<szt,2> fuse_antiparallel( const szt, const szt, const szt, const bool ) noexcept;
	std::array<szt,2> fuse_parallel( const szt, const szt, const bool ) noexcept;

	void rename_mito(const szt, const szt);
	void copy_neibs( const szt f, const szt ef, const szt t, const szt et ) noexcept;
	void update_neibs(	const szt oldn, const szt oend,
						const szt n1, const szt n2,
						const szt newn, const szt nend,
						const bool removefromneibs ) noexcept;

private:

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
CoreTransformer<Mt>::
CoreTransformer(
		const MitoD::ConfigMain& cfgMain,
		Oel& oel
	)
	: Structure<Mt> {cfgMain, oel}
{}

template<typename Mt>
void CoreTransformer<Mt>::
rename_mito( const szt f, const szt t )
{
	copy_neibs(f, 1, t, 1);
	copy_neibs(f, 2, t, 2);
	mt[t].g = mt[f].g;
	mt[t].cl = mt[f].cl;
}

template<typename Mt>
void CoreTransformer<Mt>::
copy_neibs( const szt f, const szt ef, const szt t, const szt et ) noexcept
{
	for (szt j=1; j<=mt[f].nn[ef]; j++) {
		mt[t].neib[et][j] = mt[f].neib[ef][j];
		mt[t].neen[et][j] = mt[f].neen[ef][j];
	}
	mt[t].nn[et] = mt[f].nn[ef];

	update_neibs(f, ef, 1, mt[f].nn[ef], t, et, false);					// substitute f in f's neib's neibs for t
}

template<typename Mt>
void CoreTransformer<Mt>::
update_neibs( const szt oldn, const szt oend,
			  const szt n1,	  const szt n2,
			  const szt newn, const szt nend,
			  const bool removefromneibs ) noexcept
{
	for (szt j=n1; j<=n2; j++) {
		const auto cn = mt[oldn].neib[oend][j];										// our neib currently processed
		const auto ce = mt[oldn].neen[oend][j];										// our neib currently processed

		szt i1 {0};
		while (1) {
			XASSERT(++i1 <= mt[cn].nn[ce], "A neib was not found: "+STR(oldn)+" "+STR(cn));
			if (mt[cn].neib[ce][i1] == oldn &&
				mt[cn].neen[ce][i1] == oend)
				break;
		}
		if (removefromneibs) {
			mt[cn].neib[ce][i1] = mt[cn].neib[ce][mt[cn].nn[ce]];				// remove oldn from the neib list of its j-th neib
			mt[cn].neen[ce][i1] = mt[cn].neen[ce][mt[cn].nn[ce]--];

			mt[oldn].neib[oend][j] = mt[oldn].neib[oend][mt[oldn].nn[oend]];	// remove the j-th neib from the oldn's list of neibs
			mt[oldn].neen[oend][j] = mt[oldn].neen[oend][mt[oldn].nn[oend]--];
		}
		else {
			mt[cn].neib[ce][i1] = newn;
			mt[cn].neen[ce][i1] = nend;
		}
	}
}
template<typename Mt>
std::array<szt,2> CoreTransformer<Mt>::
fuse_antiparallel( const szt end, const szt w1, const szt w2, const bool verbose ) noexcept
{
	const auto len1 = mt[w1].g.size();
	const auto len2 = mt[w2].g.size();
	const auto cl1 = mt[w1].cl;
	const auto cl2 = mt[w2].cl;
	if (verbose) { 
		this->oel.so <<"Fusion11a: "<< w1 << "(of "<< len1 << ") with " << w2 << "(of "<< len2 <<") at end "<< end << std::endl;
		this->oel.sl <<"Fusion11a: "<< w1 << "(of "<< len1 << ") with " << w2 << "(of "<< len2 <<") at end "<< end << std::endl;
		mt[w1].print(w1, "     before a: ");
		mt[w2].print(w2, "     before a: ");
	}
	XASSERT(w1 != w2, "Error during antiparallel fusion: w1 == w2: fuse_toLoop should be used instead.");
	XASSERT(!mt[w1].nn[end], "Error during antiparallel fusion: end of w1 is not free.");
	XASSERT(!mt[w2].nn[end], "Error during antiparallel fusion: end of w2 is not free.");

	const szt opend = (end==2) ? 1 : 2;
	if (end == 1)
		copy_neibs(w1, 2, w1, 1);					// copy w1's 1-end neibs to its 0-end
	copy_neibs(w2, opend, w1, 2);					// copy w2's 1-end neibs to w1's 1-end

	if (mt[w2].cl != mt[w1].cl)
		update_mtcl_fuse(w1, w2);

	if (end == 1) mt[w1].reflect_g();				//  Edge::a takes values [1:g.size()]	// for w1 particles reflect positions if 1-ends are joined;
	else 		  mt[w2].reflect_g();				//  Edge::a takes values [1:g.size()]	// for w2 particles reflect positions if 2-ends are joined;

	mt[w1].g.resize(len1+len2);
	for (szt ia=0; ia<len2; ia++)
		mt[w1].g[len1+ia] = mt[w2].g[ia];

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
	const auto len1 = mt[w1].g.size();
	const auto len2 = mt[w2].g.size();
	const auto cl1 = mt[w1].cl;
	const auto cl2 = mt[w2].cl;

	if(verbose) {
		this->oel.so <<"Fusion11p: "<< w1 << "(of "<< len1 << ") with " << w2 << "(of "<< len2 << ")" << std::endl;
		this->oel.sl <<"Fusion11p: "<< w1 << "(of "<< len1 << ") with " << w2 << "(of "<< len2 << ")" << std::endl;
		mt[w1].print(w1, "     before p: ");
		mt[w2].print(w2, "     before p: ");
	}
	XASSERT(w1 != w2, "Error during parallel fusion: w1 == w2: fuse_toLoop should be used instead.");
	XASSERT(!mt[w1].nn[1], "Error during parallel fusion: end 1 of w1 is not free.");
	XASSERT(!mt[w2].nn[2], "Error during parallel fusion: end 2 of w2 is not free.");

	copy_neibs(w2, 1, w1, 1);
	if(mt[w2].cl != mt[w1].cl)
		update_mtcl_fuse(w1, w2);

	mt[w1].g.resize( len1 + len2 );

	for (szt ia=len1+len2-1; ia>=len2; ia--) mt[w1].g[ia] = mt[w1].g[ia-len2];
	for (szt ia=0; ia<len2; ia++)			 mt[w1].g[ia] = mt[w2].g[ia];

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
	changeCl(c2, c1);
	if (c2 != clnum - 1)
		changeCl(clnum-1, c2);
	clnum--;
}

template<typename Mt> constexpr
void CoreTransformer<Mt>::
changeCl( const szt cf, const szt ct ) noexcept			// args by value
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
	szt indcl {0};
	for (szt j=1; j<=mtnum; j++)
		if (mt[j].cl == cl)
			indcl = mt[j].set_gCl(cl, indcl);
}

}

#endif /* CoreTransformer_h */
