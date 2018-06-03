#ifndef AbilityForFission_h
#define AbilityForFission_h

#include <vector>
#include <algorithm>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "CoreTransformer.h"

namespace MitoD {

template<typename> class Fission;

template<typename Mt>
class AbilityForFission : public CoreTransformer<Mt> {

public:

	typedef AbilityForFission<Mt> thisT;
	typedef Fission<thisT> ThisFission;

	using Structure<Mt>::oel;
	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::clnum;
	using Structure<Mt>::glm;
	using Structure<Mt>::verbose;
	using Structure<Mt>::basic_update;
	using Structure<Mt>::update_structure_ind;
	using CoreTransformer<Mt>::copy_neibs;
	using CoreTransformer<Mt>::update_neibs;
	using CoreTransformer<Mt>::fuse_antiparallel;
	using CoreTransformer<Mt>::fuse_parallel;
	using CoreTransformer<Mt>::update_gIndcl;

	explicit AbilityForFission( const MitoD::ConfigMain&, Oel& );

	std::array<szt,2> fiss( const szt, const szt, const bool);
	std::array<szt,2> fiss2(const szt, const szt, const bool);
	std::array<szt,2> fiss3(const szt, const szt, const bool);

private:

	std::vector<szt>  vis;	// auxiliary variable indicating a visited or not status of mitos during the graph search

	bool update_cl_fiss(const szt, const szt);
	bool dfs(const szt, const szt, const szt, const szt);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFission<Mt>::
AbilityForFission(
		const MitoD::ConfigMain& cfgMain,
		Oel& oel
	)
	: CoreTransformer<Mt> {cfgMain, oel}

{}

template<typename Mt>
bool AbilityForFission<Mt>::
update_cl_fiss( const szt w, const szt e )
{
	vis.resize(mtnum+2);
	std::fill(vis.begin()+1, vis.begin()+mtnum+1, 0);			// nothing is visited in the beginning of the search

	const szt oe {(e == 1) ? szt(2) : szt(1)};

	const bool loop = dfs(w, e, w, oe);
	if (!loop) {
		clnum++;
		szt clind {0};
		for (szt i=1; i<=mtnum; i++)
			if (vis[i])
				clind = mt[i].setCl(clnum-1, clind);
	}
	return loop;
}

// depth-first search of the network graph
template<typename Mt>
bool AbilityForFission<Mt>::
dfs( const szt w1, const szt e1,
	 const szt w2, const szt e2 )
{
	for (szt i=1; i<=mt[w1].nn[e1]; i++) {
		const auto cn = mt[w1].neib[e1][i];
		const auto ce = mt[w1].neen[e1][i]; 
		if (cn == w2) {
			if(ce == e2)
				return true;
		}
		else {
			if(!vis[cn]) {
				vis[cn] = 1;
				const szt ne {(ce == 1) ? szt(2) : szt(1)};
				if (dfs(cn, ne, w2, e2))
					return true;
			}
		}	
	}
	return false;
}

// a is counted from 1 and is the last element to remain in the old mito
template<typename Mt>
std::array<szt,2> AbilityForFission<Mt>::
fiss( const szt w, const szt a, const bool verbose )
{
	if (        a && a < mt[w].g.size())											 return fiss2(w, a, verbose);	// node 2 -> nodes 1+1; cuts between g[a-1] and g[a]
	else if (( !a && mt[w].nn[1] <= 2) || (a == mt[w].g.size() && mt[w].nn[2] <= 2)) return fiss3(w, a, verbose);	// node 3 -> nodes 2+1 and node 2 -> nodes 1+1 in pure loops
	else oel.exit("ERROR: Attempt of an unpropriate fission!");
	return {huge<szt>, huge<szt>};
}

// divides at a node of degree 2
// 'w' is a global segment index
// 'a' is the node position inside the segment 'w'
// 'a' is counted from 1 and is the last element to remain in the original segment
template<typename Mt>
std::array<szt,2> AbilityForFission<Mt>::
fiss2( const szt w, const szt a, const bool verbose )
{
	if(verbose)
		mt[w].print(w, "fission2:  ", a);					// cuts between g[a-1] and g[a]

	szt ind1 = mt[w].g[a-1].ind;
	szt ind2 = mt[w].g[a].ind;

	bool clusterisaloop = 0;
	mt[w].nn[2] ? (clusterisaloop = update_cl_fiss(w, 2))
				: clnum++;

	mt.emplace_back(oel);
	++mtnum;

	mt[mtnum].g.resize(mt[w].g.size()-a);

	for (szt ia=a; ia<mt[w].g.size(); ia++)
		mt[mtnum].g[ia-a] = mt[w].g[ia];

	mt[w].g.resize(a);
	mt[mtnum].nn[1] = 0;

	this->copy_neibs(w, 2, mtnum, 2);

	mt[mtnum].cl = clusterisaloop
				 ? mt[w].cl
				 : clnum - 1;

	if (!clusterisaloop) {
		update_gIndcl(mt[w].cl);		// renumber Edge::indcl of the remaining part of the original cluster
		update_gIndcl(clnum-1);			// renumber Edge::indcl of the newly formed cluster
	}

	mt[w].nn[2] = 0;

	bool isSelfLooped {false};
	if (1 == mt[w].nn[1] &&
		1 == mt[mtnum].nn[2] && 
		w == mt[mtnum].neib[2][mt[mtnum].singleNeibInd(2)] && 
		mtnum == mt[w].neib[1][mt[w].singleNeibInd(1)] ) {			// before fission 'w' was looped into itself
		isSelfLooped = true;
		if (verbose) {
			mt[w].print(w, "       transiently producing ", -1);
			mt[mtnum].print(mtnum, "             and ", -1);
			std::cout << std::endl;
		}
		update_neibs(w, 1, 1, 1, -1, -1, true);
		fuse_parallel(w, mtnum, verbose);
	}

	basic_update();
	const auto w1 = glm[ind1];
	const auto w2 = glm[ind2];

	if(verbose) {
		mt[w1].print(w1, "       producing ");
		if (isSelfLooped)
			oel.print("from a mito looped into itself");
		else mt[w2].print(w2, "             and ");
		std::cout << std::endl;
	}

	return {{mt[w].cl, mt[mtnum].cl}};
}

// divides at a node of degree 3
// 'w' is a global segment index
// 'at' is the node position inside the segment 'w'
// 'at' is counted from 1 and is the last element to remain in the original segment
template<typename Mt>
std::array<szt,2> AbilityForFission<Mt>::
fiss3( const szt w, const szt a, const bool verbose )
{
	if(verbose)
		mt[w].print(w, "fission3:  ", a);

	const auto clini {mt[w].cl};
	bool clusterisaloop {false};
	bool f {false};
	szt n[2], e[2];
	auto ind1 {huge<szt>};
	auto ind2 {huge<szt>};
	if (!a) {															// at end 1
		ind1 = mt[w].g.front().ind;
		ind2 = mt[mt[w].neib[1][1]].g[mt[mt[w].neib[1][1]].end2a(mt[w].neen[1][1])].ind;
		if (mt[w].nn[1] == 2) {
			const auto ninds = mt[w].doubleNeibInds(1);
			f = true;
			for (int j=0; j<2; j++) {
				n[j] = mt[w].neib[1][ninds[j]];
				e[j] = mt[w].neen[1][ninds[j]];
			}
		}
		else if (mt[w].nn[1] == 1)
			n[0] = mt[w].neib[1][1];

		clusterisaloop = update_cl_fiss(w, 1);		// if not a loop, this increments clnum and forms a new cluster from w's end 1 neibs and beyond (excluding w itself)
		if(!clusterisaloop)
			update_gIndcl(clini);					// renumber Edge::indcl of the remaining part of the original cluster

		update_neibs(w, 1, 1, mt[w].nn[1], -1, -1, true);

		if (f && n[0] != n[1]) {
			if (mt[n[0]].nn[e[0]] == 1 && 
			    mt[n[1]].nn[e[1]] == 1 ) {

				const auto indneib0 = mt[n[0]].singleNeibInd(e[0]);
				const auto indneib1 = mt[n[1]].singleNeibInd(e[1]);

				if (mt[n[0]].neib[e[0]][indneib0] == n[1] && 
				    mt[n[0]].neen[e[0]][indneib0] == e[1] && 
				    mt[n[1]].neib[e[1]][indneib1] == n[0] && 
				    mt[n[1]].neen[e[1]][indneib1] == e[0] ) {

					update_neibs( n[0], e[0], 1, 1, -1, -1, true );

					if (e[0] == e[1])
						fuse_antiparallel(e[0], n[0], n[1], verbose);
					else
						(e[0] == 1 && e[1] == 2) ? fuse_parallel(n[0], n[1], verbose)
												 : fuse_parallel(n[1], n[0], verbose);
				}
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 1 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1");
			}
			else if (mt[n[0]].nn[e[0]] == 1 && 
					 mt[n[1]].nn[e[1]] == 0 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0" );
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 0 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0" );
			}
		}
	}
	else if (a == mt[w].g.size()) {								// at end 2
		ind1 = mt[w].g.back().ind;
		ind2 = mt[mt[w].neib[2][1]].g[mt[mt[w].neib[2][1]].end2a(mt[w].neen[2][1])].ind;
		if (mt[w].nn[2] == 2 ) {
			const auto ninds = mt[w].doubleNeibInds( 2 );
			f = true;
			for (int j=0; j<2; j++) {
				n[j] = mt[w].neib[2][ninds[j]];
				e[j] = mt[w].neen[2][ninds[j]];
			}
		}
		else if(mt[w].nn[2] == 1)
			n[0] = mt[w].neib[2][1];

		clusterisaloop = update_cl_fiss(w, 2);				// if not a loop, this increments clnum and forms a new cluster from w's end 2 neibs and beyond (excluding w itself)
		if(!clusterisaloop)
			update_gIndcl( clini );			// renumber Edge::indcl of the remaining part of the original cluster

		update_neibs(w, 2, 1, mt[w].nn[2], -1, -1, true);

		if (f && n[0] != n[1] ) {
			if (mt[n[0]].nn[e[0]] == 1 && 
				mt[n[1]].nn[e[1]] == 1 ) {
				const auto indneib0 = mt[n[0]].singleNeibInd( e[0] );
				const auto indneib1 = mt[n[1]].singleNeibInd( e[1] );
				if (mt[n[0]].neib[e[0]][indneib0] == n[1] && 
				    mt[n[0]].neen[e[0]][indneib0] == e[1] && 
				    mt[n[1]].neib[e[1]][indneib1] == n[0] && 
				    mt[n[1]].neen[e[1]][indneib1] == e[0] ) {
					update_neibs( n[0], e[0], 1, 1, -1, -1, true );
					if (e[0] == e[1])
						fuse_antiparallel(e[0], n[0], n[1], verbose);
					else
						(e[0] == 1 && e[1] == 2) ? fuse_parallel( n[0], n[1], verbose )
												 : fuse_parallel( n[1], n[0], verbose );
				}
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 1 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1" );
			}
			else if (mt[n[0]].nn[e[0]] == 1 && 
					 mt[n[1]].nn[e[1]] == 0 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0" );
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 0 ) { 
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0" );
			}
		}
	}
	basic_update();
	const auto w1 = glm[ind1];
	const auto w2 = glm[ind2];
	if (mt[w1].cl != clini && mt[w2].cl != clini ) 
		oel.exit("Error in fiss3: w1 != clini && w2 != clini" );
	if (verbose ) {
		mt[w1].print( w1, "       producing " );
		if (w1 != w2 ) mt[w2].print( w2, "             and " );
		std::cout << std::endl;
	}
	return {mt[w1].cl, mt[w2].cl};
}

}
#endif /* AbilityForFission_h */
