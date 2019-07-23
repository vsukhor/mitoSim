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

#ifndef AbilityForFission_h
#define AbilityForFission_h

#include <vector>
#include <algorithm>

#include "utils/Misc.h"
#include "utils/Oel.h"

#include "CoreTransformer.h"

namespace MitoD {

template<typename> class Fission;

/**
 * The AbilityForFission class.
 * Adds node type-specific fission capability and updates the network for it.
 */

template<typename Mt>
class AbilityForFission
	: public CoreTransformer<Mt> {

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
	using CoreTransformer<Mt>::copy_neigs;
	using CoreTransformer<Mt>::update_neigs;
	using CoreTransformer<Mt>::fuse_antiparallel;
	using CoreTransformer<Mt>::fuse_parallel;
	using CoreTransformer<Mt>::update_gIndcl;

	/** Constructor
	 * @param cfg configuration object
	 * @param oel logging facility object
	 */
	explicit AbilityForFission(const Config& cfg, Oel& oel);

	/** Performs fission of a segment
	 * @param w segment index
	 * @param a division positon inside the segment
	 * @param verb verbasity
	 */
	std::array<szt,2> fiss( const szt w, const szt a, const bool verb);

	/**
	 * Divides the segment at a node of degree 2.
	 * @param w a global segment index
	 * @param a the node position inside the segment
	 * @param verb verbasity
	 */
	std::array<szt,2> fiss2(const szt w, const szt a, const bool verb);

	/**
	 * Divides the segment at a node of degree 3.
	 * @param w a global segment index
	 * @param a the node position inside the segment
	 * @param verb verbasity
	 */
	std::array<szt,2> fiss3(const szt w, const szt a, const bool verb);

private:

	std::vector<szt>  vis;	/**< auxiliary, indicating a visited or not status of segments during the graph search */

	/** Updates network component to account for changes resulting from its division.
	 * Returns a flag indicating if the division produces a pair of disconnected components
	 * @param w segment index
	 * @param e segment end
	 */
	bool update_cl_fiss(const szt w, const szt e);

	/** Depth-first search of the network graph.
	 * @param w1 initial segment index
	 * @param e1 initial segment end
	 * @param w2 final segment index
	 * @param e2 final segment end
	 * @return true if there is a connection between {w1,e1} and {w2,e2}
	 */
	bool dfs(const szt w1, const szt e1, const szt w2, const szt e2);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
AbilityForFission<Mt>::
AbilityForFission(
		const Config& cfg,
		Oel& oel
	)
	: CoreTransformer<Mt> {cfg, oel}

{}

template<typename Mt>
bool AbilityForFission<Mt>::
update_cl_fiss( const szt w, const szt e )
{
	vis.resize(mtnum+2);
	std::fill(vis.begin()+1, vis.begin()+mtnum+1, 0);		// nothing is visited in the beginning of the search

	const szt oe {(e == 1) ? szt(2) : szt(1)};

	const bool isCycle = dfs(w, e, w, oe);
	if (!isCycle) {
		clnum++;
		szt clind {};
		for (szt i=1; i<=mtnum; i++)
			if (vis[i])
				clind = mt[i].setCl(clnum-1, clind);
	}
	return isCycle;
}

template<typename Mt>
bool AbilityForFission<Mt>::
dfs( const szt w1, const szt e1,
	 const szt w2, const szt e2 )
{
	for (szt i=1; i<=mt[w1].nn[e1]; i++) {
		const auto cn {mt[w1].neig[e1][i]};
		const auto ce {mt[w1].neen[e1][i]};
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

// a is counted from 1 and is the last element to remain in the old segment
template<typename Mt>
std::array<szt,2> AbilityForFission<Mt>::
fiss( const szt w, const szt a, const bool verbose )
{
	// node 3 -> nodes 2+1 and node 2 -> nodes 1+1 in pure loops
	if (        a && a < mt[w].g.size())
		return fiss2(w, a, verbose);

	// node 2 -> nodes 1+1; cuts between g[a-1] and g[a]
	else if ((!a && mt[w].nn[1] <= 2) ||
			 (a == mt[w].g.size() && mt[w].nn[2] <= 2))
		return fiss3(w, a, verbose);

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
	if (verbose)
		mt[w].print(w, "fission2:  ", a);		// cuts between g[a-1] and g[a]

	const auto clini {mt[w].cl};

	const auto ind1 {mt[w].g[a-1].ind};
	const auto ind2 {mt[w].g[a].ind};

	bool inCycle {};
	mt[w].nn[2] ? (inCycle = update_cl_fiss(w, 2))
				: clnum++;

	mt.emplace_back(oel);
	++mtnum;

	std::move(mt[w].g.begin()+a, mt[w].g.end(), std::back_inserter(mt[mtnum].g));
	mt[w].g.erase(mt[w].g.begin()+a, mt[w].g.end());

	mt[mtnum].nn[1] = 0;

	this->copy_neigs(w, 2, mtnum, 2);

	mt[mtnum].cl = inCycle
				 ? mt[w].cl
				 : clnum - 1;

	if (!inCycle) {
		update_gIndcl(mt[w].cl);		// renumber Edge::indcl of the remaining part of the original cluster
		update_gIndcl(clnum-1);			// renumber Edge::indcl of the newly formed cluster
	}

	mt[w].nn[2] = 0;

	bool isSelfLooped {};
	if (1 == mt[w].nn[1] &&
		1 == mt[mtnum].nn[2] && 
		w == mt[mtnum].neig[2][mt[mtnum].single_neig_index(2)] && 
		mtnum == mt[w].neig[1][mt[w].single_neig_index(1)]) {		// before fission 'w' was looped into itself
		isSelfLooped = true;
		if (verbose) {
			mt[w].print(w, "       transiently producing ", -1);
			mt[mtnum].print(mtnum, "             and ", -1);
			std::cout << std::endl;
		}
		update_neigs(w, 1, 1, 1, -1, -1, true);
		fuse_parallel(w, mtnum, verbose);
	}

	basic_update();
	const auto w1 {glm[ind1]};
	const auto w2 {glm[ind2]};

	XASSERT(mt[w1].cl == clini || mt[w2].cl == clini, "Error in fiss3: mt[w1].cl != clini && mt[w2].cl != clini\n");
	if(verbose) {
		mt[w1].print(w1, "       producing ");
		if (isSelfLooped)
			oel.print("from a segment looped into itself");
		else mt[w2].print(w2, "             and ");
		std::cout << std::endl;
	}

	return {mt[w].cl, mt[mtnum].cl};
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
	bool inCycle {};
	bool f {};
	szt n[2], e[2];
	auto ind1 {huge<szt>};
	auto ind2 {huge<szt>};
	if (!a) {															// at end 1
		ind1 = mt[w].g.front().ind;
		ind2 = mt[mt[w].neig[1][1]].g[mt[mt[w].neig[1][1]].end2a(mt[w].neen[1][1])].ind;
		if (mt[w].nn[1] == 2) {
			const auto ninds = mt[w].double_neig_indexes(1);
			f = true;
			for (int j=0; j<2; j++) {
				n[j] = mt[w].neig[1][ninds[j]];
				e[j] = mt[w].neen[1][ninds[j]];
			}
		}
		else if (mt[w].nn[1] == 1)
			n[0] = mt[w].neig[1][1];

		// if not a cycle, this increments clnum and forms a new cluster from w's end 1 neigs and beyond (excluding w itself)
		inCycle = update_cl_fiss(w, 1);
		if(!inCycle)
			update_gIndcl(clini);			// renumber Edge::indcl of the remaining part of the original cluster

		update_neigs(w, 1, 1, mt[w].nn[1], -1, -1, true);

		if (f && n[0] != n[1]) {
			if (mt[n[0]].nn[e[0]] == 1 && 
			    mt[n[1]].nn[e[1]] == 1) {

				const auto indneig0 = mt[n[0]].single_neig_index(e[0]);
				const auto indneig1 = mt[n[1]].single_neig_index(e[1]);

				if (mt[n[0]].neig[e[0]][indneig0] == n[1] && 
				    mt[n[0]].neen[e[0]][indneig0] == e[1] && 
				    mt[n[1]].neig[e[1]][indneig1] == n[0] && 
				    mt[n[1]].neen[e[1]][indneig1] == e[0]) {

					update_neigs(n[0], e[0], 1, 1, -1, -1, true);

					if (e[0] == e[1])
						fuse_antiparallel(e[0], n[0], n[1], verbose);
					else
						(e[0] == 1 && e[1] == 2) ? fuse_parallel(n[0], n[1], verbose)
												 : fuse_parallel(n[1], n[0], verbose);
				}
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 1) {
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print(n[0], "       n[0] ", -1);
				mt[n[1]].print(n[1], "       n[1] ", -1);
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1");
			}
			else if (mt[n[0]].nn[e[0]] == 1 && 
					 mt[n[1]].nn[e[1]] == 0) {
				mt[w].print(w, "       w ", -1);
				mt[n[0]].print(n[0], "       n[0] ", -1);
				mt[n[1]].print(n[1], "       n[1] ", -1);
				oel.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0");
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 0) {
				mt[w].print(w, "       w ", -1);
				mt[n[0]].print( n[0], "       n[0] ", -1 );
				mt[n[1]].print( n[1], "       n[1] ", -1 );
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0");
			}
		}
	}
	else if (a == mt[w].g.size()) {								// at end 2
		ind1 = mt[w].g.back().ind;
		ind2 = mt[mt[w].neig[2][1]].g[mt[mt[w].neig[2][1]].end2a(mt[w].neen[2][1])].ind;
		if (mt[w].nn[2] == 2) {
			const auto ninds = mt[w].double_neig_indexes(2);
			f = true;
			for (int j=0; j<2; j++) {
				n[j] = mt[w].neig[2][ninds[j]];
				e[j] = mt[w].neen[2][ninds[j]];
			}
		}
		else if(mt[w].nn[2] == 1)
			n[0] = mt[w].neig[2][1];

		// if not a cycle, this increments clnum and forms a new cluster from w's end 2 neigs and beyond (excluding w itself)
		inCycle = update_cl_fiss(w, 2);
		if(!inCycle)
			update_gIndcl(clini);			// renumber Edge::indcl of the remaining part of the original cluster

		update_neigs(w, 2, 1, mt[w].nn[2], -1, -1, true);

		if (f && n[0] != n[1]) {
			if (mt[n[0]].nn[e[0]] == 1 && 
				mt[n[1]].nn[e[1]] == 1 ) {
				const auto indneig0 = mt[n[0]].single_neig_index( e[0]);
				const auto indneig1 = mt[n[1]].single_neig_index( e[1]);
				if (mt[n[0]].neig[e[0]][indneig0] == n[1] && 
				    mt[n[0]].neen[e[0]][indneig0] == e[1] && 
				    mt[n[1]].neig[e[1]][indneig1] == n[0] && 
				    mt[n[1]].neen[e[1]][indneig1] == e[0]) {
					update_neigs(n[0], e[0], 1, 1, -1, -1, true);
					if (e[0] == e[1])
						fuse_antiparallel(e[0], n[0], n[1], verbose);
					else
						(e[0] == 1 && e[1] == 2) ? fuse_parallel(n[0], n[1], verbose)
												 : fuse_parallel(n[1], n[0], verbose);
				}
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 1) {
				mt[w].print(w, "       w ", -1);
				mt[n[0]].print(n[0], "       n[0] ", -1);
				mt[n[1]].print(n[1], "       n[1] ", -1);
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 1");
			}
			else if (mt[n[0]].nn[e[0]] == 1 && 
					 mt[n[1]].nn[e[1]] == 0) {
				mt[w].print(w, "       w ", -1);
				mt[n[0]].print(n[0], "       n[0] ", -1);
				mt[n[1]].print(n[1], "       n[1] ", -1);
				oel.exit("mt[n[0]].nn[e[0]] == 1 && mt[n[1]].nn[e[1]] == 0");
			}
			else if (mt[n[0]].nn[e[0]] == 0 && 
					 mt[n[1]].nn[e[1]] == 0) {
				mt[w].print( w, "       w ", -1 );
				mt[n[0]].print(n[0], "       n[0] ", -1);
				mt[n[1]].print(n[1], "       n[1] ", -1);
				oel.exit("mt[n[0]].nn[e[0]] == 0 && mt[n[1]].nn[e[1]] == 0");
			}
		}
	}
	basic_update();
	const auto w1 = glm[ind1];
	const auto w2 = glm[ind2];
	XASSERT(mt[w1].cl == clini || mt[w2].cl == clini, "Error in fiss3: mt[w1].cl != clini && mt[w2].cl != clini\n");
	if (verbose) {
		mt[w1].print(w1, "       producing ");
		if (w1 != w2) mt[w2].print(w2, "             and ");
		std::cout << std::endl;
	}
	return {mt[w1].cl, mt[w2].cl};
}

}
#endif /* AbilityForFission_h */
