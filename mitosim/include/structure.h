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

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

namespace MitoD {

/**
 * The Structure class.
 * Encapsulates the major structure-related proterties, such a  collection of the Segments,
 * but is unaware of any dynamics.
 * Forms base for clases adding the network reconfiguration dynamics.
  * @tparam Mt type of the Edge forming the network
*/
template<typename Mt>
class Structure {

public:

	using Reticulum = std::vector<Mt>;

	vec3<szt>					clagl;		/**< edge adjacency_lists per cluster */

	std::vector<szt>			glm;		/**< mapping of the edge indexes to segment indexes */
	std::vector<szt>			gla;		/**< mapping of the edge indexes to element index inside segments */

	Reticulum					mt;			/**< the segments */

	std::array<szt,Mt::maxDeg>	nn {{}};	/**< total number of nodes by node degree */

	szt							mtnum {};	/**< actual number of segments */
	szt							clnum {};	/**< actual number of clusters */
	szt							mtmass {};	/**< current number of edges */

	vec2<szt>					clmt;		/**< segment indices segregated into clusters: clmt - total */
	std::vector<szt>			cls;		/**< cluster sizes measured in edges */

	// segment indices of the corresponding end degrees
	std::vector<szt>				mt11;	/**< indexes of disconnected segments not looped onto itself */
	std::vector<szt>				mtc11;	/**< indexes of disconnected segments not looped onto itself */

	std::vector<szt>				mt22;	/**< indexes of disconnected looped segments */
	std::vector<szt>				mtc22;	/**< indexes of disconnected looped segments */

	std::vector<szt>				mt33;	/**< indexes of segments between nodes of degree 3 and 3: all together */
	vec2<szt>						mtc33;	/**< indexes of segments between nodes of degree 3 and 3: sorted into clusters */
	
	std::vector<std::array<szt,2>>	mt13;	/**< {index,end} pairs for segments between nodes of degree 1 and 3: all together */
	vec2<std::array<szt,2>>			mtc13;	/**< {index,end} pairs for segments between nodes of degree 1 and 3: sorted into clusters */

	Msgr&		msgr;				/**< logging facility */
	static constexpr szt minLoopLength {2};	/**< minimal length of a segment that can bend into a cycle */

	/** Constructor.
	 * @param msgr logging facility object
	 */
	explicit Structure(Msgr& msgr);

	/** Updates internal data */
	void basic_update() noexcept;

	/** Updates internal data */
	void update_adjacency() noexcept;

	/** Updates internal vectors */
	void update_structure() noexcept;

	/** Initializes or updates glm and gla vectors */
	void make_indma() noexcept;

	/** Initializes or updates adjacency list.
	 * @param ic disconnected network component index
	 * @param a the adjacency list
	 */
	void make_adjacency_list_edges( const szt ic, vec2<szt>& a ) noexcept;

	/** Populates 'mt??', 'mtc??', 'nn' and 'clmt' vectors */
	void populate_cluster_vectors() noexcept;

	/** Updates 'nn' for the specific node degree
	 * @param deg node degree to consider
	 */
	void update_nn(const szt deg) noexcept;

	/** Updates 'nn' for all node degrese */
	void update_nn() noexcept;

	/** Prints the network components using prefix specified
	* @param tag prefix
	*/
	void print_mitos(const std::string& tag) const;

	/** Prints the network components to a text-foramtted stream
	* @param ofs output stream
	*/
	void print(std::ostream& ofs) const;

private:

	vec2<bool> clvisited;	/** temporary auxiliary field */
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
Structure<Mt>::
Structure(
		Msgr& msgr
	)
	: msgr {msgr}
{}

template<typename Mt> inline
void Structure<Mt>::
update_adjacency() noexcept
{
	if (clagl.size() < clnum) {
		clagl.resize(clnum);
		clvisited.resize(clnum);
	}
	for (szt c=0; c<clnum; c++)
		make_adjacency_list_edges(c, clagl[c]);
}

template<typename Mt> inline
void Structure<Mt>::
basic_update() noexcept
{
	make_indma();
	populate_cluster_vectors();
}

template<typename Mt> inline
void Structure<Mt>::
update_structure() noexcept
{
	basic_update();
	update_adjacency();
}

template<typename Mt> inline
void Structure<Mt>::
make_indma() noexcept
{
	cls.resize(clnum);
	std::fill(cls.begin(), cls.end(), 0);	// cluster size, # of edges
	for (szt j=1; j<=mtnum; j++)
		cls[mt[j].cl] += mt[j].g.size();	

	glm.resize(mtmass);
	gla.resize(mtmass);
	for (szt j=1; j<=mtnum; j++)
		for (szt k=0; k<mt[j].g.size(); k++) {
			const auto& g {mt[j].g[k]};
			glm[g.ind] = j;
			gla[g.ind] = k;
		}
}
template<typename Mt>
void Structure<Mt>::
make_adjacency_list_edges( const szt c, vec2<szt>& a ) noexcept
{
	clvisited[c].resize(cls[c]);

	a.resize(cls[c]);
	for (auto& o : a) o.clear();

	for (const auto& j : clmt[c])
		for (szt k=0; k<mt[j].g.size(); k++) {
			const auto ind {mt[j].g[k].indcl};
			if (k == 0) {
				for (szt e=1; e<=mt[j].nn[1]; e++) {		// connection backwards: only other segments might be found
					const auto w2 {mt[j].neig[1][e]};
					const auto a2 {mt[w2].end2a(mt[j].neen[1][e])};
					a[ind].push_back( mt[w2].g[a2].indcl );
				}
				if (mt[j].g.size() == 1)					// connection forwards: to other segment
					for (szt e=1; e<=mt[j].nn[2]; e++) {
						const auto w2 {mt[j].neig[2][e]};
						const auto a2 {mt[w2].end2a(mt[j].neen[2][e])};
						a[ind].push_back(mt[w2].g[a2].indcl);
					}
				else {										// connection forwards: to the same segment
					a[ind].push_back( mt[j].g[k+1].indcl );
				}
			}
			else if (k == mt[j].g.size()-1) {				// but not  a1 == 0  =>  mt[m1].g.size() > 1
				a[ind].push_back(mt[j].g[k-1].indcl);		// connection backwards: to the same segment
				for (szt e=1; e<=mt[j].nn[2]; e++) {		// connection forwards: to other segment
					const auto w2 {mt[j].neig[2][e]};
					const auto a2 {mt[w2].end2a(mt[j].neen[2][e])};
					a[ind].push_back(mt[w2].g[a2].indcl);
				}
			}
			else {											// edge in the bulk: a1 != 1 && a1 != mt[m1].g.size()
				a[ind].push_back(mt[j].g[k-1].indcl);		// connection backwards: to the same segment
				a[ind].push_back(mt[j].g[k+1].indcl);		// connection forwards: to the same segment
			}
		}
}

template<typename Mt>
void Structure<Mt>::
populate_cluster_vectors() noexcept
{
	mt11.clear();	mtc11.resize(clnum);		std::fill(mtc11.begin(),  mtc11.end(),  huge<szt>);
	mt22.clear();	mtc22.resize(clnum);		std::fill(mtc22.begin(),  mtc22.end(),  huge<szt>);
	mt33.clear();	mtc33.resize(clnum);		for (auto& o : mtc33) o.clear();
	mt13.clear();	mtc13.resize(clnum);		for (auto& o : mtc13) o.clear();

	nn = {{zero<szt>}};
	clmt.resize(clnum);		for (auto& o : clmt) o.clear();		// # of segments
	
	for (szt j=1; j<=mtnum; j++) {
		const auto& m {mt[j]};
		clmt[m.cl].push_back(j);						// mitochondrial indexes clusterwise
		nn[1] += m.num_nodes(2);

		const auto e {m.has_one_free_end()};
		if (e) {
			const szt oe {e == 1 ? static_cast<szt>(2) : static_cast<szt>(1)};
			nn[0]++;
			if (m.nn[oe] == 2) {
				const std::array<szt,2> je {j, e};
				mtc13[m.cl].emplace_back(je);			// segment index, free end index
				mt13.emplace_back(je);
				nn[2]++;
			}
		}
		else if (m.nn[1] == 0 && m.nn[2] == 0) {
			mtc11[m.cl] = j;							// it is a separate segment since it has two free ends
			mt11.push_back(j);							// it is a separate segment since it has two free ends
			nn[0] += 2;
		}
		else if (m.is_cycle()) {
			mtc22[m.cl] = j;
			mt22.push_back(j);							// it is a separate segment since it has two free ends
		}
		else if (m.nn[1] == 2 && m.nn[2] == 2) {
			mtc33[m.cl].push_back(j);
			mt33.push_back(j);
			nn[2] += 2;
		}
		else {
			; XASSERT(false, "Error in populate_cluster_vectors: failed classification for "+STR(j)+"\n");
		}
	}
	nn[2] /= 3;
}

template<typename Mt>
void Structure<Mt>::
update_nn( const szt deg ) noexcept
{
	auto count_nodes = [&](const szt deg) noexcept	{
		szt k {};
		for (szt i=1; i<=mtnum; i++)
			k += mt[i].num_nodes(deg);
		if (deg == 3) return k/3;
		else		  return k;
	};

	nn[deg-1] = count_nodes(deg);
}

template<typename Mt>
void Structure<Mt>::
update_nn() noexcept
{
	update_nn(1);
	update_nn(2);
	update_nn(3);
}

template<typename Mt>
void Structure<Mt>::
print_mitos( const std::string& tag ) const
{
	for (szt j=1; j<=mtnum; j++)
		mt[j].print(j, tag, -1);
	msgr.print("");
}
template<typename Mt>
void Structure<Mt>::
print( std::ostream& ofs ) const
{
	ofs << " X ";
	for (const auto o : nn)
		ofs << o << " ";
	ofs <<  "m11 " << mt11.size()
		<< " m22 " << mt22.size()
		<< " m33 " << mt33.size()
		<< " m13 " << mt13.size()
		<< " mtm " << mtmass
		<< " mtn " << mtnum
		<< " cln " << clnum;
}

}	// namespace MitoD

#endif // STRUCTURE_H
