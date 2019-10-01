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

#ifndef SEGMENT_H
#define SEGMENT_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "edge.h"

namespace MitoD {

/**
 * Segment class of the Network.
 * Segment is a sequence of edges linked linearly (without branches).
 * Segment ends may form branching sites, where it is connected to other segments.
 * A segment not connected to other segments or
 * a complete collection of segments connected to each other form
 * a disconnected network component (aka 'cluster').
 * The class handles the tasks and properties specific to a single segment
 * and its relation to other network components.
 */
template<szt>
class Segment {};

/**
 * Segment class specification for max node degree equal to 3.
 * Segment is a sequence of edges linked linearly (without branches).
 * Segment ends may form branching sites, where it is connected to other segments.
 * A segment not connected to other segments or
 * a complete collection of segments connected to each other form
 * a disconnected network component (aka 'cluster').
 * The class handles the tasks and properties specific to a single segment
 * and its relation to other network components.
 */
template<>
class Segment<3> {

public:

	static constexpr szt maxDeg {3};	/**< maximal node degree allowed */

	using EdgeT = Edge<maxDeg>;
	using thisT = Segment<maxDeg>;

	Msgr&								msgr;		/**< logging facility */
	std::vector<EdgeT>					g;			/**< the edges */
	szt									cl {};		/**< cluster index */
	std::array<szt,maxDeg>				nn {{}};	/**< number of neighbours */
	std::array<std::vector<szt>,maxDeg>	neig;		/**< neighbours */
	std::array<std::vector<szt>,maxDeg>	neen;		/**< neighbour ends */

	/** Constructor
	 * @param msgr logging facility object
	 */
	explicit Segment(Msgr& msgr);

	/** Constructor
	 * @param segmass segment mass
	 * @param cl index of subnetwork to which the sebment belongs
	 * @param mtmass total mass of the network
	 * @param ei index of the last edge in this segment
	 * @param msgr logging facility object
	 */
	explicit Segment(
		  const szt segmass,
		  const szt cl,
		  szt& mtmass, 
		  szt& ei, 
		  Msgr& msgr );

	/** Reflects the vector containing the segment edges */
	void reflect_g();

	/** Change disconnected network component-related indexes of the segment edges,
			keeping the segment index unoltered
	 * @param newcl new disconnected component index
	 * @param initind starting edge index in the current disconnected network component
	 */
	szt set_gCl(const szt newcl, const szt initind);

	/** Change disconnected network component-related indexes of the segment edges, and the the segment itself
	 * @param newcl new disconnected component index
	 * @param initind starting edge index in the current disconnected network component
	 * @return the last edge index in the current disconnected network component
	 */
	szt setCl(const szt newcl, const szt initind);

	/** Converts the segment end index of the boundary edge to internal position in the segment
	 * @param e segment end index
	 * @return internal position
	 * @return the last edge index in the current disconnected network component
	 */
	constexpr szt end2a(const szt& e) const;

	/** Determine if the segment has one free end
 	 * @return the end index if true
	 */
	constexpr szt has_one_free_end() const;

	/** Neigbour indexe at a segment end
	 * @param e segment end
 	 * @return the neighbour index
	 */
	szt single_neig_index(const szt& e) const;

	/** Neigbour indexes at a segment end
	 * @param e segment end
 	 * @return the neighbour indexes
	 */
	std::vector<szt> double_neig_indexes(const szt& e) const;

	/** Report if the segment is looped onto itself */
	constexpr bool is_cycle() const;

	/** Report the number of nodes of a given degree.
	 * @param deg node degree
	 * @return the number of nodes
	 */
	szt num_nodes(const szt& deg) const;

	/** Report the segment length measured in edges.
	 * @return the segment length measured in edges
	 */
	szt length() const noexcept { return g.size(); }

	ulong set_end_fin(const szt a);
	ulong set_bulk_fin(const szt a);

	void print( const szt w, const std::string& tag, const szt at=huge<szt> ) const;
	void print( std::ostream& os, const szt w, const std::string& tag, const szt at=huge<szt> ) const;

	/** Write the segment to a binary file.
	 * @param ofs std::ofstream to write to
	 */
	void write(std::ofstream& ofs) const;

private:

	/**
	* Inserts an edge imediately after g[a] making it g[a+1].
	* @return the pointer to the newly inserted edge.
	* @param a position of the edge preceding the newly inserted one relative to the segment end 1
	* @param p edge to be inserted
	*/
	EdgeT* increment_length( const long a, EdgeT p );	// inserts a particle after g[a] making it g[a+1]

	/** Initialises the neig vectors at both ends. */
	void init_ends();
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Segment<3>::
Segment( Msgr& msgr )
	: msgr {msgr}
{
	init_ends();
}

Segment<3>::
Segment(
	  const szt segmass,
	  const szt cl,
	  szt& mtmass, 		// var ref
	  szt& ei, 			// var ref
	  Msgr& msgr			// var ref
	)
	: msgr {msgr}
	, cl {cl}
{
	init_ends();

	for (szt a=0; a<segmass; a++)
		increment_length(long(a-1), {ei++, a, cl});

	mtmass += segmass;
}

inline
void Segment<3>::
init_ends()
{
	neig[1].resize(maxDeg);
	neig[2].resize(maxDeg);

	neen[1].resize(maxDeg);
	neen[2].resize(maxDeg);
}

inline // inserts a particle imediately after g[a] making it g[a+1]
typename Segment<3>::EdgeT* Segment<3>::
increment_length( const long a, Segment<3>::EdgeT p )
{
	g.insert(g.begin()+a+1, std::move(p));
	return &g[a+1];
}

inline
void Segment<3>::
reflect_g()
{
	std::reverse(g.begin(), g.end());
	for (auto& o : g)
		o.reflect();
}

inline
szt Segment<3>::
set_gCl( const szt newcl, const szt initind )
{
	for (szt i=0; i<g.size(); i++) {
		g[i].cl = newcl;
		g[i].indcl = initind + i;
	}
	
	return initind + (szt)g.size();
}

inline
szt Segment<3>::
setCl( const szt newcl, const szt initind )
{
	cl = newcl;
	return set_gCl(newcl, initind);
}

constexpr
szt Segment<3>::
end2a( const szt& e ) const 
{
	return (e == 1) ? 0 : (szt)g.size()-1;
}

constexpr
szt Segment<3>::
has_one_free_end() const 		// return the end index if true
{
	if (	!nn[1] &&  nn[2])  	return 1;
	else if (nn[1] && !nn[2])	return 2;
	else						return 0;
}

inline
szt Segment<3>::
single_neig_index( const szt& e ) const
{
	for (szt i=1; i<=nn[e]; i++) 
		if (neig[e][i])
			return i;
			
	return huge<szt>;
}

inline
std::vector<szt> Segment<3>::
double_neig_indexes( const szt& e ) const
{
	XASSERT(nn[e] == 2, "Error in Mito::double_neig_indexes: nn[e] != 2 in cluster "+STR(cl)+"\n");

	std::vector<szt> neigInds(nn[e]);
	szt j {};
	for (szt i=1; i<=nn[e]; i++)
		if (neig[e][i])
			neigInds[j++] = i;
	
	return neigInds;
}

constexpr
bool Segment<3>::
is_cycle() const
{
	return nn[1] == 1 && 
		   nn[2] == 1 && 
		   neig[1][single_neig_index(1)] == neig[2][single_neig_index(2)];
}

szt Segment<3>::
num_nodes( const szt& deg ) const				// deg = 1, 2, 3
{														
	if (deg == 1) {														// count nodes of degree 1
		if (      nn[1] &&  nn[2]) 	return 0;
		else if (!nn[1] && !nn[2]) 	return 2;
		else					 	return 1;
	}
	else if (deg == 2)
		return nn[1] && nn[2] && is_cycle() ? g.size() : g.size() - 1;	// count nodes of degree 2
	else if (deg == 3) {												// count nodes of degree 3
		if (     nn[1] == 2 && nn[2] == 2) return 2;
		else if (nn[1] == 2 || nn[2] == 2) return 1;
		else if (nn[1] != 2 && nn[2] != 2) return 0;
	}
	else msgr.exit("Error in Mito::num_nodes(). Not implemented deg", deg);
	return huge<szt>;
}

inline
ulong Segment<3>::
set_end_fin( const szt e )
{
	auto& f {g[end2a(e)].fin};
	f[e-1] = nn[e] ? 1 : 0;

	return f[e-1];
}


inline
ulong Segment<3>::
set_bulk_fin( const szt a )
{
	g[a].fin[1] = g[a+1].fin[0] = 1;

	return g[a].fin[1];
}

void Segment<3>::
print( const szt w,
	   const std::string& tag, 
	   const szt at ) const
{
	if (msgr.so) print(*msgr.so, w, tag, at);
	if (msgr.sl) print(*msgr.sl, w, tag, at);
}

void Segment<3>::
print( std::ostream& os, 
	   const szt w,
	   const std::string& tag, 
	   const szt at ) const
{
	os << "\t" << tag << w;
	if (at == huge<szt>) os << "(of ";
	else os << "(at " << at << " of ";
	os << g.size() << ") [ ";	
	for (szt i=1; i<=nn[1]; i++) os << neig[1][i] << " ";
	os << "] { ";		
	for (szt i=1; i<=nn[1]; i++) os << neen[1][i] << " ";	
	os << "} [ ";
	for (szt i=1; i<=nn[2]; i++) os << neig[2][i] << " ";
	os << "] { ";	
	for (szt i=1; i<=nn[2]; i++) os << neen[2][i] << " ";	
	os << "} " << cl;	
#ifdef PRINT_EDGES
	os << std::endl;
	for (szt i=0; i<g.size(); i++) 
		g[i].print(os, i, 1);
#else
	os << " len " << g.size();
#endif
	os << std::endl;
}

void Segment<3>::
write( std::ofstream& ofs ) const
{
	szt len = (szt)g.size();
	ofs.write(reinterpret_cast<const char*>(&len), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&cl), sizeof(szt));
	
	ofs.write(reinterpret_cast<const char*>(&nn[1]), sizeof(szt));
	
	for (szt j=1; j<=nn[1]; j++) {
		ofs.write(reinterpret_cast<const char*>(&neig[1][j]), sizeof(szt));
		ofs.write(reinterpret_cast<const char*>(&neen[1][j]), sizeof(szt));
	}
	
	ofs.write(reinterpret_cast<const char*>(&nn[2]), sizeof(int));
	
	for (szt j=1; j<=nn[2]; j++) {
		ofs.write(reinterpret_cast<const char*>(&neig[2][j]), sizeof(szt));
		ofs.write(reinterpret_cast<const char*>(&neen[2][j]), sizeof(szt));
	}
	for (szt j=0; j<len; j++)
		g[j].write(ofs);
}

}	// namespace MitoD

#endif // SEGMENT_H
