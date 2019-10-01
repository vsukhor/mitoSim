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

#ifndef EDGE_H
#define EDGE_H

#include "utils/common/misc.h"

namespace MitoD {

/**
 * The Network Edge class.
 * Edge is a minimal structural unit of the network.
 * The class handles the tasks and properties specific to a single edge and its relation to other network components.
 * @tparam ContentT slot for specifying strucutre of the internal content the Edge can hold; currently not used
 */
template<int ContentT>
class Edge {

public:

	szt					ind {huge<szt>};	/**< index network-wide: starts from 0 */
	szt					indcl {huge<szt>};	/**< index cluster-wide: starts from 0 */
	szt					cl {huge<szt>};		/**< current cluster index */
	std::array<ulong,2>	fin {{}};			/**< contribution to fission propensity at each end */

	/** Constructor
	 * @param ind index network-wide
	 * @param indcl index cluster-wide
	 * @param cl current cluster index
	 */
	Edge(const szt ind,
		 const szt indcl,
		 const szt cl);

	/** Swap the edge ends */
	void reflect();

	/** Write the edge to a file
	 * @param ofs output file stream
	*/
	void write(std::ofstream& ofs) const;

	/** Print the edge to a stream
	 * @param os output stream
	 * @param a position inside segment
	 * @param endline flag to end current line
	*/
	void print(std::ostream& os,
			   const szt a,
			   const bool endline) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int ContentT>
Edge<ContentT>::
Edge(
	const szt ind,
	const szt indcl,
	const szt cl
	)
	: ind {ind}
	, indcl {indcl}
	, cl {cl}
{}

template<int ContentT>
void Edge<ContentT>::
reflect()
{
	std::swap(fin[0], fin[1]);
}

template<int ContentT>
void Edge<ContentT>::
write( std::ofstream &ofs ) const
{
	ofs.write(reinterpret_cast<const char*>(&ind), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&indcl), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&cl), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&fin[0]), sizeof(ulong));
	ofs.write(reinterpret_cast<const char*>(&fin[1]), sizeof(ulong));
}

template<int ContentT>
void Edge<ContentT>::
print( std::ostream& os,
	   const szt a,
	   const bool endline ) const
{
	os << "[" << a << "] ";
	os << " ind " << ind; 
	os << " indcl " << indcl;
	os << " fin " << fin[0] << " " << fin[1];
	if (endline) os << "\n";
}

}	// namespace MitoD

#endif // EDGE_H

