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

#ifndef FUSION_CANDIDATES_H
#define FUSION_CANDIDATES_H

namespace MitoD {

/**
 * Container for fusion candidate nodes.
 * Should be used only the reactions not intended for fusion to a loop segment.
 */
struct FusionCandidatesXX {

	std::vector<std::array<szt,2>> u;	/**< segment and end indexes of the 1st participant */
	std::vector<std::array<szt,2>> v;	/**< segment and end indexes of the 2nd participant */

	/** Empty the container */
	void clear() noexcept { u.clear(); v.clear(); }

	/** Add a node pair
	 * @param uc segment and end indexes of the 1st participant
	 * @param vc segment and end indexes of the 2nd participant
	*/
	void add(const std::array<szt,2>& uc,
			 const std::array<szt,2>& vc ) {
		u.emplace_back(uc);
		v.emplace_back(vc);
	}

	/** Report the number of elements
	* @result current number of candidate fusion pairs
	*/
	szt size() const noexcept { return u.size(); }

	/** Print the content out */
	void print() {
		for(szt i=0; i<size(); i++)
			print(i, false);
		std::cout << "\n";
	}

	/** Print a particular element
	 * @param i element index inside the container
	 * @param nl bool true is cr is intended
	 */
	void print( const szt i,
				bool nl=true ) {
		std::cout << " [" <<  u[i][0] << " " << u[i][1] << " + "
						  <<  v[i][0] << " " << v[i][1] << "] ";
		if (nl) std::cout << "\n";
	}
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Container for fusion candidate nodes.
 * Should be used only the reactions intended for fusion to a loop segment.
 */
struct FusionCandidatesXU {

	std::vector<std::array<szt,2>> 	u;	/**< segment and end indexes of the non-looped participant */
	std::vector<szt>	 			v;	/**< segment index of the looped participant */

	/** Empty the container */
	void clear() noexcept { u.clear(); v.clear(); }

	/** Add a node pair
	 * @param uc segment and end indexes of the non-looped participant
	 * @param vc segment index of the looped participant
	*/
	void add(std::array<szt,2> uc,
			 szt vc ) {
		u.emplace_back(uc);
		v.emplace_back(vc);
	}

	/** Report the number of elements
	 * @return current number of elements
	 */
	szt size() const noexcept { return u.size(); }
};

}	// namespace MitoD

#endif // FUSION_CANDIDATES_H
