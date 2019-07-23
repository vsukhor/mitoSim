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

#ifndef FUSION_h
#define FUSION_h

#include <string>
#include <vector>

#include "../utils/Oel.h"
#include "Reaction.h"

namespace MitoD {

/**
 * Base class for fusion reaction classes.
 */
template<typename Ntw>
class Fusion
	: public Reaction {

	friend Gillespie<Reaction>;

public:

	/** Constructor.
	 * @param oel logging facility object
	 * @param ind reaction id
	 * @param netw the network
	 * @param rate rate constant
	 * @param it iteration counter
	 * @param time current time
	 * @param verbose bool if work in verbose mode
	 * @param nodeDegree1 degree of the 1st node
	 * @param nodeDegree2 degree of the 2nd node
	 * @param srt reaction name literal
	 */
	Fusion( Oel& oel,
			const szt ind,
			Ntw& netw,
			const real rate,
			const szt& it,			// const ref
			const real& time,		// const ref
			const bool verbose,
			const uint nodeDegree1,
			const uint nodeDegree2,
			const std::string& srt
		)
		: Reaction {oel, ind, rate, it, time, verbose, "fusion", srt}
		, netw {netw}
		, rnd {netw.rnd}
		, deg1(nodeDegree1)
		, deg2(nodeDegree2)
	{}

	static constexpr auto is_active( const std::unique_ptr<Reaction>& r );

	/**
	* Populates the vector of reactions that need a score update after *this has fired
	* and initializes the propensities and effective rate
	*/
	void initialize_dependencies( const vup<Reaction>& rc ) noexcept override;

	/** Prints the parameters*/
	void print(const bool) const override;

protected:

	using Reaction::oel;

	// Convenience references
	Ntw&		 		netw;		/**< ref: the host network for this reaction */
	RandFactory& 		rnd;		/**< ref: random number factory */

	std::array<szt,2>	cc;			/**< indices of the clusters */

	void update_netw_stats() override;

private:

	using Reaction::set_score;
	using Reaction::it;

	std::vector<Reaction*>  dependents;	/**< reactions that need a score update after *this has fired */
	uint					deg1, deg2;	/**< input node degrees */

	/** Sets this reaction propensity for the whole network */
	virtual void set_prop() noexcept = 0;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> constexpr
auto Fusion<Ntw>::
is_active( const std::unique_ptr<Reaction>& r )
{
	return  r->srt == "fu1L" ||
			r->srt == "fu11" ||
			r->srt == "fu12";
}

template<typename Ntw>
void Fusion<Ntw>::
initialize_dependencies( const vup<Reaction>& rc ) noexcept
{
	for (const auto& o : rc)
		if (
			Fusion<Ntw>::is_active(o) ||
			Fission<Ntw>::is_active(o)
			)
			dependents.push_back(o.get());

	set_prop();
	set_score();
}

template<typename Ntw> 
void Fusion<Ntw>::
update_netw_stats()
{
	netw.update_books();
	
	for (auto& o : dependents) {
		o->update_prop(cc[0], cc[1]);
		o->set_score();
	}
}

template<typename Ntw>
void Fusion<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	oel.print<false>(" deg1 %d", deg1);
	oel.print<false>(" deg2 %d", deg2);
	if (le) oel.print("");
}

}

#endif /* FUSION_h */
