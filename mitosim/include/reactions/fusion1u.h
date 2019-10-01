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

#ifndef FUSION1U_H
#define FUSION1U_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "fusion.h"

namespace MitoD {
	
/**
 * Reaction slot for fusion of a degree 1 node with a loop segment.
 * @tparam Ntw type of the network
 */
template<typename Ntw>
class Fusion1U
	: public Fusion<1,2,Ntw> {

	friend Gillespie<Reaction,RandFactory>;

public:

	/** Constructor.
	 * @param msgr logging facility object
	 * @param ind reaction id
	 * @param netw the network
	 * @param rate rate constant
	 * @param it iteration counter
	 * @param time current time
	 */
	Fusion1U( Msgr& msgr,
			  const szt ind,
			  Ntw& netw,
			  const real rate,
			  const ulong& it,		// const ref
			  const real& time		// const ref
			)
		: Fusion<1,2,Ntw> {msgr, ind, netw, rate, it, time, name}
	{}

	/** Sets the Gillespie score for this reaction */
	void set_score() noexcept final;

	/** Gillespie score for this reaction
	* @result total weight of this reaction in the simulation set
	*/
	real get_score() const noexcept final { return *score; };

	/** Updates propensity for a pair of network components
	* @param c1 index of the 1st component to update
	* @param c2 index of the 2nd component to update
	*/
	void update_prop(const szt c1, const szt c2) noexcept final;

	/** Executes the raction event */
	void fire() noexcept final;

	/** Returns the number of times this reaction was fired */
	szt event_count() const noexcept final { return eventCount; }

	/** Prints the parameters
	* @param le true if new line after the output
	*/
	void print(const bool le) const override;

private:

	using Reaction::rate;
	using Reaction::msgr;
	using Fusion<1,2,Ntw>::cc;
	using Fusion<1,2,Ntw>::netw;
	using Fusion<1,2,Ntw>::rnd;
	using Fusion<1,2,Ntw>::update_netw_stats;
	using Fusion<1,2,Ntw>::print;

	static const std::string name;	/**< reaction name constant */
	real*	score {};			/**< current rate as seen by the Gillespie reactor */
	szt		eventCount {};		/**< number of times this reaction was fired */
	szt		propTotal {};		/**< total propensity for this reaction over all network components */

	/** Sets this reaction propensity for the whole network */
	void set_prop() noexcept final;

	/** Attaches this score to the Gillespie mechanism
	* @param a placeholder in the Gillespie object responsible for this reaction score
	*/
	void attach_score_pointer(real* a) noexcept final { score = a; };
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename Ntw> const std::string Fusion1U<Ntw>::name {"fu1L"};

template<typename Ntw>
void Fusion1U<Ntw>::
set_score() noexcept
{
	*score = rate * propTotal;
}

template<typename Ntw>
void Fusion1U<Ntw>::
set_prop() noexcept
{
	propTotal = netw.fu1L.set_prop();
}

template<typename Ntw>
void Fusion1U<Ntw>::
update_prop( const szt, const szt ) noexcept
{
	set_prop();
}

template<typename Ntw>
void Fusion1U<Ntw>::
fire() noexcept
{
	if constexpr (verbose) print(true);

	eventCount++;
	
	cc = netw.fu1L.fire();

	update_netw_stats();
}

template<typename Ntw>
void Fusion1U<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	msgr.template print<false>(" score "+STR(*score));
	msgr.template print<false>(" eventCount "+STR(eventCount));
	if (le) msgr.print("\n");
}

}	// namespace MitoD

#endif	// FUSION1U_H 
