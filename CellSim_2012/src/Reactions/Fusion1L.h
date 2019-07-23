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

#ifndef _Fusion1L_h
#define _Fusion1L_h

#include "Fusion.h"

namespace MitoD {
	
/**
 * Reaction slot for fusion of a degree 1 node with a loop segment.
 */
template<typename Ntw>
class Fusion1L
	: public Fusion<Ntw> {

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
	 */
	Fusion1L( Oel& oel,
			  const szt ind,
			  Ntw& netw,
			  const real rate,
			  const szt& it,		// const ref
			  const real& time,		// const ref
			  const bool verbose,
			  const uint nodeDegree1,
			  const uint nodeDegree2 )
		: Fusion<Ntw> {oel, ind, netw, rate, it, time, verbose, nodeDegree1, nodeDegree2, name}
	{
		netw.fu1L.verbose = verbose;
	}

	/** Sets the Gillespie score for this reaction */
	void set_score() noexcept final;

	/** Returns the Gillespie score for this reaction */
	real get_score() const noexcept final { return *score; };

	/** Updates propensity for two network components indexed in the parameters */
	void update_prop(szt,szt) noexcept final;

	/** Executes the raction event */
	void operator()() noexcept final;

	/** Returns the number of times this reaction was fired */
	szt get_eventCount() const noexcept final { return eventCount; }

	/** Prints the parameters*/
	void print(const bool) const final;

private:

	using Reaction::rate;
	using Reaction::verbose;
	using Reaction::oel;
	using Fusion<Ntw>::cc;
	using Fusion<Ntw>::netw;
	using Fusion<Ntw>::rnd;
	using Fusion<Ntw>::update_netw_stats;
	using Fusion<Ntw>::print;

	static const std::string name;	/**< reaction name constant */
	real*	score {};			/**< current rate as seen by the Gillespie reactor */
	szt		eventCount {};		/**< number of times this reaction was fired */
	szt		propTotal {};		/**< total propensity for this reaction over all network components */

	/** Sets this reaction propensity for the whole network */
	void set_prop() noexcept final;

	/** Attaches this score to the Gillespie mechanism */
	void attach_score_pointer(real* a) noexcept final { score = a; };
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename Ntw> const std::string Fusion1L<Ntw>::name {"fu1L"};

template<typename Ntw>
void Fusion1L<Ntw>::
set_score() noexcept
{
	*score = rate * propTotal;
}

template<typename Ntw>
void Fusion1L<Ntw>::
set_prop() noexcept
{
	propTotal = netw.fu1L.set_prop();
}

template<typename Ntw>
void Fusion1L<Ntw>::
update_prop( szt, szt ) noexcept
{
	set_prop();
}

template<typename Ntw>
void Fusion1L<Ntw>::
operator()() noexcept
{
	if(verbose)	print(true);

	eventCount++;
	
	cc = netw.fu1L();

	update_netw_stats();
}

template<typename Ntw>
void Fusion1L<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	oel.print0(" score "+STR(*score));
	oel.print0(" eventCount "+STR(eventCount));
	if (le) oel.print("\n");
}

}

#endif	/* _Fusion1L_h */
