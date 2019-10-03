/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann,
   Helmholtz Center for Infection Research (Braunschweig, Germany).
   All Rights Reserved.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

============================================================================== */

#ifndef FUSION11_h
#define FUSION11_h

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "fusion.h"

namespace MitoSim {
	
/**
 * Reaction slot for fusion of two nodes of degree 1.
 * @tparam Ntw type of the network
 */
template<typename Ntw>
class Fusion11
	: public Fusion<1,1,Ntw> {

	friend Gillespie<Reaction,RandFactory>;

public:

	/**@brief Constructor.
	 * @param msgr Output message processor.
	 * @param ind reaction id
	 * @param netw the network
	 * @param rate rate constant
	 * @param it iteration counter
	 * @param time current time
	 */
	Fusion11(
			Msgr& msgr,
			const szt ind,
			Ntw& netw,
			const real rate,
			const ulong& it,	// const ref
			const real& time	// const ref
		)
		: Fusion<1,1,Ntw> {msgr, ind, netw, rate, it, time, name}
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

	/** The number of times this reaction was fired
	* @result the number of times this reaction was fired
	*/
	szt event_count() const noexcept final { return eventCount; }

	/** Prints the parameters
	* @param le true if new line after the output
	*/
	void print(const bool le) const override;

private:

	using Reaction::rate;
	using Reaction::msgr;
	using Fusion<1,1,Ntw>::netw;
	using Fusion<1,1,Ntw>::cc;
	using Fusion<1,1,Ntw>::update_netw_stats;
	using Fusion<1,1,Ntw>::print;

	static const std::string name;	/**< reaction name constant */
	real*	score {};				/**< current rate as seen by the Gillespie reactor */
	szt		eventCount {};			/**< number of times this reaction was fired */
	szt		propTotal {};			/**< total propensity for this reaction over all network components */

	/** Sets this reaction propensity for the whole network */
	void set_prop() noexcept final;

	/** Attaches this score to the Gillespie mechanism
	* @param a placeholder in the Gillespie object responsible for this reaction score
	*/
	void attach_score_pointer(real* a) noexcept final { score = a; };
};
	
// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename Ntw> const std::string Fusion11<Ntw>::name {"fu11"};

template<typename Ntw>
void Fusion11<Ntw>::
set_score() noexcept
{
	*score = rate * propTotal;
}

template<typename Ntw>
void Fusion11<Ntw>::
set_prop() noexcept
{
	propTotal = netw.fu11.set_prop();
}

template<typename Ntw>
void Fusion11<Ntw>::
update_prop( szt, szt ) noexcept
{
	set_prop();
}

template<typename Ntw>
void Fusion11<Ntw>::
fire() noexcept
{
	if constexpr (verbose) print(true);

	eventCount++;

	cc = netw.fu11.fire();

	update_netw_stats();
}

template<typename Ntw>
void Fusion11<Ntw>::
print( const bool le ) const
{
	Fusion<1,1,Ntw>::print(false);
	msgr.template print<false>(" score "+STR(*score));
	msgr.template print<false>(" eventCount "+STR(eventCount));
	if (le) msgr.print("\n");
}

}	// namespace MitoSim


#endif // FUSION11_h
