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

#ifndef Fission_h
#define Fission_h

#include <string>
#include <vector>

#include "../utils/Oel.h"
#include "Reaction.h"

namespace MitoD {
using namespace Utils;

template<typename>
class Fusion;

/**
 * Fission reaction class.
 */
template<typename Ntw>
class Fission
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
	 */
	Fission( Oel& oel,
			 const szt ind,
			 Ntw& netw,
			 const real rate,
			 const szt& it,		// const ref
			 const real& time,	// const ref
			 const bool& verbose
		)
		: Reaction {oel, ind, rate, it, time, verbose, "fission", name}
		, netw {netw}
		, rnd {netw.rnd}
	{
		netw.fis.verbose = verbose;
	}

	/** Sets the Gillespie score for this reaction */
	void set_score() noexcept override;

	/** Returns the Gillespie score for this reaction */
	real get_score() const noexcept override { return *score; };

	/** Updates propensity for two network components indexed in the parameters */
	void update_prop(szt,szt) noexcept override;

	/** Executes the raction event */
	void operator()() noexcept override;

	static constexpr auto is_active(const std::unique_ptr<Reaction>&) noexcept;

	/**
	* Populates the vector of reactions that need a score update after *this has fired
	* and initializes the propensities and effective rate
	*/
	void initialize_dependencies(const vup<Reaction>&) noexcept override;

	/** Returns the number of times this reaction was fired */
	szt get_eventCount() const noexcept override { return eventCount; }

	/** Prints the parameters*/
	void print(const bool) const override;

private:

	using Reaction::srt;
 	using Reaction::rate;
	using Reaction::oel;
	using Reaction::verbose;

	// Convenience references
	Ntw&					netw;	/**< ref: the network */
	RandFactory&			rnd;	/**< ref: random number factory */

	std::array<szt,2>		cc;
	real*					score {};		/**< current rate as seen by the Gillespie reactor */
	szt						eventCount {};	/**< number of times this reaction was fired */
	std::vector<Reaction*>	dependents;		/**< reactions that need a score update after *this has fired */
	static const std::string name;			/**< reaction name constant */

	/** sets this reaction propensity for the whole network */
	void set_prop() noexcept;

	/** attaches this score to the Gillespie mechanism */
	void attach_score_pointer(real*a) noexcept override { score = a; };

	void update_netw_stats() override;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<typename Ntw> const std::string Fission<Ntw>::name {"fiss"};

template<typename Ntw> constexpr
auto Fission<Ntw>::
is_active( const std::unique_ptr<Reaction>& r ) noexcept
{
	return r->srt == name;
};

template<typename Ntw>
void Fission<Ntw>::
initialize_dependencies( const vup<Reaction>& rc ) noexcept
{
	for (const auto& o : rc) {
		if (
			Fusion<Ntw>::is_active(o) ||
			Fission<Ntw>::is_active(o)
			)
			dependents.push_back(o.get());
	}
	set_prop();
	set_score();
}

template<typename Ntw> inline
void Fission<Ntw>::
set_score() noexcept
{
	*score = rate * netw.fis.get_prTotal();
}

template<typename Ntw> inline
void Fission<Ntw>::
set_prop() noexcept
{
	netw.fis.set_prop();
}

template<typename Ntw> inline
void Fission<Ntw>::
update_prop( szt c0, szt c1 ) noexcept
{
	netw.fis.update_prop(c0);
	if (c1 != c0 && c1 != huge<szt>)
		netw.fis.update_prop(c1);
}

template<typename Ntw>
void Fission<Ntw>::
update_netw_stats()
{
	for(auto& o : dependents) {
		o->update_prop(cc[0], cc[1]);
		o->set_score();
	}
}

template<typename Ntw>
void Fission<Ntw>::
operator()() noexcept
{
	if (verbose) print(true);

	eventCount++;

	cc = netw.fis();

	update_netw_stats();
}

template<typename Ntw>
void Fission<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	oel.print<false>(" score %f", *score);
	oel.print<false>(" eventCount %d", eventCount);
	if (le) oel.print("\n");
}

}

#endif /* Fission_h */
