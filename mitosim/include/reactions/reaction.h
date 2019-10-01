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

#ifndef REACTION_h
#define REACTION_h

#include <string>

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

namespace MitoD {

/**
* An abstract base class for all the reactions
*/
class Reaction {

public:	// Only constant parameters are public

	const szt			ind {};	 /**< index in Simulation::rc, i.e. index among all used and not used reactions */
	const real			rate {}; /**< reaction rate constant */

	// Convenience references
	const ulong&		it;		/**< ref: internal network iteration counter */
	const real&			time;	/**< ref: internal network time */

	const std::string	shortName;	/**< reaction name */
	const std::string	srt;		/**< reaction name abbreviation */

	/** Constructor */
	Reaction( Msgr& msgr,
			  const szt ind,
			  const real rate,
			  const ulong& it,		// const ref
			  const real& time,		// const ref
			  const std::string& shortName,
			  const std::string& srt
		)
		: ind {ind}
		, rate {rate}
		, it {it}
		, time {time}
		, shortName {shortName}
		, srt {srt}
		, msgr {msgr}
	{}

	/** virtual destructor */
	virtual ~Reaction() {};

	
	/** Sets the Gillespie score for this reaction */
	virtual void set_score() noexcept = 0;

	/** Returns the Gillespie score for this reaction */
	virtual real get_score() const noexcept = 0;

	/** Updates propensity for a pair of network components
	* @param c1 index of the 1st component to update
	* @param c2 index of the 2nd component to update
	*/
	virtual void update_prop(const szt c1, const szt c2) noexcept = 0;
	
	/** Executes the raction event */
	virtual void fire() noexcept = 0;

	/** Attaches this score to the Gillespie mechanism
	* @param a placeholder in the Gillespie object responsible for this reaction score
	*/
	virtual void attach_score_pointer(real* a) noexcept = 0;

	/** Populates the vector of reactions that need a score update after *this has fired
	* and initializes the propensities and effective rate
	* @param rc vector of unique pointers to all reactions taking part in the simulation
	*/
	virtual void initialize_dependencies(const vup<Reaction>& rc) noexcept = 0;

	/** The number of times this reaction was fired
	* @result the number of times this reaction was fired
	*/
	virtual szt event_count() const noexcept = 0;

	/** Prints the parameters common to all reactions
	* @param le true if new line after the output
	*/
	virtual	void print( const bool le ) const
	{
		msgr.print<false>(" it %d", it);
		msgr.print<false>(" srt "+srt);
		msgr.print<false>(" rate %f", rate);
		if (le) msgr.print("\n");
	}

protected:

	Msgr& msgr;		/**< ref: logging facility */

	/**
	* Pure virtual function: All network and reaction updates necessary after the given reaction event was executed
	*/
	virtual void update_netw_stats() = 0;

};

}	// namespace MitoD

#endif // REACTION_H
