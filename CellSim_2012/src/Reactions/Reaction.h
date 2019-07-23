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
#include "../utils/Misc.h"

namespace MitoD {

/**
* An abstract base class for all the reactions
*/

class Reaction {

public:	// Only constant parameters are public

	const szt			ind {};	 /**< index in Simulation::rc, i.e. index among all used and not used reactions */
	const real			rate {}; /** reaction rate constant */

	// Convenience references
	const szt&			it;		/**< ref: internal network iteration counter */
	const real&			time;	/**< ref: internal network time */

	const bool			verbose;	/**< verbosity of the short logs */
	const std::string	srtype;		/**< reaction name */
	const std::string	srt;		/**< reaction name abbreviation */

	/** Constructor & destructor */
	Reaction( Oel& oel,
			  const szt ind,
			  const real rate,
			  const szt& it,		// const ref
			  const real& time,		// const ref
			  const bool verbose,
			  const std::string& srtype,
			  const std::string& srt
		)
		: ind {ind}
		, rate {rate}
		, it {it}
		, time {time}
		, verbose {verbose}
		, srtype {srtype}
		, srt {srt}
		, oel {oel}
	{}
	virtual ~Reaction() {};

	
	/** Sets the Gillespie score for this reaction */
	virtual void set_score() noexcept = 0;

	/** Returns the Gillespie score for this reaction */
	virtual real get_score() const noexcept = 0;

	/** Updates propensity for two network components */
	virtual void update_prop(szt,szt) noexcept = 0;
	
	/** Executes the raction event */
	virtual void operator()() noexcept = 0;

	/** Attaches this score to the Gillespie mechanism */
	virtual void attach_score_pointer(real*) noexcept = 0;

	/**
	* Populates the vector of reactions that need a score update after *this has fired
	* and initializes the propensities and effective rate
	*/
	virtual void initialize_dependencies(const vup<Reaction>&) noexcept = 0;

	/** Returns the number of times this reaction was fired */
	virtual szt get_eventCount() const noexcept = 0;

	virtual	void print( const bool le ) const
	{
	/** Prints out parameters common to all reactions */
		oel.print<false>(" it %d", it);
		oel.print<false>(" srt "+srt);
		oel.print<false>(" rate %f", rate);
		if (le) oel.print("\n");
	}

protected:

	Oel& oel;		/**< ref: logging facility */

	/**
	* All network and reaction updates necessary after the given reaction event was executed
	*/
	virtual void update_netw_stats() = 0;

};

}

#endif /* REACTION_h */
