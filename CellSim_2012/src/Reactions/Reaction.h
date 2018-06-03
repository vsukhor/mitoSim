#ifndef REACTION_h
#define REACTION_h

#include <string>
#include <memory>
#include "../utils/Misc.h"

namespace MitoD {

class Reaction {

public:

	const szt			ind;				// index in Simulation::rc, i.e. the index among all used and not used reactions
	const real			rate {zero<real>};
	const szt&			it;
	const real&			time;
	const bool			verbose;
	const std::string	srtype;
	const std::string	srt;

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
		, rate{rate}
		, it {it}
		, time {time}
		, verbose {verbose}
		, srtype {srtype}
		, srt {srt}
		, oel {oel}
	{}

	virtual void set_score() noexcept = 0;
	virtual real get_score() const noexcept = 0;
	virtual void update_prop(szt,szt) noexcept = 0;
	virtual void make() noexcept = 0;
	virtual void update_netw_stats() = 0;
	virtual void attach_score_pointer(real*) noexcept = 0;
	virtual void initialize_dependencies( const vup<Reaction>& ) noexcept = 0;
	virtual szt event_count() const noexcept = 0;

	virtual	void print( bool le ) const
	{
		oel.print<false>(" it %d", it);
		oel.print<false>(" srt "+srt);
		oel.print<false>(" rate %f", rate);
		if (le) oel.print("\n");
	}

protected:

	Oel	&oel;
};

}

#endif /* REACTION_h */
