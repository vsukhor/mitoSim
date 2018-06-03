#ifndef _Fusion1L_h
#define _Fusion1L_h

#include "Fusion.h"

namespace MitoD {
	
template<typename Ntw>
class Fusion1L: public Fusion<Ntw> {

	friend Gillespie<Reaction>;

public:

	Fusion1L( Oel& oel,
			  const szt ind,
			  Ntw& netw,
			  const real rate,
			  const szt& it,		// const ref
			  const real& time,		// const ref
			  const bool verbose,
			  const uint nodeDegree1,
			  const uint nodeDegree2 )
		: Fusion<Ntw>(oel, ind, netw, rate, it, time, verbose, nodeDegree1, nodeDegree2, "fu1L")
	{
		netw.fu1L.verbose = verbose;
	}

	void set_score() noexcept final;
	real get_score() const noexcept final { return *score; };
	void update_prop(szt,szt) noexcept final;

	void make() noexcept final;
	szt event_count() const noexcept final { return eventCount; }

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

	real*	score {nullptr};
	szt		eventCount {0};
	szt		propTotal {0};

	void set_prop() noexcept final;
	void attach_score_pointer(real* a) noexcept final { score = a; };
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
make() noexcept
{
	if(verbose)	print(true);

	eventCount++;
	
	cc = netw.fu1L.fire();

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
