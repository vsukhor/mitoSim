#ifndef _Fusion12_h
#define _Fusion12_h

#include "Fusion.h"

namespace MitoD {
	
template<typename Ntw>
class Fusion12: public Fusion<Ntw>
{
	friend Gillespie<Reaction>;

public:

	Fusion12(
			Oel& oel,
			const szt ind,
			Ntw& netw,
			const real rate,
			const szt& it,			// const ref
			const real& time,		// const ref
			const bool verbose,
			const uint node_degree1,
			const uint node_degree2
		)
		: Fusion<Ntw>( oel, ind, netw, rate, it, time, verbose, node_degree1, node_degree2, "fu12" )
	{
		netw.fu12.verbose = verbose;
	}

	void set_score() noexcept final;
	real get_score() const noexcept final { return *score; };
	void update_prop(szt,szt) noexcept final;

	void make() noexcept override;
	szt event_count() const noexcept final { return eventCount; }

	void print(const bool) const final;

private:

	using Reaction::rate;
	using Reaction::verbose;
	using Reaction::oel;
	using Fusion<Ntw>::netw;
	using Fusion<Ntw>::cc;
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
void Fusion12<Ntw>::
set_score() noexcept
{
	*score = rate * propTotal;
}

template<typename Ntw>
void Fusion12<Ntw>::
set_prop() noexcept
{
	propTotal = netw.fu12.set_prop();
}

template<typename Ntw>
void Fusion12<Ntw>::
update_prop( szt, szt ) noexcept
{
	set_prop();
}

template<typename Ntw>
void Fusion12<Ntw>::
make() noexcept
{
	if(verbose) print(true);

	eventCount++;

	cc = netw.fu12.fire();

	update_netw_stats();
}

template<typename Ntw>
void Fusion12<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	oel.print0(" score "+STR(*score));
	oel.print0(" eventCount "+STR(eventCount));
	if (le) oel.print("\n");
}

}


#endif	/* _Fusion12_h */
