#ifndef _Fusion11_h
#define _Fusion11_h

#include "Fusion.h"

namespace MitoD {
	
template<typename Ntw>
class Fusion11: public Fusion<Ntw> {

	friend Gillespie<Reaction>;

public:

	Fusion11(
			Oel& oel,
			const szt ind,
			Ntw& netw,
			const real rate,
			const szt& it,		// const ref
			const real& time,	// const ref
			const bool verbose,
			const uint node_degree1,
			const uint node_degree2
		)
		: Fusion<Ntw>( oel, ind, netw, rate, it, time, verbose, node_degree1, node_degree2, "fu11" )
	{
		netw.fu11.verbose = verbose;
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
make() noexcept
{
//	oel.exit("ext");
//	std::cout <<"verbose_"<<verbose<<"_"<<std::endl;
	if(verbose)	print(true);

	eventCount++;

	cc = netw.fu11.fire();

	update_netw_stats();
}

template<typename Ntw>
void Fusion11<Ntw>::
print( const bool le ) const
{
	Fusion<Ntw>::print(false);
	oel.print0(" score "+STR(*score));
	oel.print0(" eventCount "+STR(eventCount));
	if (le) oel.print("\n");
}

}


#endif /* _Fusion11_h */
