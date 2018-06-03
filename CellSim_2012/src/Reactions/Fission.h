
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

template<typename Ntw>
class Fission: public Reaction
{
	friend Gillespie<Reaction>;

public:

	Fission( Oel& oel,
			 const szt ind,
			 Ntw& netw,
			 const real rate,
			 const szt& it,		// const ref
			 const real& time,	// const ref
			 const bool& verbose
		)
		: Reaction(oel, ind, rate, it, time, verbose, "fission", "fiss")
		, netw {netw}
		, rnd {netw.rnd}
	{
		netw.fis.verbose = verbose;
	}

	void set_score() noexcept override;
	real get_score() const noexcept override { return *score; };
	void update_prop(szt,szt) noexcept override;

	void make() noexcept override;
	szt event_count() const noexcept override { return eventCount; }

	static constexpr auto is_active(const std::unique_ptr<Reaction>&) noexcept;
	void initialize_dependencies(const vup<Reaction>&) noexcept override;

	void print(const bool) const override;

private:

	using Reaction::srt;
 	using Reaction::rate;
	using Reaction::oel;
	using Reaction::verbose;

	Ntw&					netw;
	RandFactory&			rnd;
	std::array<szt,2>					cc;
	real*					score {nullptr};
	szt						eventCount {0};
	std::vector<Reaction*>	dependents;	// reactions that need a score update after *this has fired

	void set_prop() noexcept;
	void attach_score_pointer(real*a) noexcept override { score = a; };
	void update_netw_stats() override;

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> constexpr
auto Fission<Ntw>::
is_active( const std::unique_ptr<Reaction>& r ) noexcept
{
	return  r->srt == "fiss";
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
	*score = rate * netw.fis.propTotal();
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
//	basic_update();						 was called from fiss
	
//	if (cc[0] != cc[1] && !( cc[0] == netw.clnum-1 || cc[1] == netw.clnum-1 ) ) 
//		Globals::eexit<string>( 1, 1, 1, logfile, "Error Fission::update_netw_stats: failed assumption for cluster-specific updates " );
//	if (cc[0] > cc[1] ) {
//		const auto tmp = cc[0];
//		cc[0] = cc[1];
//		cc[1] = tmp;
//	}
	
/*	for(auto& o : netw.clsi)						// clsi
		if(o.size() < netw.clnum)
			o.resize(netw.clnum);
	netw.populate_clsi(cc[0]);
	if (cc[0] != cc[1])
		netw.populate_clsi(cc[1]);

	if (netw.clcr.size() < netw.clnum) {			// clcr, clcrj
		netw.clcr.resize(netw.clnum);
		netw.clcrj.resize( netw.clnum);
	}

	netw.populate_clcr(cc[0]);
	if (cc[0] != cc[1])
		netw.populate_clcr( cc[1]);
*/
//	if (netw.fis.prop.size() < netw.clnum)			// fissVec
//		netw.fis.prop.resize(netw.clnum);

//	netw.update_OxPh();
//	netw.update_clq();
		
	for(auto& o : dependents) {
		o->update_prop(cc[0], cc[1]);
		o->set_score();
	}
}

template<typename Ntw>
void Fission<Ntw>::
make() noexcept
{
	if (verbose)
		print(true);

	eventCount++;

	cc = netw.fis.fire();

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
