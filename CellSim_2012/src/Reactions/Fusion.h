#ifndef FUSION_h
#define FUSION_h

#include <string>
#include <vector>

#include "../utils/Oel.h"
#include "Reaction.h"

namespace MitoD {

template<typename Ntw>
class Fusion : public Reaction {

	friend Gillespie<Reaction>;

public:

	Fusion( Oel& oel,
			const szt ind,
			Ntw& netw,
			const real rate,
			const szt& it,			// const ref
			const real& time,		// const ref
			const bool verbose,
			const uint nodeDegree1,
			const uint nodeDegree2,
			const std::string& srt
		)
		: Reaction(oel, ind, rate, it, time, verbose, "fusion", srt)
		, netw {netw}
		, rnd {netw.rnd}
		, deg1(nodeDegree1)
		, deg2(nodeDegree2)
	{}

	static constexpr auto is_active( const std::unique_ptr<Reaction>& r );
	void initialize_dependencies( const vup<Reaction>& rc ) noexcept override;

	void print(const bool) const override;

protected:

	using Reaction::oel;

	Ntw&		 netw;
	RandFactory& rnd;
	std::array<szt,2>		 cc;

	void update_netw_stats() override;

private:

	using Reaction::set_score;
	using Reaction::it;

	std::vector<Reaction*>  dependents;	// reactions that need a score update after *this has fired
	uint					deg1, deg2;	// E: 2; real: 3

//	szt w_from_score( const real* const score ) const noexcept;
	virtual void set_prop() noexcept = 0;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw> constexpr
auto Fusion<Ntw>::
is_active( const std::unique_ptr<Reaction>& r )
{
	return  r->srt == "fu1L" ||
			r->srt == "fu11" ||
			r->srt == "fu12";
}

template<typename Ntw>
void Fusion<Ntw>::
initialize_dependencies( const vup<Reaction>& rc ) noexcept
{
	for (const auto& o : rc)
		if (
			Fusion<Ntw>::is_active(o) ||
			Fission<Ntw>::is_active(o)
			)
			dependents.push_back(o.get());

	set_prop();
	set_score();
}

template<typename Ntw> 
void Fusion<Ntw>::
update_netw_stats()
{
	netw.update_books();
	
	for (auto& o : dependents) {
		o->update_prop(cc[0], cc[1]);
		o->set_score();
	}
}

template<typename Ntw>
void Fusion<Ntw>::
print( const bool le ) const
{
	Reaction::print(false);
	oel.print<false>(" deg1 %d", deg1);
	oel.print<false>(" deg2 %d", deg2);
	if (le) oel.print("");
}

}

#endif /* FUSION_h */
