
#ifndef NtwFission_h
#define NtwFission_h

namespace MitoD {

template<typename> class Fission;

template<typename Ntw>
class NtwFission {

public:

	friend Fission<Ntw>;

	explicit NtwFission(Ntw&);

	ulong set_prop()  noexcept;
	void update_prop(const szt) noexcept;
	constexpr ulong propTotal() const noexcept { return prTotal; }
	std::array<szt,2> fire() noexcept;

private:

	Ntw&						host;
	typename Ntw::Reticulum&	mt;
	const szt&					clnum;
	bool						verbose {false};

	std::vector<ulong>	pr;
	ulong				prTotal {0UL};

	void set_prop(const szt)  noexcept;
	void updateVecSmall() noexcept;

	bool find_random_node(szt&, szt&) const noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFission<Ntw>::
NtwFission( Ntw& host )
	: host {host}
	, mt {host.mt}
	, clnum {host.clnum}
{}

template<typename Ntw>
ulong NtwFission<Ntw>::
set_prop() noexcept
{
	pr.resize(clnum);
	for (szt ic=0; ic<clnum; ic++)
		set_prop(ic);

	return std::accumulate(pr.begin(), pr.end(), 0UL);
}

template<typename Ntw>
void NtwFission<Ntw>::
set_prop( const szt ic ) noexcept
{
	pr[ic] = 0UL;
	for (const auto w : host.clmt[ic]) {
		pr[ic] += mt[w].set_end_fin(1) +
				  mt[w].set_end_fin(2);
		for (szt a=0; a<mt[w].g.size()-1; a++)
			pr[ic] += 2UL * mt[w].set_bulk_fin(a);
	}
}

template<typename Ntw>
void NtwFission<Ntw>::
update_prop( const szt c ) noexcept		// incremental clnum changes are assumed
{
	if (pr.size() > clnum) {
		pr.resize(clnum);
	}
	else if (pr.size() < clnum)
		pr.resize(clnum);

	if (c < clnum)
		set_prop(c);

	prTotal = std::accumulate(pr.begin(), pr.end(), zero<real>);
}

template<typename Ntw>
std::array<szt,2> NtwFission<Ntw>::
fire() noexcept
{
	szt w {huge<szt>};
	szt a {huge<szt>};

	find_random_node(w, a);

	return host.fiss(w, a, verbose);
}

template<typename Ntw>
bool NtwFission<Ntw>::
find_random_node( szt& w, szt& a ) const noexcept
{
	auto k {host.rnd.uniform1(prTotal)};
	auto ksum {0UL};
	for (w=1; w<=host.mtnum; w++) {
		const auto& g {mt[w].g};
		a = 0;
		ksum += g[a].fin[0];
		if (k <= ksum)
			return true;
		for (; a<g.size()-1;) {
			ksum += g[a].fin[1];
			a++;
			ksum += g[a].fin[0];
			if (k <= ksum)
				return true;
		}
		ksum += g[a].fin[1];
		if (k <= ksum) {
			a++;
			return true;
		}
	}
	return false;
}

}

#endif /* NtwFission_h */
