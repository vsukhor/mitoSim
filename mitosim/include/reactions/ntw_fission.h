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

#ifndef NTW_FISSION_H
#define NTW_FISSION_H

namespace MitoD {

template<typename> class Fission;

/**
 * Base class for network-specific fission reaction slots.
 * @tparam Ntw type of the network
 */
template<typename Ntw>
class NtwFission {

public:

	friend Fission<Ntw>;

	explicit NtwFission(Ntw&);		/**< Constructor */

	/** sets this reaction propensity for the whole network */
	ulong set_prop()  noexcept;

	/**
	 * Updates this reaction propensity for the whole network.
	 * This is done after updating it for the cluster indexed.
	  * @param c cluster index that triggers the update
	*/
	void update_prop(const szt c) noexcept;

	/** prTotal getter */
	constexpr ulong get_prTotal() const noexcept { return prTotal; }

private:

	Ntw& host;		/**< ref: the host network for this reaction */

	// Convenience references to some of the host members
	typename Ntw::Reticulum&	mt;		/**< ref: the segments */
	const szt&					clnum;	/**< ref: current number fo clusters */

	// Propensities
	std::vector<ulong>	pr;				/**< propensities per cluster */
	ulong				prTotal {};		/**< total propensity */

	/**
	 * Sets this reaction propensity for the indexed cluster.
	 * Does not update the whole network propensity.
	 * @param c index of the cluster that is updateed
	 */
	void set_prop(const szt c) noexcept;

	/** Executes the reaction event */
	auto fire() noexcept;

	/** Find sa random node from those suitable for this reaction.
	 * @param w index of random segment
	 * @param a random position inside the segment 
	 */
	bool find_random_node(szt& w, szt& a) const noexcept;
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
auto NtwFission<Ntw>::
fire() noexcept
{
	szt w {huge<szt>};
	szt a {huge<szt>};

	find_random_node(w, a);

	return host.fiss(w, a);
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

}	// namespace MitoD

#endif // NTW_FISSION_H
