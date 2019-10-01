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

#ifndef NTW_FUSION11_H
#define NTW_FUSION11_H

#include "utils/common/misc.h"
#include "fusion_candidates.h"

namespace MitoD {

template<typename> class Fusion11;

/**
 * Network-specific reaction slot for fusion of two nodes of degree 1.
 * @tparam Ntw type of the network
 */
template<typename Ntw>
class NtwFusion11 {

public:

	friend Fusion11<Ntw>;

	explicit NtwFusion11(Ntw&);		/**< Constructor */

	/** sets this reaction propensity for the whole network */
	szt set_prop() noexcept;

private:

	Ntw& host;	/**< ref: the host network for this reaction */

	// Convenience references to some of the host members
	RandFactory&							rnd;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;

	FusionCandidatesXX	cnd; /**< node pairs suitable for this type of fusion */

	/** populates the vector of node pairs suitable for this type of fusion */
	void populate() noexcept;

	/** executes the reaction event */
	auto fire() noexcept;
};
// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion11<Ntw>::
NtwFusion11( Ntw& host )
	: host {host}
	, rnd {host.rnd}
	, mt11 {host.mt11}
	, mt13 {host.mt13}
{}

template<typename Ntw>
szt NtwFusion11<Ntw>::
set_prop() noexcept
{
	 populate();
	 return cnd.size();
}

template<typename Ntw>
void NtwFusion11<Ntw>::
populate() noexcept
{
	constexpr auto minLL = Structure<typename Ntw::ST>::minLoopLength;

	cnd.clear();
	const auto mtn11  {mt11.size()};
	for (szt i1=0; i1<mtn11; i1++) {					// 11 ends to ...
		const auto w1 {mt11[i1]};
		if (host.mt[w1].g.size() >= minLL)				// ... opposite end in the same segment
			cnd.add({w1,1}, {w1,2});

		for (const auto e1 : {szt(1),szt(2)}) {
			for (szt i2=i1+1; i2<mtn11; i2++)			// ... other 11 segments (both ends to both ens)
				for (const auto e2 : {szt(1),szt(2)})
					cnd.add({w1,e1}, {mt11[i2],e2});

			for (const auto& we2 : mt13)				// ... free ends of 13
				cnd.add({w1, e1}, we2);

			for (const auto& we2 : mt13)				// ... free ends of 14
				cnd.add({w1, e1}, we2);
		}
	}
	const auto mtn13  {mt13.size()};
	for (szt i1=0; i1<mtn13; i1++) 						// free ends of 13 to ...
		for (szt i2=i1+1; i2<mtn13; i2++)				// ... free ends of other 13
			cnd.add(mt13[i1], mt13[i2]);
}

template<typename Ntw>
auto NtwFusion11<Ntw>::
fire() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse11(cnd.u[r][0], cnd.u[r][1], cnd.v[r][0], cnd.v[r][1]);
}

}	// namespace MitoD

#endif // NTW_FUSION11_H
