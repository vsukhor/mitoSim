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

#ifndef NtwFusion1L_h
#define NtwFusion1L_h

#include "../FusionCandidates.h"

namespace MitoD {

template<typename> class Fusion1L;

/**
 * Network-specific reaction slot for fusion of a degree 1 node with a looped segment.
 */
template<typename Ntw>
class NtwFusion1L {

	friend Fusion1L<Ntw>;

public:

	explicit NtwFusion1L(Ntw&);		/**< Constructor */

	/** sets this reaction propensity for the whole network */
	szt set_prop() noexcept;

private:

	Ntw& host;	/**< ref: the host network for this reaction */

	// Convenience references to some of the host members
	RandFactory&							rnd;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;
	const std::vector<szt>&					mt22;

	bool verbose {};			/**< verbosity of the short logs */

	FusionCandidatesXL	cnd;	/**< node pairs suitable for this type of fusion */

	/** populates the vector of node pairs suitable for this type of fusion */
	void populate() noexcept;

	/** executes the reaction event */
	auto operator()() noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion1L<Ntw>::
NtwFusion1L( Ntw& host )
	: host {host}
	, rnd {host.rnd}
	, mt11 {host.mt11}
	, mt13 {host.mt13}
	, mt22 {host.mt22}
{}

template<typename Ntw>
szt NtwFusion1L<Ntw>::
set_prop() noexcept
{
	 populate();
	 return cnd.size();
}

template<typename Ntw>
void NtwFusion1L<Ntw>::
populate() noexcept
{
	cnd.clear();
	for (const auto w2 : mt22) {
		for (const auto w1 : mt11)
			for (const auto e1 : {szt(1),szt(2)})
				cnd.add({w1,e1}, w2);		// e2 is 1 by convention

		for (const auto& we1 : mt13)
			cnd.add(we1, w2);
	}
}

template<typename Ntw>
auto NtwFusion1L<Ntw>::
operator()() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse1L(cnd.u[r][0], cnd.u[r][1], cnd.v[r], verbose);
}


}

#endif /* NtwFusion_h */
