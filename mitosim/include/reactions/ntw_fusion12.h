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

#ifndef NTW_FUSION12_H
#define NTW_FUSION12_H

#include "../fusion_candidates.h"

namespace MitoD {

template<typename> class Fusion12;

/**
 * Network-specific reaction slot for fusion of a degree 1 node with a degree 2 node.
 * @tparam Ntw type of the network
 */
template<typename Ntw>
class NtwFusion12 {

public:

	friend Fusion12<Ntw>;

	explicit NtwFusion12(Ntw&);		/**< Constructor */

	/** sets this reaction propensity for the whole network */
	szt set_prop() noexcept;

private:

	Ntw& host;	/**< ref: the host network for this reaction */
	
	// Convenience references to some of the host members
	RandFactory&							rnd;
	const typename Ntw::Reticulum&			mt;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;
	const std::vector<szt>&					mt22;
	const std::vector<szt>&					mt33;

	FusionCandidatesXX	cnd;	/**< node pairs suitable for this type of fusion */

	/** populates the vector of node pairs suitable for this type of fusion */
	void populate() noexcept;

	/** executes the reaction event */
	auto fire() noexcept;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion12<Ntw>::
NtwFusion12( Ntw& host )
	: host {host}
	, rnd {host.rnd}
	, mt {host.mt}
	, mt11 {host.mt11}
	, mt13 {host.mt13}
	, mt22 {host.mt22}
	, mt33 {host.mt33}
{}

template<typename Ntw>
szt NtwFusion12<Ntw>::
set_prop() noexcept
{
	 populate();
	 return cnd.size();
}

template<typename Ntw>
void NtwFusion12<Ntw>::
populate() noexcept
{
	constexpr auto minLL = Structure<typename Ntw::ST>::minLoopLength;
	cnd.clear();
	for (const auto w1 : mt11) 								// 11 ends to ...
		for (const auto e1 : {szt(1),szt(2)}) {
			const std::array<szt,2> we1 {w1,e1};
			for (const auto w2 : mt11)						// ... 11 bulk
				for (szt i=1; i<mt[w2].g.size(); i++) {
					const auto skip {w1 == w2 && ((e1 == 1 && i < minLL) ||
												  (e1 == 2 && mt[w2].g.size()-i < minLL))};
					if (!skip) {
						cnd.add(we1, {w2,i});
						cnd.add(we1, {w2,i});
					}
				}
			for (const auto& we2 : mt13) 					// ... 13 bulk
				for (szt i=1; i<mt[we2[0]].g.size(); i++)
					cnd.add(we1, {we2[0],i});

			for (const auto w2 : mt33)						// ... 33 bulk
				for (szt i=1; i<mt[w2].g.size(); i++)
					cnd.add(we1, {w2,i});

			for (const auto w2 : mt22)						// ... 22 bulk
				for (szt i=1; i<mt[w2].g.size(); i++)
					cnd.add(we1, {w2,i});
		}

	for (const auto& we1 : mt13) {							// a free end of 13 to ...
		for (const auto w2 : mt11)							// ... 11 bulk
			for (szt i=1; i<mt[w2].g.size(); i++)
				cnd.add(we1, {w2,i});

		for (const auto& we2 : mt13) {						// ... 13 bulk
			for (szt i=1; i<mt[we2[0]].g.size(); i++) {
				const auto skip {we1[0] == we2[0] &&
								 ((we1[1] == 1 && i < minLL) ||
								  (we1[1] == 2 && mt[we2[0]].g.size()-i < minLL))};
				if (!skip)
					cnd.add(we1, {we2[0],i});
			}
		}
		for (const auto w2 : mt33)							// ... 33 bulk
			for (szt i=1; i<mt[w2].g.size(); i++)
				cnd.add(we1, {w2,i});

		for (const auto w2 : mt22)							// ... 22 bulk
			for (szt i=1; i<mt[w2].g.size(); i++)
				cnd.add(we1, {w2,i});
	}
}

template<typename Ntw>
auto NtwFusion12<Ntw>::
fire() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse12(cnd.u[r][0], cnd.u[r][1], cnd.v[r][0], cnd.v[r][1]);
}

}	// namespace MitoD

#endif // NTW_FUSION12_H
