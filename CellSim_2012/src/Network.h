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

#ifndef Network_h
#define Network_h

#include "AbilityForFusion.h"
#include "Reactions/NtwFission.h"
#include "Reactions/NtwFusion11.h"
#include "Reactions/NtwFusion12.h"
#include "Reactions/NtwFusion1L.h"
#include "Simulation.h"

namespace MitoD {

/**
 * The Network class.
 * Represents a fully dynamic network.
 */
template<typename SegmentT>
class Network
	:  public AbilityForFusion<SegmentT> {

public:
	
	typedef Network<SegmentT> thisT;
	typedef SegmentT ST;

	using Structure<SegmentT>::mt;
	using Structure<SegmentT>::mtnum;
	using Structure<SegmentT>::clnum;
	using Structure<SegmentT>::mtmass;
	using Structure<SegmentT>::nn;
	using Structure<SegmentT>::glm;
	using Structure<SegmentT>::oel;

	const Config& 			cfg;	/**< configuration */
	RandFactory& 			rnd;	/**< random number factory */
	real					time;	/**< current time */
	szt						it;		/**< iteration counter */

	NtwFission<thisT>		fis;	/**< slot for fission reaction */
	NtwFusion11<thisT>		fu11;	/**< slot for fusion of raction of nodes with degrees 1 and 1 */
	NtwFusion12<thisT>		fu12;	/**< slot for fusion of raction of nodes with degrees 1 and 2 */
	NtwFusion1L<thisT>		fu1L;	/**< slot for fusion of raction of nodes with degrees 1 and a loop */

	/** Constructor.
	 * @param cfg configuration object
	 * @param runIndex run index
	 * @param rnd random number factory
	 * @param oel logging facility object
	 */
	explicit Network(
			const Config& cfg,
			szt runIndex,
			RandFactory& rnd,
			Oel& oel
		);

	/** Generates the network components */
	void generate_mitos();

	/** update the network state variables */
	void update_books() noexcept;

	/** Writes network to a file.
	 * @param startnew start a new file vs. adding new data records
	 * @param last is the final writeout
	 * @param itr current simulation iteration
	 * @param t current simulation time
	 */
	void save_mitos(const bool startnew, const bool last, const szt itr, const real t) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename SegmentT>
Network<SegmentT>::
Network(
		const Config& cfg,
		szt runIndex,
		RandFactory& rnd,
		Oel& oel
	)
	: AbilityForFusion<SegmentT> {cfg, oel}
	, cfg {cfg}
	, rnd {rnd}
	, time {zero<real>}
	, it {}
	, fis {*this}
	, fu11 {*this}
	, fu12 {*this}
	, fu1L {*this}
{
	generate_mitos();
	update_books();
	Simulation<thisT> sim {*this, rnd, time, it, oel};
	sim();
}

template<typename SegmentT>
void Network<SegmentT>::
generate_mitos()
{
	mtnum = cfg.mtmassini / cfg.segmassini;
	if (mtnum < 1)
		oel.exit("The system should have at least one segment initially");

	mtmass = 0;
	szt m {1};
	szt ei {};
	mt.emplace_back(oel);			// an "empty" one
	while (m <= mtnum) {
		mt.emplace_back(cfg.segmassini, m-1, mtmass, ei, oel);
		m++;
	}
	clnum = mtnum;
	oel.print("Generated mtnum %d of mtmass: %d", mtnum, mtmass);
}

template<typename SegmentT>
void Network<SegmentT>::
update_books() noexcept
{
	this->basic_update();
}

template<typename SegmentT>
void Network<SegmentT>::
save_mitos( const bool startnew, const bool last, const szt itr, const real t ) const
{
	const auto fname = (last) ? cfg.workingDirOut+"mitos_last_"+cfg.runName
							  : cfg.workingDirOut+"mitos_"	   +cfg.runName;
	const auto flags = (startnew) ? std::ios::binary | std::ios::trunc
								  : std::ios::binary | std::ios::app;
	std::ofstream ofs {fname, flags};
	if (ofs.fail())
		oel.print("Cannot open file: "+fname);

	ofs.write((char*) &t, sizeof(real));						//	cout << t << " ";
	ofs.write((char*) &mtnum, sizeof(szt));						//	cout << mtnum << endl;

	static szt mtnummax, nn1max, nn2max;
	if (!last) {
		if (startnew) {
			mtnummax = 0;
			nn1max = 0;
			nn2max = 0;
		}
		if (mtnum > mtnummax)
			mtnummax = mtnum;
	}
	for (szt q=1; q<=mtnum; q++) {
		mt[q].write( ofs );										//	cout << endl;
		if (!last) {
			if (mt[q].nn[1] > nn1max) nn1max = mt[q].nn[1];
			if (mt[q].nn[2] > nn2max) nn2max = mt[q].nn[2];
		}
	}
	ofs.write((char*) &mtnummax, sizeof(szt));					//	cout << endl << mtnummax << " ";
	ofs.write((char*) &nn1max, sizeof(szt));					//	cout << nn1max << " ";
	ofs.write((char*) &nn2max, sizeof(szt));					//	cout << nn2max << " ";

	szt nst2save;
	nst2save = last
			 ? szt{}
			 : szt(itr/cfg.saveFrequency);
	ofs.write((char*) &nst2save, sizeof(szt));				//	cout << nst2save << " ";
}

}

#endif /* Network_h */
	
	
	
	
	
	
	
	
