/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann,
   Helmholtz Center for Infection Research (Braunschweig, Germany).
   All Rights Reserved.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

============================================================================== */

/**
* \file network.h
* \brief High-level network components and functionality.
* \author Valerii Sukhorukov
*/

#ifndef NETWORK_H
#define NETWORK_H

#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "config.h"
#include "ability_for_fusion.h"
#include "ntw_fission.h"
#include "ntw_fusion11.h"
#include "ntw_fusion12.h"
#include "ntw_fusion1u.h"
#include "simulation.h"

namespace MitoD {

/**
 * The Network class template.
 * Represents a fully dynamic network, capable for both fusion and division.
 * @tparam SegmentT type of the segment used by the network
 */
template<typename SegmentT>
class Network
	:  public AbilityForFusion<SegmentT> {

public:
	
	using thisT = Network<SegmentT>;
	using ST = SegmentT;

	using Structure<SegmentT>::mt;
	using Structure<SegmentT>::mtnum;
	using Structure<SegmentT>::clnum;
	using Structure<SegmentT>::mtmass;
	using Structure<SegmentT>::nn;
	using Structure<SegmentT>::glm;
	using Structure<SegmentT>::msgr;

	const Config& 			cfg;	/**< configuration */
	RandFactory& 			rnd;	/**< random number factory */
	double					time;	/**< current time */
	ulong					it;		/**< iteration counter */

	NtwFission<thisT>		fis;	/**< slot for fission reaction */
	NtwFusion11<thisT>		fu11;	/**< slot for fusion of raction of nodes with degrees 1 and 1 */
	NtwFusion12<thisT>		fu12;	/**< slot for fusion of raction of nodes with degrees 1 and 2 */
	NtwFusion1L<thisT>		fu1L;	/**< slot for fusion of raction of nodes with degrees 1 and a loop */

	/** Constructor.
	 * @param cfg configuration object
	 * @param runIndex run index
	 * @param rnd random number factory
	 * @param msgr logging facility object
	 */
	explicit Network(
			const Config& cfg,
			szt runIndex,
			RandFactory& rnd,
			Msgr& msgr
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
	void save_mitos(const bool startnew,
					const bool last,
					const szt itr,
					const real t) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename SegmentT>
Network<SegmentT>::
Network(
		const Config& cfg,
		szt runIndex,
		RandFactory& rnd,
		Msgr& msgr
	)
	: AbilityForFusion<SegmentT> {msgr}
	, cfg {cfg}
	, rnd {rnd}
	, time {zero<double>}
	, it {}
	, fis {*this}
	, fu11 {*this}
	, fu12 {*this}
	, fu1L {*this}
{
	generate_mitos();
	update_books();
	Simulation<thisT> sim {*this, rnd, time, it, msgr};
	sim();
}

template<typename SegmentT>
void Network<SegmentT>::
generate_mitos()
{
	mtnum = cfg.mtmassini / cfg.segmassini;
	if (mtnum < 1)
		msgr.exit("The system should have at least one segment initially");

	mtmass = 0;
	szt m {1};
	szt ei {};
	mt.emplace_back(msgr);			// an "empty" one
	while (m <= mtnum) {
		mt.emplace_back(cfg.segmassini, m-1, mtmass, ei, msgr);
		m++;
	}
	clnum = mtnum;
	msgr.print("Generated mtnum %d of mtmass: %d", mtnum, mtmass);
}

template<typename SegmentT>
void Network<SegmentT>::
update_books() noexcept
{
	this->basic_update();
}

template<typename SegmentT>
void Network<SegmentT>::
save_mitos( const bool startnew,
			const bool last,
			const szt itr,
			const real t ) const
{
	const auto fname = (last) ? cfg.workingDirOut+"mitos_last_"+cfg.runName
							  : cfg.workingDirOut+"mitos_"	   +cfg.runName;
	const auto flags = (startnew) ? std::ios::binary | std::ios::trunc
								  : std::ios::binary | std::ios::app;
	std::ofstream ofs {fname, flags};
	if (ofs.fail())
		msgr.print("Cannot open file: "+fname);

	ofs.write((char*) &t, sizeof(real));
	ofs.write((char*) &mtnum, sizeof(szt));

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
		mt[q].write( ofs );
		if (!last) {
			if (mt[q].nn[1] > nn1max) nn1max = mt[q].nn[1];
			if (mt[q].nn[2] > nn2max) nn2max = mt[q].nn[2];
		}
	}
	ofs.write(reinterpret_cast<const char*>(&mtnummax), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&nn1max), sizeof(szt));
	ofs.write(reinterpret_cast<const char*>(&nn2max), sizeof(szt));

	szt nst2save;
	nst2save = last
			 ? szt{}
			 : szt(itr/cfg.saveFrequency);
	ofs.write(reinterpret_cast<const char*>(&nst2save), sizeof(szt));
}

}	// namespace MitoD

#endif // NETWORK_H
	
	
	
	
	
	
	
	