#ifndef Network_h
#define Network_h

#include "NetworkCore.h"
#include "Netw_IO.h"
#include "Reactions/NtwFission.h"
#include "Reactions/NtwFusion11.h"
#include "Reactions/NtwFusion12.h"
#include "Reactions/NtwFusion1L.h"

namespace MitoD {

template<typename SegmentT>
class Network
	: public NetworkCore<SegmentT> {

public:
	
	typedef Network<SegmentT> thisT;
	typedef SegmentT ST;

	using Structure<SegmentT>::clnum;

	Netw_IO<false,thisT>	io;
	NtwFission<thisT>		fis;
	NtwFusion11<thisT>		fu11;
	NtwFusion12<thisT>		fu12;
	NtwFusion1L<thisT>		fu1L;

	explicit Network(
			const MitoD::ConfigMain& cfgMain,
			const std::string& runName,
			RandFactory& rnd,
			const real& time,	// const ref
			const szt& it,		// const ref
			Oel& oel
		);

	void update_books() noexcept;

private:

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename SegmentT>
Network<SegmentT>::
Network(
		const MitoD::ConfigMain& cfgMain,
		const std::string& runName,
		RandFactory& rnd,
		const real& time,	// const ref
		const szt& it,		// const ref
		Oel& oel
	)
	: NetworkCore<SegmentT> {cfgMain, runName, rnd, time, it, oel}
	, io {*this}
	, fis {*this}
	, fu11 {*this}
	, fu12 {*this}
	, fu1L {*this}
{
	this->generate_mitos();
		update_books();
}

template<typename SegmentT>
void Network<SegmentT>::
update_books() noexcept			// modified
{
	this->basic_update();
//	this->populate_clcr();			 modified
//	this->populate_clsi();			modified
}

}

#endif /* Network_h */
	
	
	
	
	
	
	
	
