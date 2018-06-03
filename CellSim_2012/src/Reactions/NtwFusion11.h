
#ifndef NtWFusion11_h
#define NtWFusion11_h

#include "../FusionCandidates.h"

namespace MitoD {

template<typename> class Fusion11;

template<typename Ntw>
class NtwFusion11 {

public:

	friend Fusion11<Ntw>;

	explicit NtwFusion11(Ntw&);

	szt set_prop() noexcept;

private:

	Ntw&									host;
	RandFactory&							rnd;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;
	const szt&								minLoopLength;
	bool									verbose {false};

	FusionCandidatesXX						cnd;

	void populate() noexcept;
	std::array<szt,2> fire() noexcept;
};
// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Ntw>
NtwFusion11<Ntw>::
NtwFusion11( Ntw& host )
	: host {host}
	, rnd {host.rnd}
	, mt11 {host.mt11}
	, mt13 {host.mt13}
	, minLoopLength {host.minLoopLength}
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
	cnd.clear();
	const auto mtn11  {mt11.size()};
	for (szt i1=0; i1<mtn11; i1++) {					// 11 ends to ...
		const auto w1 {mt11[i1]};
		if (host.mt[w1].g.size() >= minLoopLength)			// ... opposite end in the same mito
			cnd.add({w1,1}, {w1,2});

		for (const auto e1 : {szt(1),szt(2)}) {
			for (szt i2=i1+1; i2<mtn11; i2++)			// ... other 11 mitos (both ends to both ens)
				for (const auto e2 : {szt(1),szt(2)})
					cnd.add({w1,e1}, {mt11[i2],e2});

			for (const auto& we2 : mt13)				// ... free ends of 13
				cnd.add({w1, e1}, we2);

			for (const auto& we2 : mt13)				// ... free ends of 14
				cnd.add({w1, e1}, we2);
		}
	}
	const auto mtn13  {mt13.size()};
	for (szt i1=0; i1<mtn13; i1++) {					// free ends of 13 to ...
		for (szt i2=i1+1; i2<mtn13; i2++)				// ... free ends of other 13
			cnd.add(mt13[i1], mt13[i2]);
	}
}

template<typename Ntw>
std::array<szt,2> NtwFusion11<Ntw>::
fire() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse11(cnd.u[r][0], cnd.u[r][1], cnd.v[r][0], cnd.v[r][1], verbose);
}

}

#endif /* Fusion11_h */
