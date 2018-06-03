
#ifndef NtwFusion12_h
#define NtwFusion12_h

#include "../FusionCandidates.h"

namespace MitoD {

template<typename> class Fusion12;

template<typename Ntw>
class NtwFusion12 {

public:

	friend Fusion12<Ntw>;

	explicit NtwFusion12(Ntw&);

	szt set_prop() noexcept;

private:

	Ntw&									host;
	RandFactory&							rnd;
	const typename Ntw::Reticulum&			mt;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;
	const std::vector<szt>&					mt22;
	const std::vector<szt>&					mt33;
	const szt&								minLoopLength;
	bool									verbose {false};

	FusionCandidatesXX	cnd;

	void populate() noexcept;
	std::array<szt,2> fire() noexcept;
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
	, minLoopLength {host.minLoopLength}
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
	cnd.clear();
	for (const auto w1 : mt11) 								// 11 ends to ...
		for (const auto e1 : {szt(1),szt(2)}) {
			const std::array<szt,2> we1 {w1,e1};
			for (const auto w2 : mt11)						// ... 11 bulk
				for (szt i=1; i<mt[w2].g.size(); i++) {
					const auto skip {w1 == w2 && ((e1 == 1 && i < minLoopLength) ||
												  (e1 == 2 && mt[w2].g.size()-i < minLoopLength))};
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
								 ((we1[1] == 1 && i < minLoopLength) ||
								  (we1[1] == 2 && mt[we2[0]].g.size()-i < minLoopLength))};
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
std::array<szt,2> NtwFusion12<Ntw>::
fire() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse12(cnd.u[r][0], cnd.u[r][1], cnd.v[r][0], cnd.v[r][1], verbose);
}

}

#endif /* NtwFusion12_h */