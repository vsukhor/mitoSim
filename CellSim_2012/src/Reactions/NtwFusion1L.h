
#ifndef NtwFusion1L_h
#define NtwFusion1L_h

#include "../FusionCandidates.h"

namespace MitoD {

template<typename> class Fusion1L;

template<typename Ntw>
class NtwFusion1L {

	friend Fusion1L<Ntw>;

public:

	explicit NtwFusion1L(Ntw&);

	szt set_prop() noexcept;

private:

	Ntw&									host;
	RandFactory&							rnd;
	const std::vector<szt>&					mt11;
	const std::vector<std::array<szt,2>>&	mt13;
	const std::vector<szt>&					mt22;
	bool									verbose {false};

	FusionCandidatesXL						cnd;

	void populate() noexcept;
	std::array<szt,2> fire() noexcept;
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
std::array<szt,2> NtwFusion1L<Ntw>::
fire() noexcept
{
	const auto r {rnd.uniform0(cnd.size())};

	return host.fuse1L(cnd.u[r][0], cnd.u[r][1], cnd.v[r], verbose);
}


}

#endif /* NtwFusion_h */
