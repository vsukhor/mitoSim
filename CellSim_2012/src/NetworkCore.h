#ifndef NetworkCore_h
#define NetworkCore_h

#include "ConfigMito.h"

#include "AbilityForFusion.h"

namespace MitoD {

template<typename> class Structure;
template<typename> class CoreTransformer;

template<typename Mt>
class NetworkCore : public ConfigMito,
					public AbilityForFusion<Mt>
{
public:

	using Structure<Mt>::mt;
	using Structure<Mt>::mtnum;
	using Structure<Mt>::clnum;
	using Structure<Mt>::mtmass;
	using Structure<Mt>::nn;
	using Structure<Mt>::glm;
	using Structure<Mt>::oel;
	using CoreTransformer<Mt>::rename_mito;
	using CoreTransformer<Mt>::copy_neibs;


	const MitoD::ConfigMain&		cfgMain;
	const std::string				runName;
	RandFactory&					rnd;

	std::vector<std::pair<szt,szt>>	m2shorten;

	const real&						time;
	const szt&						it;

	explicit NetworkCore(
		const MitoD::ConfigMain& cfgMain,
		const std::string& runName,
		RandFactory& rnd,
		const real& time,
		const szt& it,
		Oel& oel
	);

	void generate_mitos();
	void shift_end_edge( const szt f, const szt ef, const szt t, const szt et );
	void edgeInd2mitoLink( const szt ind1, const szt ind2, szt& m1, szt& e1, szt& m2, szt& e2 ) const noexcept;
	void find_edge( szt& m, szt& a, const szt& ind ) const;
	void sort_mitos_length( szt* sorted ) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
NetworkCore<Mt>::
NetworkCore(
		const MitoD::ConfigMain& cfgMain,
		const std::string& runName,
		RandFactory& rnd,
		const real& time,
		const szt& it,
		Oel& oel
	)
	: ConfigMito {cfgMain, oel}
	, AbilityForFusion<Mt> {cfgMain, oel}
	, cfgMain {cfgMain}
	, runName {runName}
	, rnd {rnd}
	, time {time}
	, it {it}
{
	Mt::ThisEdge::length = this->edgeLength;
}


template<typename Mt>
void NetworkCore<Mt>::
generate_mitos()
{
	mtnum = this->mtmassini / this->segmassini;
	if (mtnum < 1)
		oel.exit("The system should have at least one mito initially"); 

	mtmass = 0;
	szt m {1};
	szt ei {0};
	mt.emplace_back(0, 0, mtmass, ei, oel);				// "empty" one
	while (m <= mtnum) {
		mt.emplace_back(this->segmassini, m-1, mtmass, ei, oel);
		m++;
	}
	clnum = mtnum;
	oel.print("Generated mtnum %d of mtmass: %d", mtnum, mtmass);
}


// shifts an edge from ef end of f to et end of t
template<typename Mt>
void NetworkCore<Mt>::
shift_end_edge( const szt f, const szt ef, const szt t, const szt et )
{
//	print_mitos( " SEE0 " );
	if (ef == 2) {
		if (et == 1) mt[t].g.insert(mt[t].g.begin(), mt[f].g.back());
		else 		 mt[t].g.push_back(mt[f].g.back());
	}
	else {						// ef == 1
		if (et == 1) mt[t].g.insert(mt[t].g.begin(), mt[f].g[0]);
		else 		 mt[t].g.push_back(mt[f].g[0]);

		for (szt i=1; i<mt[f].g.size(); i++)
			mt[f].g[i-1] = mt[f].g[i];
	}
	mt[f].g.pop_back();
//	print_mitos( " SEE1 " );
}

template<typename Mt>
void NetworkCore<Mt>::
edgeInd2mitoLink( const szt ind1, const szt ind2, szt& m1, szt& e1, szt& m2, szt& e2 ) const noexcept		// the first 2 arguments by value, others by ref !!!!
{
	const auto w1 = glm[ind1];
	const auto w2 = glm[ind2];

	for (const auto j : {1,2})
		if (mt[w1].endInd(j) == ind1)

			for (szt i=1; i<=mt[w1].nn[j]; i++) {
				const auto neib = mt[w1].neib[j][i];
				const auto neen = mt[w1].neen[j][i];

				if (mt[neib].endInd(neen) == ind2)
					for (szt i1=1; i1<=mt[w2].nn[neen]; i1++)

						if (mt[w2].neib[neen][i1] == w1 &&
							mt[w2].neen[neen][i1] == j) {
							m1 = w1; e1 = j;
							m2 = w2; e2 = neen;
							return;
						}
			}
}

template<typename Mt>
void NetworkCore<Mt>::
find_edge( szt& m, szt& a, const szt& ind ) const
{
	for (szt j=1; j<=mtnum; j++)
		for (szt i=0; i<mt[j].g.size(); i++)
			if (mt[j].g[i].ind == ind) {
				m = j;
				a = i;
				return;
			}
	XASSERT(true, "Error: edge not found: "+STR(ind));
}

template<typename Mt>
void NetworkCore<Mt>::
sort_mitos_length( szt* sorted ) const
{
	std::map<szt,szt> mit;
	for (szt i=1; i<=mtnum; i++)
		mit.insert(std::pair<szt,szt>(i, mt[i].g.size()));

	for (szt j=0; j<mtnum; j++) {
		szt maxi = huge<szt>;
		szt maxl = 0;
		std::map<szt,szt>::iterator e = mit.begin();
		while( e != mit.end() ) {
			if (e->second > maxl ) { maxl = e->second; maxi = e->first; }
			e++;
		}
		sorted[j] = maxi;
		mit.erase(mit.find(maxi));
	}
}

}

#endif /* NetworkCore */
	
	
	
	
	
	
	
	
