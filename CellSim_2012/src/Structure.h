#ifndef Structure_h
#define Structure_h

namespace MitoD {

template<typename Mt>
class Structure {

public:

	typedef std::vector<Mt> Reticulum;

	vec3<szt>					clagl;
	vec3<szt>					clagl_ini;

	std::vector<szt>			glm, gla;
	vec2<szt>					glc;

	std::array<szt,Mt::maxnn>	nn {{0}};		// node numbers
	Reticulum					mt;
	szt							mtnum {0};		// actual number of mitos
	szt							clnum {0};		// actual number of clusters
	szt							mtmass {0};

	vec2<szt>					clmt;			// mito indices segregated into clusters: clmt - total;
	std::vector<szt>			cls;			// cluster sizes in edges
	vec2<bool>					clvisited;

	std::vector<szt>				mt11;		std::vector<szt>		mtc11;
	std::vector<szt>				mt22;		std::vector<szt>		mtc22;
	std::vector<szt>				mt33;		vec2<szt>				mtc33;
	std::vector<std::array<szt,2>>	mt13;		vec2<std::array<szt,2>>	mtc13;

	Oel&						oel;
	const szt					minLoopLength {2};	// minimal length of a mito that can bend into a loop;
	bool						verbose;

	explicit Structure(const MitoD::ConfigMain&, Oel&);

	void basic_update() noexcept;
	void update_adjacency() noexcept;
	void update_structure() noexcept;
	void update_structure_ind( const szt ind ) noexcept;

	void make_indma() noexcept;

	void make_adjacency_list_edges( const szt ic, vec2<szt>& a ) noexcept;

	void populate_cluster_vectors() noexcept;

	szt count_nnodes(const szt) const noexcept;
	void update_nn(const szt) noexcept;
	void update_nn() noexcept;
	szt find_mtmass() const noexcept;
	void print_mitos(const std::string&) const;
	void print(std::ostream&) const;
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
Structure<Mt>::
Structure( const MitoD::ConfigMain& cfgMain,
		   Oel& oel
	)
	: oel {oel}
	, verbose {cfgMain.verbose}
//	, check {*this}
{}

template<typename Mt> inline
void Structure<Mt>::
update_adjacency() noexcept
{
	if (clagl.size() < clnum) {
		clagl.resize(clnum);
		clvisited.resize(clnum);
	}
	for (szt ic=0; ic<clnum; ic++)
		make_adjacency_list_edges(ic, clagl[ic]);
}

template<typename Mt> inline
void Structure<Mt>::
basic_update() noexcept
{
	make_indma();
	populate_cluster_vectors();
}

template<typename Mt> inline
void Structure<Mt>::
update_structure() noexcept
{
	basic_update();
	update_adjacency();
}

template<typename Mt> inline
void Structure<Mt>::
update_structure_ind( const szt ind ) noexcept
{
	if (clagl.size() < clnum) {
		clagl.resize(clnum);
		clvisited.resize(clnum);
	}
	const auto ic = mt[glm[ind]].cl;
	make_adjacency_list_edges(ic, clagl[ic]);
}

template<typename Mt> inline
void Structure<Mt>::
make_indma() noexcept
{
	cls.resize(clnum);
	std::fill(cls.begin(), cls.end(), 0);			// cluster size, # of edges
	for (szt j=1; j<=mtnum; j++)
		cls[mt[j].cl] += mt[j].g.size();	

	glc.resize(clnum);
	for (szt j=0; j<clnum; j++)
		glc[j].resize(cls[j]);

	glm.resize(mtmass);
	gla.resize(mtmass);
	for (szt j=1; j<=mtnum; j++)
		for (szt k=0; k<mt[j].g.size(); k++) {
			const auto& g {mt[j].g[k]};
			glm[g.ind] = j;
			gla[g.ind] = k;
			glc[g.cl][g.indcl] = g.ind;
		}	
}
template<typename Mt>
void Structure<Mt>::
make_adjacency_list_edges( const szt ic, vec2<szt>& a ) noexcept
{
	clvisited[ic].resize(cls[ic]);

	a.resize(cls[ic]);
	for (auto& o : a)	o.clear();

	for (const auto& j : clmt[ic])
		for (szt k=0; k<mt[j].g.size(); k++) {
			const auto ind = mt[j].g[k].indcl;
			if (k == 0) {
				for (szt ie=1; ie<=mt[j].nn[1]; ie++) {				// connection backwards: only other mitos might be found
					const auto w2 = mt[j].neib[1][ie];
					const auto a2 = mt[w2].end2a(mt[j].neen[1][ie]);
					a[ind].push_back( mt[w2].g[a2].indcl );
				}
				if (mt[j].g.size() == 1)								// connection forwards: to other mito
					for (szt ie=1; ie<=mt[j].nn[2]; ie++) {
						const auto w2 = mt[j].neib[2][ie];
						const auto a2 = mt[w2].end2a(mt[j].neen[2][ie]);
						a[ind].push_back(mt[w2].g[a2].indcl);
					}
				else {													// connection forwards: to the same mito
					a[ind].push_back( mt[j].g[k+1].indcl );
				}
			}
			else if (k == mt[j].g.size()-1) {							// but not  a1 == 0  =>  mt[m1].g.size() > 1
				a[ind].push_back(mt[j].g[k-1].indcl);					// connection backwards: to the same mito
				for (szt ie=1; ie<=mt[j].nn[2]; ie++) {				// connection forwards: to other mito
					const auto w2 = mt[j].neib[2][ie];
					const auto a2 = mt[w2].end2a( mt[j].neen[2][ie]);
					a[ind].push_back(mt[w2].g[a2].indcl);
				}
			}
			else {														// edge in the bulk: a1 != 1 && a1 != mt[m1].g.size()
				a[ind].push_back(mt[j].g[k-1].indcl);					// connection backwards: to the same mito
				a[ind].push_back(mt[j].g[k+1].indcl);					// connection forwards: to the same mito
			}
		}
}

template<typename Mt>
void Structure<Mt>::
populate_cluster_vectors() noexcept
{
	mt11.clear();	mtc11.resize(clnum);		std::fill(mtc11.begin(),  mtc11.end(),  huge<szt>);
	mt22.clear();	mtc22.resize(clnum);		std::fill(mtc22.begin(),  mtc22.end(),  huge<szt>);
	mt33.clear();	mtc33.resize(clnum);		for (auto& o : mtc33) o.clear();
	mt13.clear();	mtc13.resize(clnum);		for (auto& o : mtc13) o.clear();
	nn = {{zero<szt>}};
	clmt.resize(clnum);		for (auto& o : clmt) o.clear();		// # of segments
	
	for (szt j=1; j<=mtnum; j++) {
		const auto& m {mt[j]};
		clmt[m.cl].push_back(j);						// mitochondrial indexes clusterwise
		nn[1] += m.nnodes(2);

		const auto e {m.hasOneFreeEnd()};
		if (e) {
			const szt oe {e == 1 ? static_cast<szt>(2) : static_cast<szt>(1)};
			nn[0]++;
			if (m.nn[oe] == 2) {
				const std::array<szt,2> je {j, e};
				mtc13[m.cl].emplace_back(je);			// mito index, free end index
				mt13.emplace_back(je);
				nn[2]++;
			}
		}
		else if (m.nn[1] == 0 && m.nn[2] == 0) {
			mtc11[m.cl] = j;							// it is a separate segment since it has two free ends
			mt11.push_back(j);							// it is a separate segment since it has two free ends
			nn[0] += 2;
		}
		else if (m.isPureLoop()) {
			mtc22[m.cl] = j;
			mt22.push_back(j);							// it is a separate segment since it has two free ends
		}
		else if (m.nn[1] == 2 && m.nn[2] == 2) {
			mtc33[m.cl].push_back(j);
			mt33.push_back(j);
			nn[2] += 2;
		}
		else {
			; XASSERT(true, "Error in populate_cluster_vectors: failed classification for "+STR(j));
		}
	}
	if (false) {
		if (verbose) { std::cout << "mt11 "; for (const auto o : mt11) std::cout << o << " "; std::cout << std::endl; }
		if (verbose) { std::cout << "mt22 "; for (const auto o : mt22) std::cout << o << " "; std::cout << std::endl; }
		if (verbose) { std::cout << "mt33 "; for (const auto o : mt33) std::cout << o << " "; std::cout << std::endl; }
		if (verbose) { std::cout << "mt13 "; for (const auto& o : mt13) std::cout << o[0] << " "; std::cout << std::endl; }
	}
	nn[2] /= 3;
}

template<typename Mt>
szt Structure<Mt>::
count_nnodes( const szt deg ) const noexcept
{
	szt k {0};
	for (szt i=1; i<=mtnum; i++)
		k += mt[i].nnodes(deg);
	if (deg == 3) return k/3;
	else		       return k;
}

template<typename Mt>
void Structure<Mt>::
update_nn( const szt deg ) noexcept
{
	nn[deg-1] = count_nnodes(deg);
}

template<typename Mt>
void Structure<Mt>::
update_nn() noexcept
{
	update_nn(1);
	update_nn(2);
	update_nn(3);
}

template<typename Mt>
void Structure<Mt>::
print_mitos( const std::string& tag ) const
{
	for (szt j=1; j<=mtnum; j++)
		mt[j].print(j, tag, -1);
	oel.print("");
}
template<typename Mt>
szt Structure<Mt>::
find_mtmass() const noexcept
{
	szt mass {0};
	for (szt j=1; j<=mtnum; j++)
		mass += mt[j].length();

	return mass;
}

template<typename Mt>
void Structure<Mt>::
print( std::ostream& ofs ) const
{
	ofs << " X ";
	for (const auto o : nn)
		ofs << o << " ";
	ofs <<  "m11 " << mt11.size()
		<< " m22 " << mt22.size()
		<< " m33 " << mt33.size()
		<< " m13 " << mt13.size()
		<< " mtm " << mtmass
		<< " mtn " << mtnum
		<< " cln " << clnum;
}

}
#endif /* Structure_h */
