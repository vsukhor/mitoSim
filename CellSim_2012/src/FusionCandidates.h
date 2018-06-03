
#ifndef FusionCandidates_h
#define FusionCandidates_h

namespace MitoD {

struct FusionCandidatesXX {

	std::vector<std::array<szt,2>> u;
	std::vector<std::array<szt,2>> v;

	void clear() noexcept { u.clear(); v.clear(); }

	void add(const std::array<szt,2>& uc, const std::array<szt,2>& vc ) {
		u.emplace_back(uc);
		v.emplace_back(vc);
	}
	szt size() const noexcept { return u.size(); }

	void print() {
		for(szt i=0; i<size(); i++)
			print(i, false);
		std::cout << "\n";
	}

	void print( const szt i, bool nl=true ) {
			std::cout << " [" <<  u[i][0] << " " <<  u[i][1] << " + "
							  <<  v[i][0] << " " << v[i][1] << "] ";
			if (nl) std::cout << "\n";
	}
};

struct FusionCandidatesXL {

	std::vector<std::array<szt,2>> u;
	std::vector<szt>	 v;

	void clear() noexcept { u.clear(); v.clear(); }

	void add(std::array<szt,2> uc, szt vc ) {
		u.emplace_back(uc);
		v.emplace_back(vc);
	}

	szt size() const noexcept { return u.size(); }
};

}

#endif /* NtwFusionCandidates_h */
