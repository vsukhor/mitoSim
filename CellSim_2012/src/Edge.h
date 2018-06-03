#ifndef Edge_h
#define Edge_h

namespace MitoD {

template<int ContentT>
class Edge {

public:

	static const int ConT {ContentT};
	static real		  length;

	szt	ind {huge<szt>};		// index: starts from 0
	szt	indcl {huge<szt>};		// index: starts from 0
	szt	cl {huge<szt>};			// current cl
	szt	nnb_ini {0};			// # of neib edges at some time point
	szt	type {huge<szt>};

	std::array<ulong,2>	fin {0};
	std::array<bool,2>	cin {false};
	real				t_prev {huge<real>};	// time the last motion by pulling took place
	szt					clini {huge<szt>};

	Edge() {}
	Edge(szt,szt,szt);

	void reflect();
	void shiftContent(Edge&) {}						// shifts the edge content to t
	void shiftContent(real,Edge&) {}				// shifts fraction q of the edge content to t
	void shiftContent(std::vector<Edge*>&) {}		// shifts the edge content to a set of particles t
	bool spreadContent(Edge&) { return false; };	// returns true if some content is still left to spread

	void check(szt) const;

	void write(std::ofstream&) const;
	void print(const szt, const bool) const;
	void print(std::ostream&, const szt, const bool) const;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int ContentT> real Edge<ContentT>::length {zero<real>};

template<int ContentT>
Edge<ContentT>::
Edge(
	szt ind,
	szt indcl,
	szt cl
	)
	: ind {ind}
	, indcl {indcl}
	, cl {cl}
{}

template<int ContentT>
void Edge<ContentT>::
reflect()
{
	std::swap(fin[0], fin[1]);
	std::swap(cin[0], cin[1]);
}

template<int ContentT>
void Edge<ContentT>::
check( szt mtmass ) const
{
	if (ind >= mtmass) {
		print(huge<szt>, true);
		std::cout << "check_particles(): Edge.pm is out of range; ind, mtmass" << ind << mtmass << std::endl;
	}
}

template<int ContentT>
void Edge<ContentT>::
write( std::ofstream &ofs ) const
{
	ofs.write((char*) &ind, sizeof(szt));
}

template<int ContentT>
void Edge<ContentT>::
print( const szt a, const bool finish ) const
{
	print(std::cout, a, finish);
}

template<int ContentT>
void Edge<ContentT>::
print( std::ostream& os, const szt a, const bool finish ) const
{
	os << "[" << a << "] ";
	os << " ind " << ind; 
	os << " indcl " << indcl;
	os << " fin " << fin[0] << " " << fin[1];
	os << " cin " << cin[0] << " " << cin[1];
	if (finish) os << "\n";
}

}

#endif /* Edge_Spaceless_Empty_h */

