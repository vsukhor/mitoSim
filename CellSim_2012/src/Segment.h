#ifndef Segment_h
#define Segment_h

#include "Edge.h"

namespace MitoD {

template<int ContentT>
class Segment {

public:
	static const int ConT {ContentT};

	typedef Edge<ContentT> ThisEdge;
	typedef Segment<ContentT> thisT;


	Oel&							oel;
	std::vector<ThisEdge>			g;				// edges
	szt								cl {0};			// cluster index
	std::array<szt,3>				nn {{0}};		// number of neibours
	std::array<std::vector<szt>,3>	neib;			// niebors
	std::array<std::vector<szt>,3>	neen;			// neibour ends
	static const szt				maxnn {3};

	Segment(Oel&);
	Segment(
		  const szt segmass,
		  const szt cl,
		  szt& mtmass, 
		  szt& ei, 
		  Oel& oel );
		  
	ThisEdge* increment_length( const long a, const ThisEdge p );					// inserts a particle imediately after g[a] making it g[a+1]
	void reflect_g();

	szt set_gCl( const szt newcl, const szt initind );
	szt setCl( const szt newcl, const szt initind );

	constexpr szt end2a( const szt& e ) const;
	szt endInd( const szt& e ) const;
	constexpr szt hasOneFreeEnd() const;			// return the end index if true
	szt singleNeibInd( const szt& e ) const;
	std::vector<szt> doubleNeibInds( const szt& e ) const;
	bool hasSuchNeib( const szt& e, const szt& n ) const;
	constexpr bool isPureLoop() const;
	constexpr bool isLoop( const szt& i ) const;

	szt nnodes( const szt& deg ) const;
	szt length() const noexcept { return g.size(); }

	ulong set_end_fin( const szt a );
	ulong set_bulk_fin( const szt a );

	void print( const szt w, const std::string& tag, const szt at = huge<szt> ) const;
	void print( std::ostream& os, const szt w, const std::string& tag, const szt at = huge<szt> ) const;
	void write( std::ofstream& ) const;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template<int ContentT> const int Segment<ContentT>::ConT;
template<int ContentT> const szt Segment<ContentT>::maxnn;

template<int ContentT>
Segment<ContentT>::
Segment( Oel& oel )
	: oel(oel)
{
	neib[1].resize(maxnn+1);
	neib[2].resize(maxnn+1);
	
	neen[1].resize(maxnn+1);
	neen[2].resize(maxnn+1);
}

template<int ContentT>
Segment<ContentT>::
Segment(
	  const szt segmass,
	  const szt cl,
	  szt& mtmass, 
	  szt& ei, 
	  Oel& oel
	)
	: oel(oel),
	  cl(cl)
{
	neib[1].resize(maxnn+1);	neib[2].resize(maxnn+1);
	neen[1].resize(maxnn+1);	neen[2].resize(maxnn+1);

	for (szt a=0; a<segmass; a++)
		increment_length(long(a-1), {ei++, a, cl});

	mtmass += segmass;

	for (auto& o : g) o.clini = cl;
}

template<int ContentT>
typename Segment<ContentT>::ThisEdge* Segment<ContentT>::
increment_length( const long a, const Segment<ContentT>::ThisEdge p ) 	// inserts a particle imediately after g[a] making it g[a+1]
{
	g.insert(g.begin()+a+1, p);
	return &g[a+1];
}

template<int ContentT> inline
void Segment<ContentT>::
reflect_g()
{
	const std::vector<ThisEdge> temp = g;
	for (szt i=0; i<g.size(); i++) { 
		g[i] = temp[g.size()-i-1];
		g[i].reflect();
	}
}

template<int ContentT> inline
szt Segment<ContentT>::
set_gCl( const szt newcl, const szt initind )
{
	for (szt i=0; i<g.size(); i++) {
		g[i].cl = newcl;
		g[i].indcl = initind + i;
	}
	
	return initind + (szt)g.size();
}

template<int ContentT> inline
szt Segment<ContentT>::
setCl( const szt newcl, const szt initind )
{
	cl = newcl;
	return set_gCl(newcl, initind);
}

template<int ContentT> constexpr
szt Segment<ContentT>::
end2a( const szt& e ) const 
{
	return (e == 1) ? 0 : (szt)g.size()-1;
}

template<int ContentT> inline
szt Segment<ContentT>::
endInd( const szt& e ) const
{
	return (e == 1) ? g.front().ind
					: g.back().ind;
}

template<int ContentT> constexpr
szt Segment<ContentT>::
hasOneFreeEnd() const 			// return the end index if true
{
	if (	!nn[1] &&  nn[2])  	return 1;
	else if (nn[1] && !nn[2])	return 2;
	else						return 0;
}

template<int ContentT> inline
szt Segment<ContentT>::
singleNeibInd( const szt& e ) const 
{
	for (szt i=1; i<=nn[e]; i++) 
		if (neib[e][i])
			return i;
			
	return huge<szt>;
}

template<int ContentT> inline
std::vector<szt> Segment<ContentT>::
doubleNeibInds( const szt& e ) const
{
	XASSERT(nn[e] == 2, "Error in Mito::doubleNeibInds: nn[e] != 2 in cluster "+STR(cl));

	std::vector<szt> neibInds(nn[e]);
	szt j = 0;
	for (szt i=1; i<=nn[e]; i++)
		if (neib[e][i])
			neibInds[j++] = i;
	
	return neibInds;
}

template<int ContentT> inline
bool Segment<ContentT>::
hasSuchNeib( const szt& e, const szt& n ) const
{
	for (szt i=1; i<=nn[e]; i++)
		if (neib[e][i] == n)
			return true;
	
	return false;
}

template<int ContentT> constexpr
bool Segment<ContentT>::
isPureLoop() const 
{
	return nn[1] == 1 && 
		   nn[2] == 1 && 
		   neib[1][singleNeibInd(1)] == neib[2][singleNeibInd(2)];
}

template<int ContentT> constexpr
bool Segment<ContentT>::
isLoop( const szt& i ) const 
{
	return isPureLoop() || 
		   hasSuchNeib(1, i) ||
		   hasSuchNeib(2, i);
}

template<int ContentT>
szt Segment<ContentT>::
nnodes( const szt& deg ) const				// deg = 1, 2, 3, 4
{														
	if (deg == 1) {															// count nodes of degree 1
		if (      nn[1] &&  nn[2]) 	return 0;
		else if (!nn[1] && !nn[2]) 	return 2;
		else					 	return 1;
	}
	else if (deg == 2)
		return nn[1] && nn[2] && isPureLoop() ? g.size() : g.size() - 1;	// count nodes of degree 2
	else if (deg == 3) {													//	if (deg == 3 ) count nodes of degree 3
		if (     nn[1] == 2 && nn[2] == 2) return 2;
		else if (nn[1] == 2 || nn[2] == 2) return 1;
		else if (nn[1] != 2 && nn[2] != 2) return 0;
	}
	else if (deg == 4) {													//	if (deg == 3 ) count nodes of degree 3
		if (     nn[1] == 3 && nn[2] == 3) return 2;
		else if (nn[1] == 3 || nn[2] == 3) return 1;
		else if (nn[1] != 3 && nn[2] != 3) return 0;
	}
	else oel.exit("Error in Mito::nnodes. Not implemented deg", deg);
	return huge<szt>;
}

template<int ContentT> inline
ulong Segment<ContentT>::
set_end_fin( const szt e )
{
	auto& f {g[end2a(e)].fin};
	f[e-1] = nn[e] ? 1 : 0;

	return f[e-1];
}


template<int ContentT> inline
ulong Segment<ContentT>::
set_bulk_fin( const szt a )
{
	g[a].fin[1] = g[a+1].fin[0] = 1;

	return g[a].fin[1];
}

template<int ContentT>
void Segment<ContentT>::
print( const szt w,
	   const std::string& tag, 
	   const szt at ) const
{
	print(oel.so, w, tag, at);
	print(oel.sl, w, tag, at);
}

template<int ContentT>
void Segment<ContentT>::
print( std::ostream& os, 
	   const szt w,
	   const std::string& tag, 
	   const szt at ) const
{
	os << "\t" << tag << w;
	if (at == huge<szt>) os << "(of ";
	else os << "(at " << at << " of ";
	os << g.size() << ") [ ";	
	for (szt i=1; i<=nn[1]; i++) os << neib[1][i] << " ";	
	os << "] { ";		
	for (szt i=1; i<=nn[1]; i++) os << neen[1][i] << " ";	
	os << "} [ ";
	for (szt i=1; i<=nn[2]; i++) os << neib[2][i] << " ";	
	os << "] { ";	
	for (szt i=1; i<=nn[2]; i++) os << neen[2][i] << " ";	
	os << "} " << cl;	
#ifdef PRINT_EDGES
	os << std::endl;
	for (szt i=0; i<g.size(); i++) 
		g[i].print(os, i, 1);
#else
	os << " len " << g.size();
#endif
	os << std::endl;
}

template<int ContentT>
void Segment<ContentT>::
write( std::ofstream& ofs ) const
{
	szt len = (szt)g.size();
	ofs.write((char*) &len, sizeof(szt));
	ofs.write((char*) &cl, sizeof(szt));
	
	ofs.write((char*) &nn[1], sizeof(szt));
	
	for (szt j=1; j<=nn[1]; j++) {
		ofs.write((char*) &neib[1][j], sizeof(szt));
		ofs.write((char*) &neen[1][j], sizeof(szt));
	}
	
	ofs.write((char*) &nn[2], sizeof(int));
	
	for (szt j=1; j<=nn[2]; j++) {
		ofs.write((char*) &neib[2][j], sizeof(szt));
		ofs.write((char*) &neen[2][j], sizeof(szt));
	}
//	ofs.write( (char*) &dpsi, sizeof(real));
	for (szt j=0; j<len; j++)
		g[j].write(ofs);
}

}

#endif /* Spaceless_Empty_Segment_h */
