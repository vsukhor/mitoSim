/* ==============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann.
   All Rights Reserved.
   Developed at Helmholtz Center for Infection Research, Braunschweig, Germany.
   Please see Readme file for further information

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

============================================================================== */

#ifndef Gillespie_h
#define Gillespie_h

#include "Misc.h"

namespace Utils {

template<typename Reaction> 
class Gillespie {

public :

	typedef real realT;
	
	Gillespie(RandFactory&) noexcept;

	void add_reaction(std::unique_ptr<Reaction>);
	void initialize() noexcept;
	bool set_asum() noexcept;
	void make() noexcept;
	realT tau() const noexcept { return tau_; }
	void printScores(std::ostream &) const;
	void log_data(std::ostream&) const;
	constexpr szt num_reactions() const noexcept;
	std::string srt(const szt) const noexcept;

private:

	vup<Reaction>			rc;
	std::vector<realT>		a;
	szt						rind {huge<szt>};
	realT					tau_ {};
	RandFactory&			rnd;
	szt						nreact {};
	std::vector<szt>		rtype;
	std::vector<realT>		csums;
	realT					asum {};
	std::vector<realT>		auxf;
	std::vector<szt>		auxi;
	std::vector<szt>		rinds;

	void set_tau() noexcept;
	void set_rind() noexcept;
	constexpr Reaction* currRc() const noexcept;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Reaction> 
Gillespie<Reaction>::
Gillespie( RandFactory& rnd ) noexcept
	: rnd(rnd)
{}

template<typename Reaction>
void Gillespie<Reaction>::
add_reaction( std::unique_ptr<Reaction> rup )
{
	rc.push_back(std::move(rup));
}

template<typename Reaction>
void Gillespie<Reaction>::
initialize() noexcept
{
	for (szt ri=0; ri<rc.size(); ri++)
		rtype.push_back(ri);

	nreact = rtype.size();
	a.resize(nreact);
	for (szt i=0; i<nreact; i++) {
		rc[rtype[i]]->attach_score_pointer(&a[i]);
		rc[rtype[i]]->initialize_dependencies(rc);
	}
	
	auxf.resize(nreact);
	csums.resize(nreact);
	auxi.resize(nreact);
	rinds.resize(nreact);
}

template<typename Reaction> 
bool Gillespie<Reaction>::
set_asum() noexcept
{
	asum = std::accumulate(a.begin(), a.end(), zero<realT>);
	
	if (asum == zero<realT>)
		return false;
	return true;
}

template<typename Reaction> 
void Gillespie<Reaction>::set_rind() noexcept
{
	for (szt i=0; i<nreact; i++)
		auxf[i] = a[i] / asum;
	std::partial_sum(auxf.begin(), auxf.end(), csums.begin(), std::plus<realT>());

	realT ran;
	do ran = rnd.r01u(); 
	while (ran >= csums[nreact-1]);

	for (szt i=0; i<nreact; i++)
		auxi[i] = szt(ran < csums[i]);

	const szt rindnum = Utils::find(auxi, rinds);
	rind = *std::min_element(rinds.begin(), rinds.begin()+rindnum);
}

template<typename Reaction> 
void Gillespie<Reaction>::
set_tau() noexcept
{
	realT ran;
	do ran = rnd.r01u(); 
	while (ran <= zero<realT> || ran >= one<realT>);

	tau_ = std::log(one<realT>/ran) / asum;
}

template<typename Reaction> 
void Gillespie<Reaction>::
make() noexcept
{
	set_rind();
	set_tau();
	(*currRc())();
}

template<typename Reaction> constexpr
Reaction* Gillespie<Reaction>::
currRc() const noexcept
{
	return rind < huge<szt> ? rc[rtype[rind]].get()
						    : nullptr;
}

template<typename Reaction> constexpr
szt Gillespie<Reaction>::
num_reactions() const noexcept
{
	return rc.size();
}

template<typename Reaction>
std::string Gillespie<Reaction>::
srt( const szt i ) const noexcept
{
	return rc[i]->srt;
}

template<typename Reaction>
void Gillespie<Reaction>::
log_data( std::ostream& ofs ) const
{
	ofs << " tau " << tau_
		<< " rt " << (rind!=huge<szt> ? "" : "000") << padZeros3(rind);
	if (rind != huge<szt>)
		ofs << " " << currRc()->srt;  

}

template<typename Reaction>
void Gillespie<Reaction>::
printScores( std::ostream& ofs ) const
{
	for (szt i=0; i<nreact; i++)
		rc[i]->set_score();
	ofs << " scores: "; 
	for (szt i=0; i<nreact; i++)
		ofs << rc[i]->srt << " "
			<< rc[i]->get_eventCount() << " "
			<< rc[i]->get_score() << " ";
}

}

#endif /* Gillespie_h */
