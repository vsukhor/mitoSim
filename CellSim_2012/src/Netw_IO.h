
#ifndef Netw_IO_h
#define Netw_IO_h

namespace MitoD {

template<bool Src, typename Ntw>
class Netw_IO {

public:

	explicit Netw_IO(Ntw&);

	void save_mitos(const bool, const bool, const szt, const real) const;
	void save_edges(const bool, const std::string&, const bool, const real&, std::vector<typename Ntw::ST::ThisEdge>) const;

	void print_clusters() const;

private:

	Ntw&							host;
	Oel&							oel;
	const MitoD::ConfigMain&		cfgMain;
	const std::string				runName;
	const typename Ntw::Reticulum&	mt;
	szt&							mtnum;
	szt&							clnum;
	szt&							mtmass;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<bool Src, typename Ntw>
Netw_IO<Src,Ntw>::
Netw_IO( Ntw& host )
	: host {host}
	, oel {host.oel}
	, cfgMain {host.cfgMain}
	, runName {host.runName}
	, mt {host.mt}
	, mtnum {host.mtnum}
	, clnum {host.clnum}
	, mtmass {host.mtmass}
{}

template<bool Src, typename Ntw>
void Netw_IO<Src,Ntw>::
save_mitos( const bool startnew, const bool last, const szt itr, const real t ) const
{
	const auto fname = (last) ? cfgMain.workingDirOut+"mitos_last_"+runName
							  : cfgMain.workingDirOut+"mitos_"	   +runName;
	const auto flags = (startnew) ? std::ios::out | std::ios::binary | std::ios::trunc
								  : std::ios::out | std::ios::binary | std::ios::app;
	std::ofstream ofs {fname, flags};
	if (ofs.fail())
		oel.print("Cannot open file: "+fname);

	ofs.write((char*) &t, sizeof(real));						//	cout << t << " ";
	ofs.write((char*) &mtnum, sizeof(szt));						//	cout << mtnum << endl;

	static szt mtnummax, nn1max, nn2max;
	if (!last) {
		if (startnew) {
			mtnummax = 0;
			nn1max = 0;
			nn2max = 0;
		}
		if (mtnum > mtnummax)
			mtnummax = mtnum;
	}
	for (szt q=1; q<=mtnum; q++) {
		mt[q].write( ofs );										//	cout << endl;
		if (!last) {
			if (mt[q].nn[1] > nn1max) nn1max = mt[q].nn[1];
			if (mt[q].nn[2] > nn2max) nn2max = mt[q].nn[2];
		}
	}
	ofs.write((char*) &mtnummax, sizeof(szt));					//	cout << endl << mtnummax << " ";
	ofs.write((char*) &nn1max, sizeof(szt));					//	cout << nn1max << " ";
	ofs.write((char*) &nn2max, sizeof(szt));					//	cout << nn2max << " ";

	szt nst2save;
	nst2save = last ? szt(0) : szt(itr/host.saveFrequency);
	ofs.write((char*) &nst2save, sizeof(szt));				//	cout << nst2save << " ";
}

template<bool Src, typename Ntw>
void Netw_IO<Src,Ntw>::
save_edges(
	const bool startnew,
	const std::string& name,
	const bool last,
	const real& t,
	std::vector<typename Ntw::ST::ThisEdge> p ) const
{
	std::string fname = (last) ? cfgMain.workingDirOut+name +"_last_"+runName
							   : cfgMain.workingDirOut+name +"_"	 +runName;
	const auto flags = (startnew) ? std::ios::out | std::ios::binary | std::ios::trunc
								  : std::ios::out | std::ios::binary | std::ios::app;
	std::ofstream ofs {fname, flags};
	if (ofs.fail())
		oel.exit("Cannot open file: "+fname);

	ofs.write((char*)&t, sizeof(real));
	szt k = p.size();
	ofs.write((char*)&k, sizeof(szt));

	for (const auto& o : p)
		o.write(ofs);
}

template<bool Src, typename Ntw>
void Netw_IO<Src,Ntw>::
print_clusters() const
{
	for (szt i=1; i<=mtnum; i++)
		std::cout << i << " " << mt[i].cl << std::endl;
}

}

#endif /* Netw_IO_h */
