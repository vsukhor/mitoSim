#ifndef ConfigMain_h
#define ConfigMain_h

#include <array>
#include <string>


#include "utils/Misc.h"
#include "utils/Par.h"
#include "utils/Errors.h"

namespace MitoD {

class ConfigMain;

class ConfigReader {

public:

	ConfigReader( const std::string& cfgFname,
				  Oel& oel
				)
		: cfgFname {cfgFname}
		, oel {oel}
	{}
	bool operator()( std::string s, const std::vector<bool>& range )
	{
		return Par<bool,true>(s, cfgFname, oel, range)();
	}

	real operator()( std::string s, const std::vector<real>& range )
	{
		return Par<real,false>(s, cfgFname, oel, range)();
	}
	szt operator()( std::string s, const std::vector<szt>& range )
	{
		return Par<szt,false>(s, cfgFname, oel, range)();
	}

private:

	const std::string	cfgFname;
	Oel&				oel;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ConfigMain {

	ConfigMain( const ConfigMain& ) = delete;				// non construction-copyable: single parameter class per simulation
	ConfigMain& operator=( const ConfigMain& ) = delete;	// non copyable: single parameter class per simulation

public:

	const std::string		workingDirIn;
	const std::string		workingDirOut;
	const std::string		configSuffix;
	const std::string		cfgFname;
	const bool				verbose;
	const szt				runIni;
	const szt				runEnd;

private:

	std::string check_fname( const std::string fname )
	{
		if (file_exists(fname))
			return fname;
		else
			exit(simple_error("Error: no config file provided: "+fname));
	}

public:

	std::string confname( const std::string& s ) const noexcept
	{ 
		return workingDirOut+"config"+s +"_"+configSuffix+".tcl";
	}

	ConfigMain( const std::string& workingDir,
				const std::string& configSuffix )
		: workingDirIn {workingDir}
		, workingDirOut {workingDir}
		, configSuffix {configSuffix}
		, cfgFname {check_fname(confname("Main"))}
		, verbose {Par<bool,true>::readin("verbose", cfgFname, {false, true})}
		, runIni {Par<szt,false>::readin("runIni", cfgFname, {0,huge<uint>})}
		, runEnd {Par<szt,false>::readin("runEnd", cfgFname, {runIni,huge<uint>})}
{}

	void print( Oel& oel )
	{
		oel.print("workingDirIn = "+workingDirIn);
		oel.print("workingDirOut = "+workingDirOut);
		oel.print("verbose = %d", verbose);
		oel.print("runIni = %d ", runIni);
		oel.print("runEnd = %d ", runEnd);
}

};

}

#endif /* ConfigMain_h */
