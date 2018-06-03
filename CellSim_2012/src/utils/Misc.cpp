#include "Misc.h"

namespace Utils {

#ifdef __LINUX__
const std::string SLASH {"/"};
#else
const std::string SLASH {"\\"};
#endif

bool file_exists( const std::string& name )
{
	if (FILE *file = fopen(name.c_str(), "r") ) {
		fclose( file );
		return true;
	} 
	else return false;
}

bool fileExists( const std::string& name )
{
  class stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

std::string get_current_time()
{
	time_t rawtime;
	time(&rawtime);
	
	return std::string(ctime(&rawtime));
}

std::string padZeros3( szt n )
{
	return std::string(n<100 ? (n<10 ? "00" : "0") : "")+STR(n);
}

}
