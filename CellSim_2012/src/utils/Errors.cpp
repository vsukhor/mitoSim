
#include "Errors.h"

namespace Utils {

long long assert_fun( const char* EX, 
				 const char *file, 
				 int line, 
				 const std::string& msg ) 
{
	std::cerr << "Assertion (" + std::string(EX) + ") failed! \n" +
			   "File " + file + ", Line " + std::to_string(line) + "\n" +
			   "Reason: " + msg;
	std::abort();
}

int simple_error( const std::string& s )
{
	std::cout << s << std::endl;
	return EXIT_FAILURE;
}

}
