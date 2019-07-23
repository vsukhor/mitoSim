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
