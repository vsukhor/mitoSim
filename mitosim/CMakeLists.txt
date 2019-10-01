cmake_minimum_required(VERSION 3.15)

project(mitosim
    VERSION 0.73
    DESCRIPTION "Simulator of mitochoondria reticulum"
    LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the main CMakeLists.
    include(CTest)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND)
#		add_subdirectory(docs)
		# set input and output files
    	set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/doxyfile.in)
    	set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    	# request to configure the file
    	configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    	message(${CMAKE_CURRENT_SOURCE_DIR})
    	message("Doxygen build started")

    	# note the option ALL which allows to build the docs together with the application
    	add_custom_target(doc_doxygen ALL
        	COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        	COMMENT "Generating API documentation with Doxygen"
        	VERBATIM)
    else()
		message(STATUS "Doxygen not found, not building docs")
    endif()
endif()

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14
if(${CMAKE_VERSION} VERSION_LESS 3.14)
    include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()

#SET(CMAKE_CXX_COMPILER /Library/Developer/CommandLineTools/usr/bin)
#SET(CMAKE_ASM_COMPILER /usr/bin/as)
#add_compile_options(-O3 -S -msse4 -march=native -Wa,-q)
add_compile_options(-O3 -msse4)

# Accumulator library
# This is header only, so could be replaced with git submodules or FetchContent
find_package(Boost REQUIRED)

#find_library(LIB_UTILS_LOCATION /Users/vs/Mitofraps/utils/lib/libutils_CPU.a)

# The compiled library code is here
#add_subdirectory(lib/utils/src)

# The executable code is here
add_subdirectory(app)

# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
###if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
###    add_subdirectory(tests)
###endif()




