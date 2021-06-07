cmake_minimum_required(VERSION 3.15)

project(mitosim
    VERSION 0.76.0
    DESCRIPTION "Simulator of mitochoondria reticulum"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O3 -msse4 -Wall)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_VERBOSE_MAKEFILE OFF)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    option(BUILD_DOC "Build documentation" OFF)
    if(BUILD_DOC)
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
endif()

find_package(Boost REQUIRED)

# Dependence on utils:
add_subdirectory(external)

# The executable code is here
add_subdirectory(app)

option(BUILD_TESTS "Build test files" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

