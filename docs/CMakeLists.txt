# Docs with Doxygen
find_package(Doxygen)
if(Doxygen_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started: ${CMAKE_CURRENT_SOURCE_DIR}")

    # note the option ALL which allows to build the docs together with
    # the application
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
#        VERBATIM
    )
	add_custom_command(
		TARGET doc_doxygen
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/html/docs/imgs
	)
	add_custom_command(
		TARGET doc_doxygen
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/imgs ${CMAKE_CURRENT_SOURCE_DIR}/html/docs/imgs
	)
else()
    message(STATUS "Doxygen not found, not building docs")
endif()
