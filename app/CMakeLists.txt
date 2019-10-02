add_executable(mitosim_CPU main)
target_compile_features(mitosim_CPU PRIVATE cxx_std_17)

target_include_directories(mitosim_CPU PUBLIC ../include
										  	../include/reactions
										  	../external)

# This depends on (header only) boost
target_link_libraries(mitosim_CPU PRIVATE Boost::boost)
target_link_libraries(mitosim_CPU PRIVATE "${LIB_UTILS}")

add_custom_command(TARGET mitosim_CPU
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mitosim_CPU> ${CMAKE_CURRENT_SOURCE_DIR})
