add_executable(mitosim_cpu main.cpp)
target_compile_features(mitosim_cpu PRIVATE cxx_std_20)

target_include_directories(mitosim_cpu PUBLIC
                           ../include
                           ../include/reactions
                           ../external)

target_link_libraries(mitosim_cpu PRIVATE Boost::boost)
target_link_libraries(mitosim_cpu PRIVATE $<TARGET_FILE:utils>)

add_custom_command(
    TARGET mitosim_cpu
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mitosim_cpu> ${CMAKE_CURRENT_SOURCE_DIR}
)
