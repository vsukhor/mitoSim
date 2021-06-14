add_executable(mitosim main.cpp)
target_compile_features(mitosim PRIVATE cxx_std_20)

target_include_directories(mitosim PUBLIC
                           ../include
                           ../include/reactions
                           ../external)

target_link_libraries(mitosim PRIVATE Boost::boost)
target_link_libraries(mitosim PRIVATE $<TARGET_FILE:utils>)

add_custom_command(
    TARGET mitosim
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mitosim> ${CMAKE_CURRENT_SOURCE_DIR}
)
