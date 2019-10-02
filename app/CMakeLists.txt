add_executable(mitosim main)
target_compile_features(mitosim PRIVATE cxx_std_17)

target_include_directories(mitosim PUBLIC ../include
										  ../include/reactions
										  /Users/vs/Mitofraps/utils)

# This depends on (header only) boost
target_link_libraries(mitosim PRIVATE Boost::boost)
target_link_libraries(mitosim PRIVATE /Users/vs/Mitofraps/utils/lib/libutils_CPU.a)

