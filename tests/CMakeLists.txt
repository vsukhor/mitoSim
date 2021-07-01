include(FetchContent)
set(FETCHCONTENT_QUIET off)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

include(GoogleTest)

add_definitions(-DUSE_UTILS_XASSERT)

add_executable(unittest1
  test_config.cpp
  test_edge.cpp
  test_segment.cpp
  test_structure.cpp
)
add_executable(unittests2
  test_core_transformer.cpp
)
target_compile_features(unittests1 PRIVATE cxx_std_20)
target_compile_features(unittests2 PRIVATE cxx_std_20)

include_directories(unittests1 PUBLIC
    ../include
    ../include/reactions
    ../external
)
include_directories(unittests2 PUBLIC
    ../include
    ../include/reactions
    ../external
)

target_link_libraries(unittests1 PRIVATE $<TARGET_FILE:utils>)
target_link_libraries(unittests1 PRIVATE gtest gmock)
target_link_libraries(unittests1 PRIVATE gtest_main)
target_link_libraries(unittests2 PRIVATE $<TARGET_FILE:utils>)
target_link_libraries(unittests2 PRIVATE gtest gmock)
target_link_libraries(unittests2 PRIVATE gtest_main)

gtest_discover_tests(unittests1)
gtest_discover_tests(unittests2)

add_custom_command(
    TARGET unittests1
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unittests1> ${CMAKE_CURRENT_SOURCE_DIR}
)
add_custom_command(
    TARGET unittests2
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unittests2> ${CMAKE_CURRENT_SOURCE_DIR}
)

