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

add_executable(unittests
  test_config.cpp
  test_edge.cpp
  test_segment.cpp
  test_structure.cpp
  test_core_transformer.cpp
  test_ability_for_fusion.cpp
  test_ability_for_fission.cpp
  test_network.cpp
  test_ntw_fusion_12.cpp
)

target_compile_features(unittests PRIVATE cxx_std_20)

include_directories(unittests PUBLIC
    ../include
    ../include/reactions
    ../external
)

target_link_libraries(unittests PRIVATE $<TARGET_FILE:utils>)
target_link_libraries(unittests PRIVATE gtest)
target_link_libraries(unittests PRIVATE gtest_main)

gtest_discover_tests(unittests)

add_custom_command(
    TARGET unittests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unittests> ${CMAKE_CURRENT_SOURCE_DIR}
)

