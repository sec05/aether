add_executable(aether_tests test_version.cpp)

target_include_directories(aether_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_tests PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(aether_tests)

if(AETHER_CLANG_TIDY)
  set_target_properties(aether_tests PROPERTIES CXX_CLANG_TIDY "${AETHER_CLANG_TIDY}")
endif()
