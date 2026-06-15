file(GLOB AETHER_TEST_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

include(GoogleTest)

foreach(test_src ${AETHER_TEST_SOURCES})
  get_filename_component(test_name ${test_src} NAME_WE)
  add_executable(${test_name} ${test_src})
  target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries(${test_name} PRIVATE aether_mesh aether_elements aether_assembly GTest::gtest_main Eigen3::Eigen)
  if(AETHER_CLANG_TIDY)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${AETHER_CLANG_TIDY}")
  endif()
  gtest_discover_tests(${test_name})
endforeach()