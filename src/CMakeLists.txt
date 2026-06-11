add_executable(aether main.cpp)

target_include_directories(aether PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether PRIVATE spdlog::spdlog_header_only Eigen3::Eigen)

if(AETHER_CLANG_TIDY)
  set_target_properties(aether PROPERTIES CXX_CLANG_TIDY "${AETHER_CLANG_TIDY}")
endif()
