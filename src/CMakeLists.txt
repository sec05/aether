add_subdirectory(mesh)
add_subdirectory(elements)
add_subdirectory(assembly)
add_subdirectory(solver)

add_executable(aether main.cpp)
target_include_directories(aether PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether PRIVATE aether_mesh aether_solver aether_elements aether_assembly spdlog::spdlog_header_only Eigen3::Eigen)