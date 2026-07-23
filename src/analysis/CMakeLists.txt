add_library(aether_analysis error_norms.cpp)
target_include_directories(aether_analysis PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_analysis PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed