add_library(aether_input input.cpp)
target_include_directories(aether_input PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_input PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed