add_library(aether_mesh rectangle.cpp quad.cpp)
target_include_directories(aether_mesh PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_mesh PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed