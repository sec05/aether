add_library(aether_elements p1_element.cpp)
target_include_directories(aether_elements PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_elements PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed