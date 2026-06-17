add_library(aether_output output.cpp)
target_include_directories(aether_output PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_output PUBLIC Eigen3::Eigen pugixml::pugixml)   # if Vec2 is Eigen-backed