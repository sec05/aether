add_library(aether_assembly assembler.cpp)
target_include_directories(aether_assembly PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_assembly PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed