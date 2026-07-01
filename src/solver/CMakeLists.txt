add_library(aether_solver solver.cpp time_stepper.cpp)
target_include_directories(aether_solver PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(aether_solver PUBLIC Eigen3::Eigen)   # if Vec2 is Eigen-backed