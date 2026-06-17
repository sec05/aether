// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <string>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "../core/geometry.hpp"
#include "../core/types.hpp"
#include <pugixml.hpp>

namespace aether::output {
class Output {
 public:
  void write_solution_vtk(
      const std::string& path, const std::vector<Vec2>& points,
      const std::vector<std::array<NodeIndex, 3>>& triangles,
      const std::vector<std::pair<std::string, Eigen::VectorXd>>& point_fields) const;
};

}  // namespace aether::output
