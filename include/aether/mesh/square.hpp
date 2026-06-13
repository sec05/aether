// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <vector>

#include "./mesh.hpp"

namespace aether::mesh {
class Square : public Mesh {
 public:
  Square(int degree, int num_nodes);
  ~Square() override = default;
  std::vector<std::array<NodeIndex, 3>> triangles() const { return triangles_; }

 private:
  // Square-specific data and methods can be added here.
  std::vector<std::array<NodeIndex, 3>> triangles_;
};
}  // namespace aether::mesh
