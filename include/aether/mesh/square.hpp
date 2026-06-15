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
  std::array<NodeIndex, 3> element_node_indices(int element_index) const override;
  std::vector<Vec2> element_nodes(int element_index) const override;

 private:
  // Square-specific data and methods can be added here.
  std::vector<std::array<NodeIndex, 3>> triangles_;
};
}  // namespace aether::mesh
