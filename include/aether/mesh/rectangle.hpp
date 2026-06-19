// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <vector>

#include "./mesh.hpp"

namespace aether::mesh {
class Rectangle : public Mesh {
 public:
  // degree is not used for this simple mesh, but we include it for consistency with the base class.
  // num nodes is the total number of nodes along the horizontal side of the quad
  // domain is a 2 element list specifcying an x interval and a y interval and the
  // mesh is generated as a uniform grid on the Cartesian product of these intervals.
  Rectangle(int degree, int num_nodes, std::vector<Vec2> domain);
  ~Rectangle() override = default;
  std::vector<std::array<NodeIndex, 3>> triangles() const { return triangles_; }
  std::array<NodeIndex, 3> element_node_indices(int element_index) const override;
  std::vector<Vec2> element_nodes(int element_index) const override;
  std::vector<NodeIndex> boundary_nodes() const override;
  std::vector<Vec2> domain() const { return domain_; }

 private:
  // Rectangle-specific data and methods can be added here.
  std::vector<std::array<NodeIndex, 3>> triangles_;
  std::vector<Vec2> domain_;
};
}  // namespace aether::mesh
