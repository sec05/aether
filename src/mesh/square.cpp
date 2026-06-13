// Copyright 2026 Spencer Evans-Cole

#include "aether/mesh/square.hpp"

namespace aether::mesh {
Square::Square(int degree, int num_nodes) : Mesh(degree, num_nodes) {
  int num_nodes_per_side = static_cast<int>(std::sqrt(num_nodes));
  if (num_nodes_per_side * num_nodes_per_side != num_nodes) {
    throw std::invalid_argument("num_nodes must be a perfect square for a square mesh.");
  }
  // Generate nodes for a uniform square mesh on [0,1] x [0,1].
  nodes_.reserve(num_nodes);
  for (int j = 0; j < num_nodes_per_side; ++j)
    for (int i = 0; i < num_nodes_per_side; ++i) {
      nodes_.emplace_back(static_cast<Real>(i) / (num_nodes_per_side - 1),
                          static_cast<Real>(j) / (num_nodes_per_side - 1));
    }

  // Generate triangle connectivity for the square mesh
  triangles_.reserve(
      static_cast<std::size_t>(2 * (num_nodes_per_side - 1) * (num_nodes_per_side - 1)));
  for (Index j = 0; j < num_nodes_per_side - 1; ++j) {
    for (Index i = 0; i < num_nodes_per_side - 1; ++i) {
      const Index n0 = j * num_nodes_per_side + i;
      const Index n1 = n0 + 1;
      const Index n2 = n0 + num_nodes_per_side;
      const Index n3 = n2 + 1;
      triangles_.push_back({NodeIndex{n0}, NodeIndex{n1}, NodeIndex{n2}});
      triangles_.push_back({NodeIndex{n1}, NodeIndex{n3}, NodeIndex{n2}});
    }
  }
}
}  // namespace aether::mesh
