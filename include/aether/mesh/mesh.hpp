// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <variant>
#include <vector>

#include "../core/geometry.hpp"
#include "../core/types.hpp"

// Abstract mesh class. Meshes are responsible for storing the geometry and topology of the domain,
// and providing access to this information for the discretization and assembly layers. Meshes are
// immutable after construction. We assume the domain is [0,1] x [0,1] (2D) or [0,1] x [0,1] x [0,1]
// (3D) for simplicity. Meshes are not responsible for storing solution fields.

namespace aether::mesh {
class Mesh {
 public:
  Mesh(int degree, int num_nodes) : degree_(degree), num_nodes_(num_nodes) {}
  virtual ~Mesh() = default;
  int degree() const { return degree_; }
  int num_nodes() const { return num_nodes_; }
  const std::vector<Vec2> &nodes() const { return nodes_; }

 protected:
  int degree_;
  int num_nodes_;
  std::vector<Vec2> nodes_;
};
}  // namespace aether::mesh
