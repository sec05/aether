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

struct BoundaryFacet {
  std::array<NodeIndex, 2> nodes;  // P1 edge endpoints, triangle-CCW order
  int marker;                      // 1..4 = outer sides, >=5 = hole rims
};

class Mesh {
 public:
  Mesh(int degree, int num_nodes) : degree_(degree), num_nodes_(num_nodes) {}
  virtual ~Mesh() = default;
  int degree() const { return degree_; }
  int num_nodes() const { return num_nodes_; }
  int num_elements() const { return num_elements_; }
  const std::vector<Vec2>& nodes() const { return nodes_; }
  virtual std::array<NodeIndex, 3> element_node_indices(int element_index) const = 0;
  virtual std::vector<Vec2> element_nodes(int element_index) const = 0;
  virtual std::vector<NodeIndex> boundary_nodes() const = 0;
  virtual Vec2 node(int node_index) const { return nodes_[node_index]; }
  virtual const std::vector<BoundaryFacet>& boundary_facets() const { return boundary_facets_; }

 protected:
  int degree_;
  int num_nodes_;
  std::vector<Vec2> nodes_;
  int num_elements_;
  std::vector<BoundaryFacet> boundary_facets_;
};
}  // namespace aether::mesh
