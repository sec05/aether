// Copyright 2026 Spencer Evans-Cole

#include <vector>
#include <algorithm>

#include "aether/mesh/rectangle.hpp"

namespace aether::mesh {
Rectangle::Rectangle(int degree, int num_nodes, std::vector<Vec2> domain)
    : Mesh(degree, num_nodes), domain_(domain) {  // NOLINT
  const Real lx = domain_[0][1] - domain_[0][0];
  const Real ly = domain_[1][1] - domain_[1][0];
  const int nx = num_nodes;                                  // nodes along x
  const Real h = lx / (nx - 1);                              // target spacing
  const int ny = static_cast<int>(std::lround(ly / h)) + 1;  // nodes along y
  const int total = nx * ny;
  num_nodes_ = total;

  nodes_.reserve(total);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      nodes_.emplace_back(domain_[0][0] + i * h, domain_[1][0] + j * h);
    }
  }

  triangles_.reserve(2 * (nx - 1) * (ny - 1));
  for (Index j = 0; j < ny - 1; ++j) {
    for (Index i = 0; i < nx - 1; ++i) {
      const Index n0 = j * nx + i;
      const Index n1 = n0 + 1;
      const Index n2 = n0 + nx;
      const Index n3 = n2 + 1;
      triangles_.push_back({NodeIndex{n0}, NodeIndex{n1}, NodeIndex{n2}});
      triangles_.push_back({NodeIndex{n1}, NodeIndex{n3}, NodeIndex{n2}});
    }
  }

  num_elements_ = static_cast<int>(triangles_.size());
}

std::array<NodeIndex, 3> Rectangle::element_node_indices(int element_index) const {
  const auto& tri = triangles_[element_index];
  return {tri[0], tri[1], tri[2]};
}

std::vector<Vec2> Rectangle::element_nodes(int element_index) const {
  const auto idx = element_node_indices(element_index);
  return {nodes_[idx[0].get()], nodes_[idx[1].get()], nodes_[idx[2].get()]};
}

std::vector<NodeIndex> Rectangle::boundary_nodes() const {
  const Real x0 = domain_[0][0], x1 = domain_[0][1];
  const Real y0 = domain_[1][0], y1 = domain_[1][1];
  const Real eps = 1e-12 * std::max({std::abs(x1 - x0), std::abs(y1 - y0), Real(1)});

  std::vector<NodeIndex> boundary;
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    const auto& node = nodes_[i];
    if (std::abs(node.x() - x0) < eps || std::abs(node.x() - x1) < eps ||
        std::abs(node.y() - y0) < eps || std::abs(node.y() - y1) < eps) {
      boundary.push_back(NodeIndex{static_cast<Index>(i)});
    }
  }
  return boundary;
}
}  // namespace aether::mesh
