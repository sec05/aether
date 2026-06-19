// Copyright 2026 Spencer Evans-Cole

#include <vector>
#include <algorithm>

#include "aether/mesh/quad.hpp"

namespace aether::mesh {

namespace {
// Bilinear map: reference [0,1]^2 -> physical quad.
// (xi,eta)=(0,0)->c0, (1,0)->c1, (1,1)->c2, (0,1)->c3.
Vec2 bilinear(const std::array<Vec2, 4>& c, Real xi, Real eta) {
  const Real w0 = (1 - xi) * (1 - eta);
  const Real w1 = xi * (1 - eta);
  const Real w2 = xi * eta;
  const Real w3 = (1 - xi) * eta;
  return Vec2(w0 * c[0].x() + w1 * c[1].x() + w2 * c[2].x() + w3 * c[3].x(),
              w0 * c[0].y() + w1 * c[1].y() + w2 * c[2].y() + w3 * c[3].y());
}
}  // namespace

Quad::Quad(int degree, int nx, int ny, std::array<Vec2, 4> corners)
    : Mesh(degree, nx * ny), corners_(corners), nx_(nx), ny_(ny) {  // NOLINT
  const int total = nx * ny;
  num_nodes_ = total;

  nodes_.reserve(total);
  for (int j = 0; j < ny; ++j) {
    const Real eta = static_cast<Real>(j) / (ny - 1);
    for (int i = 0; i < nx; ++i) {
      const Real xi = static_cast<Real>(i) / (nx - 1);
      nodes_.push_back(bilinear(corners_, xi, eta));
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

std::array<NodeIndex, 3> Quad::element_node_indices(int element_index) const {
  const auto& tri = triangles_[element_index];
  return {tri[0], tri[1], tri[2]};
}

std::vector<Vec2> Quad::element_nodes(int element_index) const {
  const auto idx = element_node_indices(element_index);
  return {nodes_[idx[0].get()], nodes_[idx[1].get()], nodes_[idx[2].get()]};
}

// Topological boundary: ring of the structured index grid. No geometry,
// no float comparisons — correct for any quad shape.
std::vector<NodeIndex> Quad::boundary_nodes() const {
  std::vector<NodeIndex> boundary;
  for (int j = 0; j < ny_; ++j)
    for (int i = 0; i < nx_; ++i)
      if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1)
        boundary.push_back(NodeIndex{static_cast<Index>(j * nx_ + i)});
  return boundary;
}
}  // namespace aether::mesh
