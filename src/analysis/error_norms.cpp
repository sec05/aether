// Copyright 2026 Spencer Evans-Cole
#include "aether/analysis/error_norms.hpp"
#include <array>
#include <cmath>
#include <algorithm>
#include <vector>

namespace aether::analysis {

ErrorNorms compute_errors(const mesh::Mesh& mesh, const elements::ReferenceElement& ref,
                          const Eigen::VectorX<Real>& u_h, const BoundaryFn& u_exact, Real t) {
  const double td = static_cast<double>(t);
  Real l2_sq = Real(0);
  Real linf = Real(0);

  const int num_elements = mesh.num_elements();
  for (int e = 0; e < num_elements; ++e) {
    const std::vector<Vec2> nodes = mesh.element_nodes(e);
    const auto idx = mesh.element_node_indices(e);

    Mat2 J;
    J.col(0) = nodes[1] - nodes[0];
    J.col(1) = nodes[2] - nodes[0];
    const Real absdet = std::abs(J.determinant());

    for (const auto& [xi, w] : ref.quadrature_points(2)) {
      const std::array<Real, 3> N = ref.shape(xi);

      Vec2 x = Vec2::Zero();
      Real uh = Real(0);
      for (int i = 0; i < 3; ++i) {
        x += N[i] * nodes[i];
        uh += N[i] * u_h(idx[i].get());
      }
      const Real diff = std::abs(uh - static_cast<Real>(u_exact(x, td)));
      l2_sq += w * diff * diff * absdet;
      linf = std::max(linf, diff);  // sup sampled at quadrature points
    }
  }
  return {std::sqrt(l2_sq), linf};
}

}  // namespace aether::analysis
