// Copyright 2026 Spencer Evans-Cole

#include "aether/assembly/assembler.hpp"

#include <array>
#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace aether::assembly {

void Assembler::assemble() {
  const int num_elements = mesh_.num_elements();
  const int n = mesh_.num_nodes();  // # DOFs for P1 (DOF == node)

  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.reserve(static_cast<std::size_t>(9 * num_elements));

  for (int e = 0; e < num_elements; ++e) {
    Mat3 local = local_stiffness_matrix(e);
    scatter(e, local, &triplets);
  }

  global_stiffness_matrix_.resize(n, n);
  global_stiffness_matrix_.setFromTriplets(triplets.begin(), triplets.end());  // sums duplicates

  rhs_.resize(n);
  rhs_.setZero();
}

Mat3 Assembler::local_stiffness_matrix(int element_index) const {
  const std::vector<Vec2> nodes = mesh_.element_nodes(element_index);

  Mat2 J;
  J.col(0) = nodes[1] - nodes[0];  // columns, not brace-init
  J.col(1) = nodes[2] - nodes[0];

  const Real detJ = J.determinant();
  const Mat2 JinvT = J.inverse().transpose();  // J^{-T}, not J^{-1}
  const Real absdet = std::abs(detJ);

  Mat3 local = Mat3::Zero();
  for (const auto& [xi, w] : ref_element_.quadrature_points(1)) {
    const std::array<Vec2, 3> gref = ref_element_.shape_grad(xi);  // call once

    std::array<Vec2, 3> grad;
    for (int i = 0; i < 3; ++i) grad[i] = JinvT * gref[i];  // physical gradients

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) local(i, j) += w * grad[i].dot(grad[j]) * absdet;
  }
  return local;
}

void Assembler::scatter(int element_index, const Mat3& local,
                        std::vector<Eigen::Triplet<Real>>* triplets) const {
  const auto idx = mesh_.element_node_indices(element_index);
  for (int a = 0; a < 3; ++a) {
    const int I = idx[a].get();  // TypedIndex -> raw row
    for (int b = 0; b < 3; ++b) {
      const int J = idx[b].get();
      triplets->emplace_back(I, J, local(a, b));
    }
  }
}

void Assembler::homogeneous_dirichlet_bc(const std::vector<NodeIndex>& boundary_nodes) {
  for (const auto& node : boundary_nodes) {
    const int i = node.get();
    global_stiffness_matrix_.coeffRef(i, i) = 1.0;  // K_ii = 1
    for (Eigen::SparseMatrix<Real>::InnerIterator it(global_stiffness_matrix_, i); it; ++it) {
      if (it.row() != i) it.valueRef() = 0.0;  // K_ij = 0 for j != i
    }
    rhs_(i) = 0.0;  // f_i = 0
  }
}
}  // namespace aether::assembly
