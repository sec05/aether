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
  // x^2 - y^2 test solution
  std::vector<std::pair<NodeIndex, Real>> bcs;
  for (const auto& node : mesh_.boundary_nodes()) {
    Vec2 coord = mesh_.nodes()[node.get()];
    Real value = std::pow(coord.x(), 2) - std::pow(coord.y(), 2);
    bcs.emplace_back(node, value);
  }
  dirichlet_bc(bcs);
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

void Assembler::dirichlet_bc(const std::vector<std::pair<NodeIndex, Real>>& bcs) {
  const int n = global_stiffness_matrix_.rows();
  std::vector<char> fixed(n, 0);
  Eigen::VectorX<Real> g = Eigen::VectorX<Real>::Zero(n);
  for (const auto& [node, value] : bcs) {
    fixed[node.get()] = 1;
    g[node.get()] = value;
  }

  // 1) Lift constrained columns to the RHS, using the ORIGINAL matrix.
  //    (g is zero on free dofs, so this only moves the K_FD * g_D terms.)
  rhs_ -= global_stiffness_matrix_ * g;

  // 2) Prescribe the value on each constrained row.
  for (const auto& [node, value] : bcs) rhs_(node.get()) = value;

  // 3) Symmetric elimination in one storage-order-agnostic sweep:
  //    zero any entry whose row OR column is constrained; unit diagonal.
  for (int c = 0; c < global_stiffness_matrix_.outerSize(); ++c)
    for (Eigen::SparseMatrix<Real>::InnerIterator it(global_stiffness_matrix_, c); it; ++it)
      if (fixed[it.row()] || fixed[it.col()])
        it.valueRef() = (it.row() == it.col()) ? Real(1) : Real(0);
  global_stiffness_matrix_.prune(aether::Real(0));
}
}  // namespace aether::assembly
