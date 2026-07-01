// Copyright 2026 Spencer Evans-Cole

#include "aether/assembly/assembler.hpp"

#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <utility>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "aether/core/geometry.hpp"
#include "aether/mesh/quad.hpp"

namespace aether::assembly {

void Assembler::assemble_stiffness() {
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

void Assembler::apply_boundary_conditions(const std::vector<BoundaryCondition>& conditions,
                                          Real t) {
  // Index conditions by marker. Last write wins if a marker is listed twice.
  std::unordered_map<int, const BoundaryCondition*> by_marker;
  for (const auto& bc : conditions) by_marker[bc.marker] = &bc;

  // BoundaryFn is double-domain; pin the conversions to Real explicitly here so
  // nothing narrows silently if Real ever becomes float.
  const double td = static_cast<double>(t);
  auto X = [&](NodeIndex i) -> Vec2 { return mesh_.nodes()[i.get()]; };

  // 2-point Gauss on the reference edge [0,1]. sqrt isn't constexpr before
  // C++23, so these are static const rather than constexpr.
  static const std::array<Real, 2> qs{Real(0.5) - Real(0.5) / std::sqrt(Real(3)),
                                      Real(0.5) + Real(0.5) / std::sqrt(Real(3))};
  static const std::array<Real, 2> qw{Real(0.5), Real(0.5)};

  std::vector<std::pair<NodeIndex, Real>> dirichlet;  // deferred to the end

  for (const auto& f : mesh_.boundary_facets()) {
    auto it = by_marker.find(f.marker);
    if (it == by_marker.end()) continue;  // unassigned marker => homogeneous Neumann
    const BoundaryCondition& bc = *it->second;

    const NodeIndex a = f.nodes[0];
    const NodeIndex b = f.nodes[1];

    // --- Essential: collect endpoint values, apply after the loop. ---
    if (bc.type == BCType::Dirichlet) {
      dirichlet.emplace_back(a, static_cast<Real>(bc.value(X(a), td)));
      dirichlet.emplace_back(b, static_cast<Real>(bc.value(X(b), td)));
      continue;
    }

    // --- Natural: integrate over this edge. ---
    const Vec2 xa = X(a);
    const Vec2 xb = X(b);
    const Real L = (xb - xa).norm();

    for (int k = 0; k < 2; ++k) {
      const Real s = qs[k];
      const Real N0 = Real(1) - s;
      const Real N1 = s;
      const Real jxw = qw[k] * L;        // weight * edge Jacobian
      const Vec2 x = N0 * xa + N1 * xb;  // physical quadrature point

      // Neumann flux q, or Robin's g — both land on the RHS as +∫ g φ ds.
      const Real g = static_cast<Real>(bc.value(x, td));
      rhs_(a.get()) += jxw * g * N0;
      rhs_(b.get()) += jxw * g * N1;

      // Robin adds the boundary-mass block +∫ α φ_i φ_j ds to the matrix.
      if (bc.type == BCType::Robin) {
        const Real alpha = static_cast<Real>(bc.robin_coeff(x, td));
        const Real c = jxw * alpha;
        auto& K = global_stiffness_matrix_;  // entries already in the pattern
        K.coeffRef(a.get(), a.get()) += c * N0 * N0;
        K.coeffRef(a.get(), b.get()) += c * N0 * N1;
        K.coeffRef(b.get(), a.get()) += c * N1 * N0;
        K.coeffRef(b.get(), b.get()) += c * N1 * N1;
      }
    }
  }

  // Apply Dirichlet last, so it wins on any corner shared with a natural edge.
  if (!dirichlet.empty()) dirichlet_bc(dirichlet);
}

void Assembler::assemble_load(const BoundaryFn& f, Real t) {
  const int num_elements = mesh_.num_elements();
  const double td = static_cast<double>(t);

  for (int e = 0; e < num_elements; ++e) {
    const std::vector<Vec2> nodes = mesh_.element_nodes(e);
    const auto idx = mesh_.element_node_indices(e);

    // Same affine geometry as the stiffness routine.
    Mat2 J;
    J.col(0) = nodes[1] - nodes[0];
    J.col(1) = nodes[2] - nodes[0];
    const Real absdet = std::abs(J.determinant());

    // Higher order than stiffness on purpose: P1 gradients are constant so
    // 1-point nails stiffness, but f·φ is non-constant, so under-integrating
    // here would inject quadrature error into your convergence rate.
    for (const auto& [xi, w] : ref_element_.quadrature_points(2)) {
      const std::array<Real, 3> N = ref_element_.shape(xi);

      // Map the reference quadrature point to physical space.
      Vec2 x = Vec2::Zero();
      for (int i = 0; i < 3; ++i) x += N[i] * nodes[i];

      const Real fx = static_cast<Real>(f(x, td));
      const Real jxw = w * fx * absdet;
      for (int i = 0; i < 3; ++i) rhs_(idx[i].get()) += jxw * N[i];
    }
  }
}

Mat3 Assembler::local_mass_matrix(int element_index) const {
  const std::vector<Vec2> nodes = mesh_.element_nodes(element_index);
  Mat2 J;
  J.col(0) = nodes[1] - nodes[0];
  J.col(1) = nodes[2] - nodes[0];
  const Real absdet = std::abs(J.determinant());
  Mat3 m;
  m << Real(2), Real(1), Real(1), Real(1), Real(2), Real(1), Real(1), Real(1), Real(2);
  return (absdet / Real(24)) * m;  // (absdet/2)/12 = absdet/24
}

void Assembler::assemble_mass() {
  const int num_elements = mesh_.num_elements();
  const int n = mesh_.num_nodes();
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.reserve(static_cast<std::size_t>(9 * num_elements));
  for (int e = 0; e < num_elements; ++e) {
    Mat3 local = local_mass_matrix(e);
    scatter(e, local, &triplets);  // reuses your existing scatter unchanged
  }
  global_mass_matrix_.resize(n, n);
  global_mass_matrix_.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace aether::assembly
