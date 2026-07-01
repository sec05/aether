// Copyright 2026 Spencer Evans-Cole
// aether/solver/time_stepper.cpp
#include "aether/solver/time_stepper.hpp"

#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace aether::solver {

TimeStepper::TimeStepper(const mesh::Mesh& mesh, const Eigen::SparseMatrix<Real>& stiffness,
                         const Eigen::SparseMatrix<Real>& mass,
                         const std::vector<BoundaryCondition>& bcs, Real theta, Real dt)
    : mesh_(mesh), K_(stiffness), M_(mass), theta_(theta), dt_(dt) {  // NOLINT
  const int n = mesh_.num_nodes();

  // Dirichlet set: endpoints of facets whose marker carries a Dirichlet BC.
  // Last facet wins on a shared corner node.
  std::unordered_map<int, const BoundaryCondition*> by_marker;
  for (const auto& bc : bcs)
    if (bc.type == BCType::Dirichlet) by_marker[bc.marker] = &bc;

  std::unordered_map<int, BoundaryFn> dof_to_g;
  for (const auto& f : mesh_.boundary_facets()) {
    auto it = by_marker.find(f.marker);
    if (it == by_marker.end()) continue;
    dof_to_g[f.nodes[0].get()] = it->second->value;
    dof_to_g[f.nodes[1].get()] = it->second->value;
  }
  dirichlet_.reserve(dof_to_g.size());
  for (const auto& [dof, g] : dof_to_g) dirichlet_.push_back({dof, mesh_.node(dof), g});

  // Constant system matrix; keep A_ pristine for the per-step lift.
  A_ = M_ + theta_ * dt_ * K_;

  // Symmetric elimination of the Dirichlet rows/cols on a copy, then factor.
  std::vector<char> fixed(n, 0);
  for (const auto& d : dirichlet_) fixed[d.dof] = 1;

  Eigen::SparseMatrix<Real> A_bc = A_;
  for (int c = 0; c < A_bc.outerSize(); ++c)
    for (Eigen::SparseMatrix<Real>::InnerIterator it(A_bc, c); it; ++it)
      if (fixed[it.row()] || fixed[it.col()])
        it.valueRef() = (it.row() == it.col()) ? Real(1) : Real(0);
  A_bc.prune(Real(0));

  chol_.compute(A_bc);  // A_bc is SPD: M SPD + theta*dt*K PSD, identity on fixed
  if (chol_.info() != Eigen::Success)
    throw std::runtime_error("TimeStepper: Cholesky factorization failed");

  // Defaults: zero load, zero field. set_initial_condition overrides u_.
  load_ = [n](Real) { return Eigen::VectorX<Real>::Zero(n); };
  u_ = Eigen::VectorX<Real>::Zero(n);
}

void TimeStepper::set_initial_condition(const InitFn& u0) {
  const int n = mesh_.num_nodes();
  u_.resize(n);
  for (int i = 0; i < n; ++i) u_[i] = u0(mesh_.node(i));
  // Enforce Dirichlet data at t=0 so the IC is consistent with the BCs.
  for (const auto& d : dirichlet_) u_[d.dof] = static_cast<Real>(d.g(d.coord, 0.0));
  t_ = Real(0);
  started_ = false;
}

void TimeStepper::set_load(LoadFn load) {
  load_ = std::move(load);
  started_ = false;  // re-cache b^n on the next step
}

void TimeStepper::step() {
  const int n = mesh_.num_nodes();
  const Real t_new = t_ + dt_;

  if (!started_) {
    load_prev_ = load_(t_);
    started_ = true;
  }
  const Eigen::VectorX<Real>& b_old = load_prev_;   // b^n  (cached)
  const Eigen::VectorX<Real> b_new = load_(t_new);  // b^{n+1}

  // theta-method RHS from the previous solution.
  Eigen::VectorX<Real> r = M_ * u_ - (Real(1) - theta_) * dt_ * (K_ * u_) +
                           dt_ * ((Real(1) - theta_) * b_old + theta_ * b_new);

  // Dirichlet lift + overwrite at t_new, against the PRISTINE A_.
  Eigen::VectorX<Real> g = Eigen::VectorX<Real>::Zero(n);
  for (const auto& d : dirichlet_)
    g[d.dof] = static_cast<Real>(d.g(d.coord, static_cast<double>(t_new)));
  r -= A_ * g;
  for (const auto& d : dirichlet_) r[d.dof] = g[d.dof];

  u_ = chol_.solve(r);
  t_ = t_new;
  load_prev_ = std::move(b_new);  // becomes b^n next step
}

void TimeStepper::advance(int num_steps) {
  for (int s = 0; s < num_steps; ++s) step();
}

}  // namespace aether::solver
