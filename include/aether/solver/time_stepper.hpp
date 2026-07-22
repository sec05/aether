// Copyright 2026 Spencer Evans-Cole
// aether/solver/time_stepper.hpp
#pragma once

#include <functional>
#include <vector>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include "aether/core/geometry.hpp"  // Vec2, Real, BoundaryFn, BoundaryCondition, BCType
#include "aether/mesh/mesh.hpp"

namespace aether::solver {

// Transient heat equation  M u' + K u = b(t)  via the theta-method:
//   (M + theta*dt*K) u^{n+1}
//       = (M - (1-theta)*dt*K) u^n + dt[(1-theta) b^n + theta b^{n+1}].
// theta = 1.0 -> backward Euler (1st order, unconditionally stable, strong damping).
// theta = 0.5 -> Crank-Nicolson (2nd order, unconditionally stable, can ring).
class TimeStepper {
 public:
  using LoadFn = std::function<Eigen::VectorX<Real>(Real)>;
  using InitFn = std::function<Real(const Vec2&)>;

  TimeStepper(const mesh::Mesh& mesh, const Eigen::SparseMatrix<Real>& stiffness,
              const Eigen::SparseMatrix<Real>& mass, const std::vector<BoundaryCondition>& bcs,
              Real theta, Real dt, Real alpha);

  void set_initial_condition(const InitFn& u0);  // call before stepping
  void set_load(LoadFn load);                    // optional; default zero

  void step();                  // advance one dt
  void advance(int num_steps);  // advance num_steps

  const Eigen::VectorX<Real>& solution() const { return u_; }
  Real time() const { return t_; }

 private:
  struct DirichletDof {
    int dof;
    Vec2 coord;
    BoundaryFn g;
  };

  const mesh::Mesh& mesh_;
  Eigen::SparseMatrix<Real> K_;                           // pristine spatial operator (copy)
  Eigen::SparseMatrix<Real> M_;                           // pristine mass (copy)
  Eigen::SparseMatrix<Real> A_;                           // M + theta*dt*K, pristine (for the lift)
  Eigen::SimplicialLLT<Eigen::SparseMatrix<Real>> chol_;  // factored A w/ BCs

  std::vector<DirichletDof> dirichlet_;
  LoadFn load_;
  Eigen::VectorX<Real> load_prev_;  // cached b^n
  bool started_ = false;

  Real theta_;
  Real dt_;
  Real t_ = Real(0);
  Real alpha_ = Real(1);
  Eigen::VectorX<Real> u_;
};

}  // namespace aether::solver
