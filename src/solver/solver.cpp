// Copyright 2026 Spencer Evans-Cole

#include "aether/solver/solver.hpp"

namespace aether::solver {
Eigen::VectorXd Solver::solve(const Eigen::SparseMatrix<Real>& A, const Eigen::VectorXd& b) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  x = solver.solve(b);
  return x;
}
}  // namespace aether::solver
