// Copyright 2026 Spencer Evans-Cole

#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "../core/types.hpp"

namespace aether::solver {
class Solver {
 public:
  Eigen::VectorXd solve(const Eigen::SparseMatrix<Real>& A, const Eigen::VectorXd& b);
};
}  // namespace aether::solver
