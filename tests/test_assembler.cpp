// Copyright 2026 Spencer Evans-Cole

#include <algorithm>

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

#include "aether/assembly/assembler.hpp"
#include "aether/mesh/square.hpp"
#include "aether/elements/p1_element.hpp"

namespace aether {
namespace {

constexpr int kDegree = 1;
constexpr int kNumNodes = 9;

class AssemblerTest : public ::testing::Test {
 protected:
  void SetUp() override {
    assembler_.assemble();

    const auto& K = *assembler_.stiffness_matrix();  // const SparseMatrix<Real>&
    n_ = mesh_.num_nodes();
    Kd_ = K.cast<double>();  // works for Real=float|double
    nnz_ = K.nonZeros();
    rows_ = K.rows();
    cols_ = K.cols();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Kd_);
    eig_ = es.eigenvalues();                       // ascending
    tol_ = 1e-9 * std::max(eig_.maxCoeff(), 1.0);  // scaled by spectral norm
  }

  mesh::Square mesh_{kDegree, kNumNodes};
  elements::P1Element ref_;
  assembly::Assembler assembler_{mesh_, ref_};

  int n_ = 0;
  Eigen::Index nnz_ = 0, rows_ = 0, cols_ = 0;
  Eigen::MatrixXd Kd_;
  Eigen::VectorXd eig_;
  double tol_ = 0.0;
};

// K must be num_nodes x num_nodes.
TEST_F(AssemblerTest, Dimensions) {
  EXPECT_EQ(rows_, n_);
  EXPECT_EQ(cols_, n_);
}

// assemble() must actually populate the matrix.
TEST_F(AssemblerTest, NonEmpty) { EXPECT_GT(nnz_, 0); }

// a(u,v) = integral grad(u).grad(v) is symmetric; local matrices are
// symmetric and scatter is index-symmetric, so K = K^T to machine eps.
TEST_F(AssemblerTest, Symmetry) {
  const double err = (Kd_ - Kd_.transpose()).cwiseAbs().maxCoeff();
  EXPECT_LT(err, 1e-12) << "max|K - K^T| = " << err;
}

// Partition of unity + grad(const)=0  =>  K * 1 = 0. The strongest single
// assembly check: nearly every scatter/index mistake gives a nonzero row sum.
TEST_F(AssemblerTest, ConstantNullSpace) {
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n_);
  const double err = (Kd_ * ones).cwiseAbs().maxCoeff();
  EXPECT_LT(err, 1e-10) << "max|K * 1| = " << err;
}

// K_ii = integral |grad(phi_i)|^2 > 0 for every node.
TEST_F(AssemblerTest, PositiveDiagonal) {
  const double min_diag = Kd_.diagonal().minCoeff();
  EXPECT_GT(min_diag, 0.0) << "min diagonal entry = " << min_diag;
}

// No BCs + connected mesh => K is positive semi-definite.
TEST_F(AssemblerTest, PositiveSemiDefinite) {
  EXPECT_GT(eig_.minCoeff(), -tol_) << "min eigenvalue = " << eig_.minCoeff();
}

// The constant is the only null vector: lambda_0 ~ 0, lambda_1 > 0. A second
// zero eigenvalue means a disconnected component, duplicated node, or dropped
// element.
TEST_F(AssemblerTest, NullSpaceIsOneDimensional) {
  EXPECT_NEAR(eig_(0), 0.0, tol_) << "lambda_0 = " << eig_(0);
  EXPECT_GT(eig_(1), tol_) << "lambda_1 = " << eig_(1);
}

}  // namespace
}  // namespace aether
