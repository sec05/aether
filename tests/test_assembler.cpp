// Copyright 2026 Spencer Evans-Cole
//
// Tests for the P1 assembler + solve pipeline on the unit-square mesh
// produced by aether::mesh::Square(1, 16) (a 4x4 node grid, 18 triangles).
//
// These encode the Laplace problem your current assembler produces:
//   -Δu = 0  on (0,1)^2,  u = g on ∂Ω,  with g(x,y) = x^2 - y^2.
// Because g is harmonic, the exact solution is u = x^2 - y^2, and the P1
// nodal values reproduce it exactly on this structured mesh. If you change
// the source term or boundary data, update kExactSolution / the RHS test.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "aether/core/types.hpp"
#include "aether/mesh/square.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/assembly/assembler.hpp"
#include "aether/solver/solver.hpp"

namespace {

// Harmonic exact solution for the configured boundary data.
double ExactSolution(double x, double y) { return x * x - y * y; }

class AssemblerTest : public ::testing::Test {
 protected:
  AssemblerTest() : mesh_(1, 16), assembler_(mesh_, ref_) {}

  void SetUp() override { assembler_.assemble(); }

  // Dense copy of the (small) stiffness matrix for convenient inspection.
  Eigen::MatrixXd Dense() const { return Eigen::MatrixXd(*assembler_.stiffness_matrix()); }

  bool IsBoundary(int i) const {
    constexpr double kTol = 1e-9;
    const double x = mesh_.node(i).x();
    const double y = mesh_.node(i).y();
    return std::abs(x) < kTol || std::abs(x - 1.0) < kTol || std::abs(y) < kTol ||
           std::abs(y - 1.0) < kTol;
  }

  aether::mesh::Square mesh_;
  aether::elements::P1Element ref_;
  aether::assembly::Assembler assembler_;
};

TEST_F(AssemblerTest, MeshTopology) {
  EXPECT_EQ(mesh_.num_nodes(), 16);
  EXPECT_EQ(mesh_.num_elements(), 18);
}

TEST_F(AssemblerTest, StiffnessDimensions) {
  const auto& K = *assembler_.stiffness_matrix();
  EXPECT_EQ(K.rows(), mesh_.num_nodes());
  EXPECT_EQ(K.cols(), mesh_.num_nodes());
}

// The assembled pattern stores 24 entries: 16 diagonals + 2 * 4 edges
// (4 horizontal + 4 vertical + 4 hypotenuse edges). nonZeros() counts
// stored coefficients, including the structural zeros from the hypotenuse
// couplings and the Dirichlet elimination.
TEST_F(AssemblerTest, StoredPatternHas24Entries) {
  EXPECT_EQ(assembler_.stiffness_matrix()->nonZeros(), 24);
}

TEST_F(AssemblerTest, StiffnessIsSymmetric) {
  const Eigen::MatrixXd K = Dense();
  EXPECT_LT((K - K.transpose()).norm(), 1e-12);
}

TEST_F(AssemblerTest, StiffnessIsPositiveDefinite) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<aether::Real>> llt(*assembler_.stiffness_matrix());
  EXPECT_EQ(llt.info(), Eigen::Success);
}

// Every interior node carries the standard 5-point Laplacian diagonal (4),
// independent of node numbering.
TEST_F(AssemblerTest, InteriorNodesHaveStandardDiagonal) {
  const Eigen::MatrixXd K = Dense();
  std::vector<int> interior;
  for (int i = 0; i < mesh_.num_nodes(); ++i) {
    if (!IsBoundary(i)) interior.push_back(i);
  }
  EXPECT_EQ(interior.size(), 4u);
  for (int i : interior) {
    EXPECT_NEAR(K(i, i), 4.0, 1e-9);
  }
}

// Explicit regression check of the known interior block for Square(1, 16).
// Interior nodes are 5, 6, 9, 10; edge neighbors couple at -1 and the
// hypotenuse (diagonal) neighbors cancel to exactly 0.
TEST_F(AssemblerTest, KnownInteriorStencilValues) {
  const Eigen::MatrixXd K = Dense();
  EXPECT_NEAR(K(5, 5), 4.0, 1e-9);
  EXPECT_NEAR(K(6, 6), 4.0, 1e-9);
  EXPECT_NEAR(K(9, 9), 4.0, 1e-9);
  EXPECT_NEAR(K(10, 10), 4.0, 1e-9);

  EXPECT_NEAR(K(5, 6), -1.0, 1e-9);
  EXPECT_NEAR(K(5, 9), -1.0, 1e-9);
  EXPECT_NEAR(K(6, 10), -1.0, 1e-9);
  EXPECT_NEAR(K(9, 10), -1.0, 1e-9);

  // Diagonal-neighbor couplings cancel for right-isosceles P1 elements.
  EXPECT_NEAR(K(5, 10), 0.0, 1e-9);
  EXPECT_NEAR(K(6, 9), 0.0, 1e-9);
}

// Dirichlet rows are reduced to identity, and boundary RHS entries hold the
// prescribed data g(x,y) = x^2 - y^2.
TEST_F(AssemblerTest, BoundaryRowsAndRhs) {
  const Eigen::MatrixXd K = Dense();
  const Eigen::VectorXd rhs = *assembler_.rhs();
  ASSERT_EQ(rhs.size(), mesh_.num_nodes());

  for (int i = 0; i < mesh_.num_nodes(); ++i) {
    if (!IsBoundary(i)) continue;
    EXPECT_NEAR(K(i, i), 1.0, 1e-9) << "boundary node " << i;
    const double off_diag = K.row(i).cwiseAbs().sum() - std::abs(K(i, i));
    EXPECT_NEAR(off_diag, 0.0, 1e-9) << "boundary node " << i;

    const double x = mesh_.node(i).x();
    const double y = mesh_.node(i).y();
    EXPECT_NEAR(rhs(i), ExactSolution(x, y), 1e-9) << "boundary node " << i;
  }
}

// End-to-end: the solve reproduces the harmonic exact solution at the nodes.
TEST_F(AssemblerTest, SolveMatchesExactSolution) {
  aether::solver::Solver solver;
  const Eigen::VectorXd solution = solver.solve(*assembler_.stiffness_matrix(), *assembler_.rhs());
  ASSERT_EQ(solution.size(), mesh_.num_nodes());

  Eigen::VectorXd u_exact(mesh_.num_nodes());
  for (int i = 0; i < mesh_.num_nodes(); ++i) {
    u_exact(i) = ExactSolution(mesh_.node(i).x(), mesh_.node(i).y());
  }

  const double err_inf = (solution - u_exact).cwiseAbs().maxCoeff();
  EXPECT_LT(err_inf, 1e-9);
}

}  // namespace
