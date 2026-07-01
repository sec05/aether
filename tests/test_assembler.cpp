// Copyright 2026 Spencer Evans-Cole
//
// Tests for the P1 assembler on the unit square (Rectangle(1, 4) => 4x4 node
// grid, 18 triangles). The assembler separates concerns:
//   * assemble_stiffness()        -> pristine stiffness K, zeroed rhs
//   * assemble_mass()             -> consistent mass M
//   * assemble_load(f)            -> interior load into rhs
//   * apply_boundary_conditions() -> natural terms + symmetric Dirichlet
//   * dirichlet_bc(node->value)   -> the symmetric-elimination primitive
// So operator tests inspect the PRISTINE matrix; BC tests run elimination first.
//
// Harmonic MMS: -Δu = 0, u = g = x^2 - y^2 on ∂Ω. Because g is harmonic the
// exact solution is u = x^2 - y^2 and the P1 nodal values reproduce it exactly
// on this structured mesh.

#include <cmath>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "aether/assembly/assembler.hpp"
#include "aether/core/types.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/mesh/rectangle.hpp"
#include "aether/solver/solver.hpp"

namespace {

double ExactSolution(double x, double y) { return x * x - y * y; }

class AssemblerTest : public ::testing::Test {
 protected:
  AssemblerTest() : mesh_(1, 4, {{0.0, 1.0}, {0.0, 1.0}}), assembler_(mesh_, ref_) {}

  void SetUp() override { assembler_.assemble_stiffness(); }  // pristine K, zeroed rhs

  Eigen::MatrixXd Dense() const { return Eigen::MatrixXd(*assembler_.stiffness_matrix()); }
  Eigen::MatrixXd DenseMass() const { return Eigen::MatrixXd(*assembler_.mass_matrix()); }

  bool IsBoundary(int i) const {
    constexpr double kTol = 1e-9;
    const double x = mesh_.node(i).x();
    const double y = mesh_.node(i).y();
    return std::abs(x) < kTol || std::abs(x - 1.0) < kTol || std::abs(y) < kTol ||
           std::abs(y - 1.0) < kTol;
  }

  // Rectangle exposes no facet markers, so apply the harmonic data as a node-set
  // Dirichlet map through the elimination primitive directly.
  void ApplyHarmonicDirichlet() {
    std::vector<std::pair<aether::NodeIndex, aether::Real>> bc;
    for (int i = 0; i < mesh_.num_nodes(); ++i) {
      if (!IsBoundary(i)) continue;
      bc.emplace_back(aether::NodeIndex(i), ExactSolution(mesh_.node(i).x(), mesh_.node(i).y()));
    }
    assembler_.dirichlet_bc(bc);
  }

  aether::mesh::Rectangle mesh_;
  aether::elements::P1Element ref_;
  aether::assembly::Assembler assembler_;
};

// ---- Mesh / dimensions ----------------------------------------------------

TEST_F(AssemblerTest, MeshTopology) {
  EXPECT_EQ(mesh_.num_nodes(), 16);
  EXPECT_EQ(mesh_.num_elements(), 18);
}

TEST_F(AssemblerTest, StiffnessDimensions) {
  const auto& K = *assembler_.stiffness_matrix();
  EXPECT_EQ(K.rows(), mesh_.num_nodes());
  EXPECT_EQ(K.cols(), mesh_.num_nodes());
}

// ---- Pristine stiffness operator (no BCs applied) -------------------------

TEST_F(AssemblerTest, StiffnessIsSymmetric) {
  const Eigen::MatrixXd K = Dense();
  EXPECT_LT((K - K.transpose()).norm(), 1e-12);
}

// The bare Galerkin stiffness is only PSD: the constant vector is a null mode
// (grad of a constant is zero), so K * 1 = 0 and row sums vanish. It becomes
// SPD only after Dirichlet elimination -- checked separately below.
TEST_F(AssemblerTest, PristineStiffnessHasConstantNullSpace) {
  const Eigen::MatrixXd K = Dense();
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(mesh_.num_nodes());
  EXPECT_LT((K * ones).cwiseAbs().maxCoeff(), 1e-9);
}

// Every interior node carries the standard 5-point Laplacian diagonal (4),
// independent of node numbering.
TEST_F(AssemblerTest, InteriorNodesHaveStandardDiagonal) {
  const Eigen::MatrixXd K = Dense();
  std::vector<int> interior;
  for (int i = 0; i < mesh_.num_nodes(); ++i)
    if (!IsBoundary(i)) interior.push_back(i);
  EXPECT_EQ(interior.size(), 4u);
  for (int i : interior) EXPECT_NEAR(K(i, i), 4.0, 1e-9);
}

// Explicit regression on the known interior block for Rectangle(1, 4).
// Interior nodes are 5, 6, 9, 10; edge neighbors couple at -1 and the
// hypotenuse (diagonal) neighbors cancel to exactly 0 for right-isosceles P1.
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

  EXPECT_NEAR(K(5, 10), 0.0, 1e-9);  // hypotenuse couplings cancel
  EXPECT_NEAR(K(6, 9), 0.0, 1e-9);
}

// ---- Load path ------------------------------------------------------------

// For constant f == 1, sum_i ∫ f φ_i = ∫ f (sum_i φ_i) = ∫ 1 = area. On the
// unit square that is 1. Exact at any quadrature order (constant integrand),
// so this isolates shape() and the geometry factor from quadrature choice.
TEST_F(AssemblerTest, ConstantLoadSumsToDomainArea) {
  auto one = [](const aether::Vec2&, aether::Real) { return aether::Real(1); };
  assembler_.assemble_load(one);  // adds into the rhs zeroed by SetUp
  EXPECT_NEAR(assembler_.rhs()->sum(), 1.0, 1e-9);
}

// ---- Post-BC stiffness (Dirichlet applied) --------------------------------

// Dirichlet rows reduce to identity; boundary rhs entries hold the data g.
TEST_F(AssemblerTest, DirichletReducesBoundaryRowsToIdentity) {
  ApplyHarmonicDirichlet();
  const Eigen::MatrixXd K = Dense();
  const Eigen::VectorXd rhs = *assembler_.rhs();
  ASSERT_EQ(rhs.size(), mesh_.num_nodes());

  for (int i = 0; i < mesh_.num_nodes(); ++i) {
    if (!IsBoundary(i)) continue;
    EXPECT_NEAR(K(i, i), 1.0, 1e-9) << "boundary node " << i;
    const double off = K.row(i).cwiseAbs().sum() - std::abs(K(i, i));
    EXPECT_NEAR(off, 0.0, 1e-9) << "boundary node " << i;
    EXPECT_NEAR(rhs(i), ExactSolution(mesh_.node(i).x(), mesh_.node(i).y()), 1e-9)
        << "boundary node " << i;
  }
}

// SPD holds only AFTER elimination -- this is the property the TimeStepper's
// one-time Cholesky factorization of (M + theta*dt*K) depends on.
TEST_F(AssemblerTest, StiffnessIsPositiveDefiniteAfterBC) {
  ApplyHarmonicDirichlet();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<aether::Real>> llt(*assembler_.stiffness_matrix());
  EXPECT_EQ(llt.info(), Eigen::Success);
}

// End-to-end: the solve reproduces the harmonic exact solution at the nodes.
TEST_F(AssemblerTest, SolveMatchesExactSolution) {
  ApplyHarmonicDirichlet();
  aether::solver::Solver solver;
  const Eigen::VectorXd solution = solver.solve(*assembler_.stiffness_matrix(), *assembler_.rhs());
  ASSERT_EQ(solution.size(), mesh_.num_nodes());

  Eigen::VectorXd u_exact(mesh_.num_nodes());
  for (int i = 0; i < mesh_.num_nodes(); ++i)
    u_exact(i) = ExactSolution(mesh_.node(i).x(), mesh_.node(i).y());

  EXPECT_LT((solution - u_exact).cwiseAbs().maxCoeff(), 1e-9);
}

// ---- Consistent mass matrix -----------------------------------------------

class MassTest : public AssemblerTest {
 protected:
  void SetUp() override {
    assembler_.assemble_stiffness();
    assembler_.assemble_mass();
  }
};

TEST_F(MassTest, MassDimensions) {
  const auto& M = *assembler_.mass_matrix();
  EXPECT_EQ(M.rows(), mesh_.num_nodes());
  EXPECT_EQ(M.cols(), mesh_.num_nodes());
}

TEST_F(MassTest, MassIsSymmetric) {
  const Eigen::MatrixXd M = DenseMass();
  EXPECT_LT((M - M.transpose()).norm(), 1e-12);
}

// Consistent mass is strictly SPD (no null space, unlike the bare stiffness).
TEST_F(MassTest, MassIsPositiveDefinite) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<aether::Real>> llt(*assembler_.mass_matrix());
  EXPECT_EQ(llt.info(), Eigen::Success);
}

// Total mass = 1^T M 1 = ∫∫_Ω 1 dA = domain area. Each element contributes
// absdet/24 * (sum of the local pattern = 12) = absdet/2 = |T|, summing to the
// area. This is THE check on the absdet/24 normalization: /12 would read 2.0,
// /48 would read 0.5. Getting 1.0 pins the constant.
TEST_F(MassTest, MassTotalEqualsDomainArea) {
  const Eigen::MatrixXd M = DenseMass();
  EXPECT_NEAR(M.sum(), 1.0, 1e-9);
}

// The matrix must be the CONSISTENT mass, not a lumped (diagonal) one: the
// off-diagonal coupling from ∫ φ_i φ_j (i != j) is nonzero. This distinguishes
// a correct consistent-mass assembly from an accidental row-sum lumping.
TEST_F(MassTest, MassIsConsistentNotLumped) {
  const Eigen::MatrixXd M = DenseMass();
  const double off_diag_mass = M.cwiseAbs().sum() - M.diagonal().cwiseAbs().sum();
  EXPECT_GT(off_diag_mass, 1e-6);
}

TEST_F(MassTest, MassDiagonalIsPositive) {
  const Eigen::MatrixXd M = DenseMass();
  for (int i = 0; i < mesh_.num_nodes(); ++i) EXPECT_GT(M(i, i), 0.0);
}

}  // namespace
