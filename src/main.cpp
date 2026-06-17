// Copyright 2026 Spencer Evans-Cole
#include <iostream>

#include <spdlog/spdlog.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "aether/mesh/square.hpp"
#include "aether/version.hpp"
#include "aether/core/types.hpp"
#include "aether/assembly/assembler.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/solver/solver.hpp"
int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';
  aether::mesh::Square square_mesh(1, 16);
  spdlog::info("Created square mesh with {} nodes", square_mesh.num_nodes());
  spdlog::info("Num mesh nodes: {}", square_mesh.num_nodes());
  spdlog::info("Num mesh elements: {}", square_mesh.num_elements());
  aether::elements::P1Element ref_element;
  aether::assembly::Assembler assembler(square_mesh, ref_element);
  assembler.assemble();
  spdlog::info("Assembled stiffness matrix with {} nonzeros",
               assembler.stiffness_matrix()->nonZeros());
  // print stiffness matrix for debugging
  const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
  std::ostringstream oss;
  oss << Eigen::MatrixXd(*assembler.stiffness_matrix()).format(fmt);
  spdlog::info("Stiffness matrix:\n{}", oss.str());
  const Eigen::IOFormat rhs_fmt(4, 0, ", ", "\n", "[", "]");
  std::ostringstream rhs_oss;
  rhs_oss << assembler.rhs()->format(rhs_fmt);
  spdlog::info("RHS vector:\n{}", rhs_oss.str());

  // Check if matrix is positive definite before solving
  Eigen::SimplicialLLT<Eigen::SparseMatrix<aether::Real>> cholesky(*assembler.stiffness_matrix());
  if (cholesky.info() != Eigen::Success) {
    spdlog::error("Stiffness matrix is not positive definite. Cannot solve.");
    return 1;
  }
  spdlog::info("Stiffness matrix is positive definite. Proceeding to solve.");

  aether::solver::Solver solver;
  Eigen::VectorXd solution = solver.solve(*assembler.stiffness_matrix(), *assembler.rhs());
  std::ostringstream sol_oss;
  sol_oss << solution.format(rhs_fmt);
  spdlog::info("Solution vector:\n{}", sol_oss.str());

  // Fine L_infty error
  Eigen::VectorXd exact_solution = Eigen::VectorXd::Zero(solution.size());
  const int n = square_mesh.num_nodes();
  Eigen::VectorXd u_exact(n);
  for (int i = 0; i < n; ++i)
    u_exact[i] = std::pow(square_mesh.node(i).x(), 2) - std::pow(square_mesh.node(i).y(), 2);

  const double err_inf = (solution - u_exact).cwiseAbs().maxCoeff();
  spdlog::info("L_infty error: {:.2e}", err_inf);
  return 0;
}
