// Copyright 2026 Spencer Evans-Cole
#include <iostream>
#include <vector>

#include <spdlog/spdlog.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "aether/mesh/rectangle.hpp"
#include "aether/version.hpp"
#include "aether/core/types.hpp"
#include "aether/assembly/assembler.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/solver/solver.hpp"
#include "aether/output/output.hpp"
#include "aether/mesh/quad.hpp"
#include "aether/input/input.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';

  aether::input::Input input = aether::input::Input::parse("input.conf");
  int m = input.nx();
  int mm = input.ny();
  std::vector<std::array<aether::Vec2, 4>> holes = input.holes();
  std::array<aether::Vec2, 4> corners = input.corners();
  aether::mesh::Quad quad(1, m, mm, corners, holes);
  spdlog::info("Created quad mesh with {} nodes", quad.num_nodes());
  spdlog::info("Num mesh nodes: {}", quad.num_nodes());
  spdlog::info("Num mesh elements: {}", quad.num_elements());
  aether::elements::P1Element ref_element;
  aether::assembly::Assembler assembler(quad, ref_element);
  assembler.assemble();
  spdlog::info("Assembled {} x {} stiffness matrix with {} nonzeros",
               assembler.stiffness_matrix()->rows(), assembler.stiffness_matrix()->cols(),
               assembler.stiffness_matrix()->nonZeros());

  // Check if matrix is positive definite before solving
  Eigen::SimplicialLLT<Eigen::SparseMatrix<aether::Real>> cholesky(*assembler.stiffness_matrix());
  if (cholesky.info() != Eigen::Success) {
    spdlog::error("Stiffness matrix is not positive definite. Cannot solve.");
    return 1;
  }
  spdlog::info("Stiffness matrix is positive definite. Proceeding to solve.");

  aether::solver::Solver solver;
  Eigen::VectorXd solution = solver.solve(*assembler.stiffness_matrix(), *assembler.rhs());

  // Fine L_infty error
  Eigen::VectorXd exact_solution = Eigen::VectorXd::Zero(solution.size());
  const int n = quad.num_nodes();
  Eigen::VectorXd u_exact(n);
  for (int i = 0; i < n; ++i)
    u_exact[i] = std::pow(quad.node(i).x(), 2) - std::pow(quad.node(i).y(), 2);

  const double err_inf = (solution - u_exact).cwiseAbs().maxCoeff();
  spdlog::info("L_infty error: {:.2e}", err_inf);

  aether::output::Output output;
  output.write_solution_vtk("solution.vtu", quad.nodes(), quad.triangles(), {{"u", solution}});
  spdlog::info("Wrote solution to solution.vtu");
  return 0;
}
