// Copyright 2026 Spencer Evans-Cole
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "aether/assembly/assembler.hpp"
#include "aether/core/types.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/input/input.hpp"
#include "aether/mesh/quad.hpp"
#include "aether/output/output.hpp"
#include "aether/solver/solver.hpp"
#include "aether/version.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';

  // --- Mesh ---
  const aether::input::Input input = aether::input::Input::parse("input.conf");
  const aether::mesh::Quad quad(1, input.nx(), input.ny(), input.corners(), input.holes());
  spdlog::info("Created quad mesh: {} nodes, {} elements", quad.num_nodes(), quad.num_elements());

  // --- Assemble stiffness ---
  aether::elements::P1Element ref_element;
  aether::assembly::Assembler assembler(quad, ref_element);
  assembler.assemble();
  spdlog::info("Assembled {} x {} stiffness matrix with {} nonzeros",
               assembler.stiffness_matrix()->rows(), assembler.stiffness_matrix()->cols(),
               assembler.stiffness_matrix()->nonZeros());

  // --- Manufactured solution: u = sin(pi x) sin(pi y), so -Laplacian(u) = 2 pi^2 u. ---
  // u vanishes on the unit-square boundary, giving homogeneous Dirichlet data.
  const aether::Real pi = std::acos(aether::Real(-1));
  auto u_exact = [pi](const aether::Vec2& p, aether::Real) {
    return std::sin(pi * p.x()) * std::sin(pi * p.y());
  };
  auto source = [pi](const aether::Vec2& p, aether::Real) {
    return 2 * pi * pi * std::sin(pi * p.x()) * std::sin(pi * p.y());
  };

  assembler.assemble_load(source);

  std::set<int> markers;
  for (const auto& f : quad.boundary_facets()) markers.insert(f.marker);
  std::vector<aether::BoundaryCondition> bcs;
  for (int marker : markers) bcs.push_back({aether::BCType::Dirichlet, marker, u_exact, {}});
  assembler.apply_boundary_conditions(bcs);

  // --- Solve ---
  aether::solver::Solver solver;
  const Eigen::VectorXd solution = solver.solve(*assembler.stiffness_matrix(), *assembler.rhs());

  // --- L_infty error against the same manufactured field. ---
  const int n = quad.num_nodes();
  Eigen::VectorXd exact(n);
  for (int i = 0; i < n; ++i) exact[i] = u_exact(quad.node(i), aether::Real(0));

  const double err_inf = (solution - exact).cwiseAbs().maxCoeff();
  spdlog::info("L_infty error: {:.2e}", err_inf);

  // --- Output ---
  aether::output::Output output;
  output.write_solution_vtk("solution.vtu", quad.nodes(), quad.triangles(), {{"u", solution}});
  spdlog::info("Wrote solution to solution.vtu");
  return 0;
}
