// Copyright 2026 Spencer Evans-Cole
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include <string>

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
#include "aether/solver/time_stepper.hpp"

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
  assembler.assemble_stiffness();  // K pristine
  assembler.assemble_mass();       // M pristine

  const aether::Real pi = std::acos(aether::Real(-1));
  auto u_exact_t = [pi](const aether::Vec2& p, aether::Real t) {
    return std::exp(-2 * pi * pi * t) * std::sin(pi * p.x()) * std::sin(pi * p.y());
  };

  std::set<int> markers;
  for (const auto& f : quad.boundary_facets()) markers.insert(f.marker);
  const aether::Real T_cold = 0.0;  // or a literal, e.g. 0.0
  auto cold = [T_cold](const aether::Vec2&, aether::Real) { return T_cold; };

  std::vector<aether::BoundaryCondition> bcs;
  bcs.push_back({aether::BCType::Dirichlet, 4, u_exact_t, cold});  // left
  // marker 3 (top) intentionally omitted => homogeneous Neumann => insulated
  aether::solver::TimeStepper stepper(quad, *assembler.stiffness_matrix(), *assembler.mass_matrix(),
                                      bcs,
                                      /*theta=*/input.theta(), /*dt=*/input.time_step());
  const aether::Real T_hot = 100.0;
  stepper.set_initial_condition([T_hot](const aether::Vec2&) { return T_hot; });
  std::vector<std::pair<double, std::string>> series;
  int steps = input.num_steps();
  aether::output::Output output;
  for (int step = 1; step <= steps; ++step) {
    stepper.step();
    spdlog::info("Step {}: t = {:.6f}", step, stepper.time());
    const std::string file = "solution_" + std::to_string(step) + ".vtu";
    output.write_solution_vtk(file, quad.nodes(), quad.triangles(),
                              {{"u", stepper.solution().cast<double>()}});
    series.emplace_back(static_cast<double>(stepper.time()), file);
  }
  output.write_pvd("solution.pvd", series);
  // open solution.pvd in ParaView, not the individual .vtu files

  return 0;
}
