// Copyright 2026 Spencer Evans-Cole
#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "aether/analysis/error_norms.hpp"
#include "aether/assembly/assembler.hpp"
#include "aether/core/types.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/input/input.hpp"
#include "aether/mesh/quad.hpp"
#include "aether/output/output.hpp"
#include "aether/solver/solver.hpp"
#include "aether/solver/time_stepper.hpp"
#include "aether/version.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {} {}", aether::kProjectName, aether::version_string());

  // --- Mesh ---
  const aether::input::Input input = aether::input::Input::parse("input.conf");
  const aether::mesh::Quad quad(1, input.nx(), input.ny(), input.corners(), input.holes());
  spdlog::info("Created quad mesh: {} nodes, {} elements", quad.num_nodes(), quad.num_elements());

  // --- Assemble pristine K (both paths need it) ---
  aether::elements::P1Element ref_element;
  aether::assembly::Assembler assembler(quad, ref_element);
  assembler.assemble_stiffness();

  // --- Manufactured solution: u = e^{-2 α π² t} sin(πx) sin(πy).
  //     Transient: solves u_t = α Δu. Its t=0 spatial profile is what the
  //     steady branch reuses as its manufactured target. ---
  constexpr aether::Real pi = std::numbers::pi_v<aether::Real>;
  const aether::Real alpha = input.alpha();
  auto u_exact = [alpha](const aether::Vec2& p, aether::Real t) {
    return std::exp(-2 * alpha * pi * pi * t) * std::sin(pi * p.x()) * std::sin(pi * p.y());
  };

  // --- Impose the exact trace as Dirichlet on every boundary marker. ---
  std::set<int> markers;
  for (const auto& f : quad.boundary_facets()) markers.insert(f.marker);

  auto zero = [](const aether::Vec2&, aether::Real) { return aether::Real(0); };
  std::vector<aether::BoundaryCondition> bcs;
  for (int m : markers) bcs.push_back({aether::BCType::Dirichlet, m, u_exact, zero});

  aether::output::Output output;
  const int steps = input.num_steps();

  // ==================================================================
  //  Steady state (steps == 0): one direct solve, no time marching.
  // ==================================================================
  if (steps == 0) {
    // Sourceless steady limit of u_t = α Δu is u ≡ 0 — a useless target. So
    // manufacture u*(x,y) = sin(πx) sin(πy) as a steady solution of -Δu = f,
    // giving the consistent source f = -Δu* = 2π² sin(πx) sin(πy). α cancels
    // between operator and source, so it does not appear here.
    auto f_steady = [](const aether::Vec2& p, aether::Real) {
      return 2 * pi * pi * std::sin(pi * p.x()) * std::sin(pi * p.y());
    };
    assembler.assemble_load(f_steady);                          // rhs <- ∫ f φ
    assembler.apply_boundary_conditions(bcs, aether::Real(0));  // Dirichlet trace u*(·)

    aether::solver::Solver solver;
    const auto u = solver.solve(assembler.stiffness_matrix(), assembler.rhs());

    const auto err =
        aether::analysis::compute_errors(quad, ref_element, u, u_exact, aether::Real(0));
    spdlog::info("Steady solve: L2 = {:.3e}  Linf = {:.3e}", err.l2, err.linf);

    output.write_solution_vtk("solution_steady.vtu", quad.nodes(), quad.triangles(),
                              {{"u", u.cast<double>()}});
    return 0;
  }

  // ==================================================================
  //  Transient: march u_t = α Δu, checking the error each step.
  // ==================================================================
  assembler.assemble_mass();

  aether::solver::TimeStepper stepper(quad, assembler.stiffness_matrix(), assembler.mass_matrix(),
                                      bcs,
                                      /*theta=*/input.theta(),
                                      /*dt=*/input.time_step(),
                                      /*scale=*/alpha);

  // IC = exact solution at t = 0, i.e. sin(πx) sin(πy).
  stepper.set_initial_condition(
      [](const aether::Vec2& p) { return std::sin(pi * p.x()) * std::sin(pi * p.y()); });

  std::vector<std::pair<double, std::string>> series;
  series.reserve(static_cast<std::size_t>(steps));

  aether::analysis::ErrorNorms worst;  // max over the whole run
  for (int step = 1; step <= steps; ++step) {
    stepper.step();
    const aether::Real t = stepper.time();

    const auto err =
        aether::analysis::compute_errors(quad, ref_element, stepper.solution(), u_exact, t);
    worst.l2 = std::max(worst.l2, err.l2);
    worst.linf = std::max(worst.linf, err.linf);

    spdlog::info("Step {}: t = {:.6f}  L2 = {:.3e}  Linf = {:.3e}", step, t, err.l2, err.linf);

    const std::string file = "solution_" + std::to_string(step) + ".vtu";
    output.write_solution_vtk(file, quad.nodes(), quad.triangles(),
                              {{"u", stepper.solution().cast<double>()}});
    series.emplace_back(static_cast<double>(t), file);
  }
  output.write_pvd("solution.pvd", series);  // open this in ParaView, not the .vtu files

  spdlog::info("Worst over run: L2 = {:.3e}  Linf = {:.3e}", worst.l2, worst.linf);
  return 0;
}
