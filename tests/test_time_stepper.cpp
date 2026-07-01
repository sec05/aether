// Copyright 2026 Spencer Evans-Cole
//
// Tests for aether::solver::TimeStepper on the transient heat equation
//   M u' + K u = b(t),  theta-method time integration.
//
// Mesh: Quad(1, m, m) on the unit square, corners CCW from bottom-left, so the
// four sides carry markers 1..4 (see Quad::classify_edge). The stepper resolves
// Dirichlet BCs by marker, which is why Quad (not Rectangle) is used here.
//
// Analytic MMS: u(x,t) = exp(-2 pi^2 t) sin(pi x) sin(pi y) solves du/dt = Lap u
// with f == 0, homogeneous Dirichlet, IC sin(pi x) sin(pi y). The IC vanishes on
// the boundary, so it is BC-compatible and Crank-Nicolson does not ring.

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include <gtest/gtest.h>

#include <Eigen/SparseCore>

#include "aether/assembly/assembler.hpp"
#include "aether/core/geometry.hpp"
#include "aether/core/types.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/mesh/quad.hpp"
#include "aether/solver/time_stepper.hpp"

namespace {

using aether::Real;
using aether::Vec2;

const double kPi = std::acos(-1.0);

const std::array<Vec2, 4> kUnitSquare{{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};

// Owns mesh + assembler and assembles K and M so a stepper can be built against
// them. The mesh must outlive any TimeStepper (which stores a const Mesh&), so
// declare a Harness before the stepper in each test.
struct Harness {
  aether::mesh::Quad mesh;
  aether::elements::P1Element ref;
  aether::assembly::Assembler assembler;

  explicit Harness(int m) : mesh(1, m, m, kUnitSquare), assembler(mesh, ref) {
    assembler.assemble_stiffness();  // pristine K, zeroed rhs
    assembler.assemble_mass();       // consistent M
  }
};

// Dirichlet on all four outer sides (markers 1..4) with a single functor.
std::vector<aether::BoundaryCondition> AllDirichlet(aether::BoundaryFn g) {
  std::vector<aether::BoundaryCondition> bcs;
  for (int m = 1; m <= 4; ++m) bcs.push_back({aether::BCType::Dirichlet, m, g, {}});
  return bcs;
}

Real UExact(const Vec2& p, Real t) {
  return std::exp(-2 * kPi * kPi * t) * std::sin(kPi * p.x()) * std::sin(kPi * p.y());
}

Real LinfError(const Harness& h, const aether::solver::TimeStepper& s,
               const std::function<Real(const Vec2&, Real)>& exact) {
  const Real T = s.time();
  Real e = 0;
  for (int i = 0; i < h.mesh.num_nodes(); ++i)
    e = std::max(e, std::abs(s.solution()[i] - exact(h.mesh.node(i), T)));
  return e;
}

Real RmsError(const Harness& h, const aether::solver::TimeStepper& s,
              const std::function<Real(const Vec2&, Real)>& exact) {
  const Real T = s.time();
  const int n = h.mesh.num_nodes();
  Real sum = 0;
  for (int i = 0; i < n; ++i) {
    const Real d = s.solution()[i] - exact(h.mesh.node(i), T);
    sum += d * d;
  }
  return std::sqrt(sum / n);
}

// ---- Exactness / consistency (independent of mesh resolution) -------------

// A spatially and temporally constant field is an exact solution of the heat
// equation: du/dt = 0 and Lap(const) = 0. The theta-method must preserve it to
// machine precision -- and it only can if the passed-in K is PRISTINE (K*1 = 0).
// This test therefore also guards the "pass pristine K" contract.
TEST(TimeStepper, ConstantFieldIsPreservedExactly) {
  Harness h(6);
  auto five = [](const Vec2&, Real) { return Real(5); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(five),
                                /*theta=*/1.0, /*dt=*/0.1);
  s.set_initial_condition([](const Vec2&) { return Real(5); });
  s.advance(20);

  for (int i = 0; i < h.mesh.num_nodes(); ++i)
    EXPECT_NEAR(s.solution()[i], 5.0, 1e-10) << "node " << i;
}

TEST(TimeStepper, ZeroStaysZero) {
  Harness h(6);
  auto zero = [](const Vec2&, Real) { return Real(0); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/1.0, /*dt=*/0.05);
  s.set_initial_condition([](const Vec2&) { return Real(0); });
  s.advance(10);
  EXPECT_LT(s.solution().cwiseAbs().maxCoeff(), 1e-12);
}

// Dirichlet dofs are overwritten to g(x, t_new) each step, so with a spatially
// uniform, time-dependent g(x,t) = t every boundary node equals the current
// time exactly. Verifies g is evaluated at t_new, not t_0.
TEST(TimeStepper, TimeDependentDirichletHeldExactly) {
  Harness h(6);
  auto ramp = [](const Vec2&, Real t) { return t; };
  const Real dt = 0.01;
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(ramp),
                                /*theta=*/1.0, dt);
  s.set_initial_condition([](const Vec2&) { return Real(0); });
  s.advance(5);
  const Real expected = 5 * dt;

  for (const auto& node : h.mesh.boundary_nodes())
    EXPECT_NEAR(s.solution()[node.get()], expected, 1e-10) << "boundary node " << node.get();
}

TEST(TimeStepper, TimeAdvancesByDt) {
  Harness h(4);
  auto zero = [](const Vec2&, Real) { return Real(0); };
  const Real dt = 0.02;
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/1.0, dt);
  s.set_initial_condition([](const Vec2&) { return Real(0); });
  EXPECT_NEAR(s.time(), 0.0, 1e-12);
  s.step();
  EXPECT_NEAR(s.time(), dt, 1e-12);
  s.advance(4);
  EXPECT_NEAR(s.time(), 5 * dt, 1e-12);
}

// ---- Physical behavior (max principle, stability) -------------------------

// Homogeneous Dirichlet, nonnegative IC, no source: the discrete max principle
// (P1 on the non-obtuse right-isosceles triangulation, backward Euler) forbids
// new extrema, so the peak decays and no undershoot below the boundary appears.
TEST(TimeStepper, MaxDecaysNoUndershootBackwardEuler) {
  Harness h(9);  // (0.5,0.5) is a node; IC peak = 1 there
  auto zero = [](const Vec2&, Real) { return Real(0); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/1.0, /*dt=*/0.005);
  s.set_initial_condition(
      [](const Vec2& p) { return std::sin(kPi * p.x()) * std::sin(kPi * p.y()); });
  const Real max0 = s.solution().maxCoeff();
  s.advance(20);
  const Real max1 = s.solution().maxCoeff();

  EXPECT_LT(max1, max0);                      // peak decays
  EXPECT_GT(max1, 0.0);                       // but stays positive
  EXPECT_GT(s.solution().minCoeff(), -1e-9);  // no undershoot below 0
}

// Backward Euler is unconditionally stable: even an enormous dt must not blow
// up. The solution should stay bounded and still decay toward zero.
TEST(TimeStepper, BackwardEulerStableForHugeDt) {
  Harness h(9);
  auto zero = [](const Vec2&, Real) { return Real(0); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/1.0, /*dt=*/1.0);  // dt >> any CFL bound
  s.set_initial_condition(
      [](const Vec2& p) { return std::sin(kPi * p.x()) * std::sin(kPi * p.y()); });
  const Real max0 = s.solution().maxCoeff();
  s.advance(3);
  const Real maxT = s.solution().maxCoeff();

  EXPECT_TRUE(std::isfinite(maxT));
  EXPECT_LT(maxT, max0);
  EXPECT_GT(s.solution().minCoeff(), -1e-9);
}

// A positive source with homogeneous Dirichlet and zero IC must heat the
// interior. Confirms the load functor is actually wired into the RHS each step.
TEST(TimeStepper, PositiveSourceHeatsInterior) {
  Harness h(9);
  // Assemble a constant load (f == 1) vector; rhs_ was zeroed by assemble_stiffness.
  auto one = [](const Vec2&, Real) { return Real(1); };
  h.assembler.assemble_load(one);
  const Eigen::VectorX<Real> b = *h.assembler.rhs();

  auto zero = [](const Vec2&, Real) { return Real(0); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/1.0, /*dt=*/0.01);
  s.set_initial_condition([](const Vec2&) { return Real(0); });
  s.set_load([b](Real) { return b; });
  s.advance(10);

  EXPECT_GT(s.solution().maxCoeff(), 1e-6);         // interior heated up
  for (const auto& node : h.mesh.boundary_nodes())  // boundary still clamped
    EXPECT_NEAR(s.solution()[node.get()], 0.0, 1e-10);
}

// ---- Verification against the analytic decaying solution ------------------

// Crank-Nicolson, tiny dt so temporal error is negligible; the remaining error
// is the P1 spatial error. A wrong mass normalization would corrupt the decay
// rate and blow this well past the bound.
TEST(TimeStepper, DecayMmsErrorIsSmall) {
  Harness h(17);  // h = 1/16
  auto zero = [](const Vec2&, Real) { return Real(0); };
  aether::solver::TimeStepper s(h.mesh, *h.assembler.stiffness_matrix(), *h.assembler.mass_matrix(),
                                AllDirichlet(zero),
                                /*theta=*/0.5, /*dt=*/0.0005);
  s.set_initial_condition([](const Vec2& p) { return UExact(p, 0.0); });
  s.advance(40);  // T = 0.02

  EXPECT_LT(LinfError(h, s, &UExact), 2e-2);
}

// Spatial O(h^2): halving h should cut the discrete-L2 error by ~4. Temporal
// error is driven far below spatial via CN + tiny dt so the rate is clean.
TEST(TimeStepper, SpatialRefinementConvergesSecondOrder) {
  auto zero = [](const Vec2&, Real) { return Real(0); };
  auto ic = [](const Vec2& p) { return UExact(p, 0.0); };
  const Real dt = 0.0005;
  const int steps = 20;  // T = 0.01

  Harness h1(9);  // h = 1/8
  aether::solver::TimeStepper s1(h1.mesh, *h1.assembler.stiffness_matrix(),
                                 *h1.assembler.mass_matrix(), AllDirichlet(zero), 0.5, dt);
  s1.set_initial_condition(ic);
  s1.advance(steps);
  const Real e1 = RmsError(h1, s1, &UExact);

  Harness h2(17);  // h = 1/16
  aether::solver::TimeStepper s2(h2.mesh, *h2.assembler.stiffness_matrix(),
                                 *h2.assembler.mass_matrix(), AllDirichlet(zero), 0.5, dt);
  s2.set_initial_condition(ic);
  s2.advance(steps);
  const Real e2 = RmsError(h2, s2, &UExact);

  EXPECT_LT(e2, e1);  // refinement helps
  const Real ratio = e1 / e2;
  EXPECT_GT(ratio, 2.5);  // ~4 expected for O(h^2); wide band
  EXPECT_LT(ratio, 6.0);
}

}  // namespace
