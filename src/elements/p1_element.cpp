// Copyright 2026 Spencer Evans-Cole

#include "aether/elements/p1_element.hpp"

namespace aether::elements {

std::array<Real, 3> P1Element::shape(const Vec2& point) const {
  const Real x = point.x();
  const Real y = point.y();
  return {Real(1) - x - y, x, y};  // N1, N2, N3
}

std::array<Vec2, 3> P1Element::shape_grad(const Vec2& /*point*/) const {
  // Constant for P1: grad(N1), grad(N2), grad(N3).
  return {Vec2(-1, -1), Vec2(1, 0), Vec2(0, 1)};
}

std::span<const QuadPoint> P1Element::quadrature_points(int order) const {
  // Weights sum to 1/2 = area of the reference triangle. Static storage: the
  // rules are fixed, so return a view instead of allocating per call.
  static const std::array<QuadPoint, 1> kDegree1{{
      {Vec2(Real(1) / Real(3), Real(1) / Real(3)), Real(0.5)},
  }};
  static const std::array<QuadPoint, 3> kDegree2{{
      {Vec2(Real(1) / Real(6), Real(1) / Real(6)), Real(1) / Real(6)},
      {Vec2(Real(2) / Real(3), Real(1) / Real(6)), Real(1) / Real(6)},
      {Vec2(Real(1) / Real(6), Real(2) / Real(3)), Real(1) / Real(6)},
  }};
  return (order <= 1) ? std::span<const QuadPoint>{kDegree1} : std::span<const QuadPoint>{kDegree2};
}

}  // namespace aether::elements
