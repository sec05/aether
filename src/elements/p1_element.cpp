// Copyright 2026 Spencer Evans-Cole

#include "aether/elements/p1_element.hpp"

#include <vector>

namespace aether::elements {
std::array<Real, 3> P1Element::shape(const Vec2 &point) const {
  // Implement shape functions for P1 element (linear triangle)
  Real x = point.x();
  Real y = point.y();
  return {1 - x - y, x, y};  // N1, N2, N3
}

std::array<Vec2, 3> P1Element::shape_grad(const Vec2 &point) const {
  // Gradients of shape functions are constant for P1 elements
  return {Vec2(-1, -1), Vec2(1, 0), Vec2(0, 1)};  // grad(N1), grad(N2), grad(N3)
}

std::vector<std::pair<Vec2, Real>> P1Element::quadrature_points(int order) const {
  // Weights sum to 1/2 = area of the reference triangle.
  if (order <= 1) {
    // 1-point centroid rule: exact to degree 1 (enough for P1 stiffness).
    return {{Vec2(1.0 / 3.0, 1.0 / 3.0), 0.5}};
  }
  // 3-point rule: exact to degree 2 (mass matrix, linear load vectors).
  const Real w = 1.0 / 6.0;
  return {{Vec2(1.0 / 6.0, 1.0 / 6.0), w},
          {Vec2(2.0 / 3.0, 1.0 / 6.0), w},
          {Vec2(1.0 / 6.0, 2.0 / 3.0), w}};
}

}  // namespace aether::elements
