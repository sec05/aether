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
  // Since the shape functions are linear, a single quadrature point at the centroid
  // is sufficient for exact integration of polynomials up to degree 2.
  return {{Vec2(1.0 / 3.0, 1.0 / 3.0), 0.5}};  // (point, weight)
}

}  // namespace aether::elements
