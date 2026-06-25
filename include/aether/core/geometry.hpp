// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <functional>

#include <cmath>
#include <Eigen/Core>

#include "./types.hpp"  // Real

namespace aether {

// ---- Fixed-size geometric primitives --------------------------------------

using Vec2 = Eigen::Matrix<Real, 2, 1>;
using Vec3 = Eigen::Matrix<Real, 3, 1>;

using Mat2 = Eigen::Matrix<Real, 2, 2>;  // e.g. 2D map Jacobians
using Mat3 = Eigen::Matrix<Real, 3, 3>;  // e.g. 3D map Jacobians

// Aliases for intent at call sites — a Point is just a position vector.
using Point2 = Vec2;
using Point3 = Vec3;

// ---- 2D scalar cross ------------------------------------------------------

// z-component of the 3D cross product of (a, b) lying in the plane.
// > 0  =>  b is counter-clockwise from a. Equals twice the signed area of
// the triangle (origin, a, b). The workhorse for 2D area and orientation.
inline Real cross(const Vec2& a, const Vec2& b) noexcept { return a.x() * b.y() - a.y() * b.x(); }

// ---- Simplex measures -----------------------------------------------------

// Signed area of triangle (a, b, c) in 2D. Positive iff a, b, c are CCW.
// Sign is useful for mesh orientation checks; take std::abs for the measure.
inline Real signed_area(const Vec2& a, const Vec2& b, const Vec2& c) noexcept {
  return Real(0.5) * cross(b - a, c - a);
}

inline Real triangle_area(const Vec2& a, const Vec2& b, const Vec2& c) noexcept {
  return std::abs(signed_area(a, b, c));
}

// Triangle area in 3D: half the magnitude of (b-a) x (c-a). Cross is written
// out so this header only needs <Eigen/Core>, not <Eigen/Geometry>.
inline Real triangle_area(const Vec3& a, const Vec3& b, const Vec3& c) noexcept {
  const Vec3 u = b - a;
  const Vec3 v = c - a;
  const Vec3 n(u.y() * v.z() - u.z() * v.y(), u.z() * v.x() - u.x() * v.z(),
               u.x() * v.y() - u.y() * v.x());
  return Real(0.5) * n.norm();
}

// Signed volume of tetrahedron (a, b, c, d). Positive iff d is on the
// positive side of the CCW-oriented base triangle (a, b, c).
inline Real signed_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) noexcept {
  const Vec3 u = b - a;
  const Vec3 v = c - a;
  const Vec3 w = d - a;
  const Vec3 uxv(u.y() * v.z() - u.z() * v.y(), u.z() * v.x() - u.x() * v.z(),
                 u.x() * v.y() - u.y() * v.x());
  return uxv.dot(w) / Real(6);
}

inline Real tetrahedron_volume(const Vec3& a, const Vec3& b, const Vec3& c,
                               const Vec3& d) noexcept {
  return std::abs(signed_volume(a, b, c, d));
}

// ---- Centroids ------------------------------------------------------------
// Templated on the vector type so they serve Vec2 and Vec3 alike.

template <typename Vec>
inline Vec centroid(const Vec& a, const Vec& b, const Vec& c) noexcept {
  return (a + b + c) / Real(3);
}

template <typename Vec>
inline Vec centroid(const Vec& a, const Vec& b, const Vec& c, const Vec& d) noexcept {
  return (a + b + c + d) / Real(4);
}

// ---- Edge normal (2D) -----------------------------------------------------

// Normal to the directed edge a->b, obtained by rotating the edge direction
// -90 degrees (clockwise). For a CCW-oriented cell boundary this points
// outward. NOT normalized — the length equals the edge length, which is often
// exactly what you want for FVM flux integrals (n * |edge| in one shot).
inline Vec2 edge_normal(const Vec2& a, const Vec2& b) noexcept {
  const Vec2 d = b - a;
  return Vec2(d.y(), -d.x());  // flip sign for the opposite orientation
}

enum class BCType { Dirichlet, Neumann, Robin };

// g(x, t) — scalar for the heat problem; widen to a small vector later if needed
using BoundaryFn = std::function<double(const Vec2&, double)>;

struct BoundaryCondition {
  BCType type;
  int marker;              // boundary facet tag: outer wall, hole rim, etc.
  BoundaryFn value;        // Dirichlet g, Neumann flux q, or Robin g
  BoundaryFn robin_coeff;  // alpha in (alpha*u + du/dn = g); unused otherwise
};

}  // namespace aether
