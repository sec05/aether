// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <array>
#include <span>

#include "../core/geometry.hpp"
#include "../core/types.hpp"

namespace aether::elements {

struct QuadPoint {
  Vec2 point;
  Real weight;
};
// Reference elements for 2D. These are used to define the shape functions and quadrature
// rules for the finite element method.
// They are defined on a reference domain (e.g. the unit square or unit triangle)
// and can be mapped to physical elements in the mesh.
class ReferenceElement {
 public:
  explicit ReferenceElement(int num_dofs) : num_dofs_(num_dofs) {}
  virtual std::array<Real, 3> shape(const Vec2& point) const = 0;
  virtual std::array<Vec2, 3> shape_grad(const Vec2& point) const = 0;
  virtual std::span<const QuadPoint> quadrature_points(int order) const = 0;

 protected:
  ReferenceElement(const ReferenceElement&) = default;
  ReferenceElement& operator=(const ReferenceElement&) = default;
  ReferenceElement(ReferenceElement&&) = default;
  ReferenceElement& operator=(ReferenceElement&&) = default;

 private:
  int num_dofs_;
};
}  // namespace aether::elements
