// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <array>
#include <span>

#include "./reference_elements.hpp"

namespace aether::elements {
class P1Element final : public ReferenceElement {
 public:
  P1Element() : ReferenceElement(3) {}
  std::array<Real, 3> shape(const Vec2 &point) const override;
  std::array<Vec2, 3> shape_grad(const Vec2 &point) const override;
  std::span<const QuadPoint> quadrature_points(int order) const override;
};
}  // namespace aether::elements
