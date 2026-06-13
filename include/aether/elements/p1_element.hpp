// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <array>
#include <vector>

#include "./reference_elements.hpp"

namespace aether::elements {
class P1Element : public ReferenceElement {
 public:
  P1Element() : ReferenceElement(3) {}
  ~P1Element() override = default;
  std::array<Real, 3> shape(const Vec2 &point) const override;
  std::array<Vec2, 3> shape_grad(const Vec2 &point) const override;
  std::vector<std::pair<Vec2, Real>> quadrature_points(int order) const override;
};
}  // namespace aether::elements
