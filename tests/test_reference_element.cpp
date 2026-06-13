// Copyright 2026 Spencer Evans-Cole
#include <cstdint>
#include <type_traits>

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "aether/core/types.hpp"
#include "aether/elements/p1_element.hpp"
#include "aether/version.hpp"
TEST(P1Element, ShapeFunctions) {
  aether::elements::P1Element p1;
  auto shape_at_center = p1.shape({0.3333333, 0.3333333});
  EXPECT_NEAR(shape_at_center[0], 0.3333333, 1e-6);
  EXPECT_NEAR(shape_at_center[1], 0.3333333, 1e-6);
  EXPECT_NEAR(shape_at_center[2], 0.3333333, 1e-6);
}

TEST(P1Element, ShapeGradients) {
  aether::elements::P1Element p1;
  auto grads = p1.shape_grad({0.0, 0.0});
  EXPECT_NEAR(grads[0].x(), -1.0, 1e-6);
  EXPECT_NEAR(grads[0].y(), -1.0, 1e-6);
  EXPECT_NEAR(grads[1].x(), 1.0, 1e-6);
  EXPECT_NEAR(grads[1].y(), 0.0, 1e-6);
  EXPECT_NEAR(grads[2].x(), 0.0, 1e-6);
  EXPECT_NEAR(grads[2].y(), 1.0, 1e-6);
}

TEST(P1Element, QuadraturePoints) {
  aether::elements::P1Element p1;
  auto qps = p1.quadrature_points(1);
  ASSERT_EQ(qps.size(), 1);
  EXPECT_NEAR(qps[0].first.x(), 0.3333333, 1e-6);
  EXPECT_NEAR(qps[0].first.y(), 0.3333333, 1e-6);
  EXPECT_NEAR(qps[0].second, 0.5, 1e-6);
}
