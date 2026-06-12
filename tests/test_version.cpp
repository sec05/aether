#include <type_traits>
#include <variant>

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "aether/core/types.hpp"
#include "aether/version.hpp"

TEST(Version, ReturnsExpectedString) { EXPECT_EQ(aether::version_string(), "0.1.0"); }

TEST(Conventions, UsesConfiguredTypes) {
  static_assert(std::is_same_v<aether::Index, std::int32_t>);
  static_assert(std::is_same_v<aether::Numeric<>::ScalarType, double>);
  static_assert(std::is_same_v<aether::Numeric<float>::ScalarType, float>);
  static_assert(std::is_same_v<aether::Expected<int, int>, std::variant<int, int>>);
}

TEST(Eigen, MatrixMultiplicationWorks) {
  Eigen::Matrix2d left;
  left << 1.0, 2.0, 3.0, 4.0;

  Eigen::Matrix2d right;
  right << 5.0, 6.0, 7.0, 8.0;

  const Eigen::Matrix2d product = left * right;

  EXPECT_DOUBLE_EQ(product(0, 0), 19.0);
  EXPECT_DOUBLE_EQ(product(0, 1), 22.0);
  EXPECT_DOUBLE_EQ(product(1, 0), 43.0);
  EXPECT_DOUBLE_EQ(product(1, 1), 50.0);
}
