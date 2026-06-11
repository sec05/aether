#include <type_traits>
#include <variant>

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
