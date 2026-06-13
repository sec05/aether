// Copyright 2026 Spencer Evans-Cole
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#if AETHER_HAS_STD_EXPECTED
#include <expected>
#endif

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "aether/core/types.hpp"
#include "aether/version.hpp"

TEST(Version, ReturnsExpectedString) { EXPECT_EQ(aether::version_string(), "0.1.0"); }

TEST(Conventions, UsesConfiguredTypes) {
  static_assert(std::is_same_v<aether::Index, std::int32_t>);
  static_assert(std::is_same_v<aether::Real, double>);

  // Strong index types share Index's layout but are mutually distinct,
  // so a NodeIndex can't be silently passed where a CellIndex is wanted.
  static_assert(std::is_same_v<aether::NodeIndex::value_type, aether::Index>);
  static_assert(!std::is_same_v<aether::NodeIndex, aether::CellIndex>);
  static_assert(sizeof(aether::NodeIndex) == sizeof(aether::Index));

#if AETHER_HAS_STD_EXPECTED
  static_assert(std::is_same_v<aether::Expected<int>, std::expected<int, aether::Error>>);
  static_assert(std::is_same_v<aether::Expected<int, int>, std::expected<int, int>>);
#endif
}

TEST(TypedIndex, DefaultsToInvalid) {
  constexpr aether::NodeIndex node;
  EXPECT_FALSE(node.valid());
  EXPECT_EQ(node.get(), aether::kInvalidIndex);
}

TEST(TypedIndex, WrapsAndComparesByValue) {
  constexpr aether::CellIndex small{2};
  constexpr aether::CellIndex large{5};

  EXPECT_TRUE(small.valid());
  EXPECT_EQ(small.get(), 2);
  EXPECT_LT(small, large);
  EXPECT_EQ(small, aether::CellIndex{2});
}

TEST(TypedIndex, IncrementsLikeAnIndex) {
  aether::DofIndex dof{0};

  EXPECT_EQ((dof++).get(), 0);  // post: returns old value
  EXPECT_EQ(dof.get(), 1);
  EXPECT_EQ((++dof).get(), 2);  // pre: returns new value
}

TEST(TypedIndex, HashesForUnorderedContainers) {
  std::unordered_map<aether::NodeIndex, int> seen;
  seen[aether::NodeIndex{7}] = 1;

  EXPECT_TRUE(seen.contains(aether::NodeIndex{7}));
  EXPECT_FALSE(seen.contains(aether::NodeIndex{8}));
}

TEST(NearlyEqual, ComparesWithinTolerance) {
  EXPECT_TRUE(aether::nearly_equal(1.0, 1.0 + 1e-13));
  EXPECT_FALSE(aether::nearly_equal(1.0, 1.1));
}

#if AETHER_HAS_STD_EXPECTED
TEST(Expected, CarriesValueOrError) {
  const aether::Expected<int> ok = 42;
  ASSERT_TRUE(ok.has_value());
  EXPECT_EQ(*ok, 42);

  const aether::Expected<int> bad = aether::fail(aether::ErrorCode::kSingularMatrix, "boom");
  ASSERT_FALSE(bad.has_value());
  EXPECT_EQ(bad.error().code, aether::ErrorCode::kSingularMatrix);
  EXPECT_EQ(bad.error().message, "boom");
}
#endif

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
