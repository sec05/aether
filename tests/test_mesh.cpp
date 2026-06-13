// Copyright 2026 Spencer Evans-Cole
#include <cstdint>
#include <type_traits>

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "aether/core/types.hpp"
#include "aether/mesh/square.hpp"
#include "aether/version.hpp"

TEST(SquareMesh, ConstructsCorrectly) {
  int degree = 1;
  int num_nodes = 16;  // 4x4 grid
  aether::mesh::Square square_mesh(degree, num_nodes);

  EXPECT_EQ(square_mesh.degree(), degree);
  EXPECT_EQ(square_mesh.num_nodes(), num_nodes);
  EXPECT_EQ(square_mesh.nodes().size(), num_nodes);

  // Check that the nodes are correctly spaced in a uniform grid
  int num_nodes_per_side = static_cast<int>(std::sqrt(num_nodes));
  for (int j = 0; j < num_nodes_per_side; ++j) {
    for (int i = 0; i < num_nodes_per_side; ++i) {
      const auto &node = square_mesh.nodes()[j * num_nodes_per_side + i];
      EXPECT_NEAR(node.x(), static_cast<aether::Real>(i) / (num_nodes_per_side - 1), 1e-6);
      EXPECT_NEAR(node.y(), static_cast<aether::Real>(j) / (num_nodes_per_side - 1), 1e-6);
    }
  }
}

TEST(SquareMesh, GeneratesCorrectTriangles) {
  const int degree = 1;
  const int num_nodes = 16;  // 4x4 grid
  aether::mesh::Square square_mesh(degree, num_nodes);

  const auto &triangles = square_mesh.triangles();

  const int nps = 4;  // num_nodes_per_side
  const std::size_t expected_count = 2 * (nps - 1) * (nps - 1);
  ASSERT_EQ(triangles.size(), expected_count);  // 18

  // Per-triangle sanity: every index in range, and no degenerate triangles.
  for (const auto &tri : triangles) {
    for (int v = 0; v < 3; ++v) {
      EXPECT_GE(tri[v].get(), 0);
      EXPECT_LT(tri[v].get(), num_nodes);
    }
    EXPECT_NE(tri[0].get(), tri[1].get());
    EXPECT_NE(tri[1].get(), tri[2].get());
    EXPECT_NE(tri[0].get(), tri[2].get());
  }

  // Exact connectivity: regenerate with the same scheme and compare in order.
  std::size_t k = 0;
  for (int j = 0; j < nps - 1; ++j) {
    for (int i = 0; i < nps - 1; ++i) {
      const int n0 = j * nps + i;
      const int n1 = n0 + 1;
      const int n2 = n0 + nps;
      const int n3 = n2 + 1;

      // lower triangle {n0, n1, n2}
      EXPECT_EQ(triangles[k][0].get(), n0);
      EXPECT_EQ(triangles[k][1].get(), n1);
      EXPECT_EQ(triangles[k][2].get(), n2);
      ++k;

      // upper triangle {n1, n3, n2}
      EXPECT_EQ(triangles[k][0].get(), n1);
      EXPECT_EQ(triangles[k][1].get(), n3);
      EXPECT_EQ(triangles[k][2].get(), n2);
      ++k;
    }
  }
}
