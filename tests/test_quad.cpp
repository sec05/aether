// Copyright 2026 Spencer Evans-Cole
#include <array>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <utility>

#include <gtest/gtest.h>

#include "aether/core/types.hpp"
#include "aether/mesh/quad.hpp"

namespace {

// Diamond (rotated square) centered at (2,2), corners ordered CCW:
//   c0 bottom (2,0), c1 right (4,2), c2 top (2,4), c3 left (0,2).
// Deliberately non-axis-aligned so the bilinear map is actually exercised.
const std::array<aether::Vec2, 4> kDiamond{{{2.0, 0.0}, {4.0, 2.0}, {2.0, 4.0}, {0.0, 2.0}}};

// Unit square, same corner convention, for the reduce-to-Rectangle check.
const std::array<aether::Vec2, 4> kUnitSquare{{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};

}  // namespace

TEST(QuadMesh, ConstructsCorrectly) {
  const int degree = 1;
  const int nx = 4;
  const int ny = 4;
  aether::mesh::Quad quad_mesh(degree, nx, ny, kDiamond);

  EXPECT_EQ(quad_mesh.degree(), degree);
  EXPECT_EQ(quad_mesh.num_nodes(), nx * ny);
  EXPECT_EQ(quad_mesh.nodes().size(), static_cast<std::size_t>(nx * ny));
}

// On the unit square the bilinear map must collapse to the same uniform grid
// the Rectangle class produces. This pins the (xi=i/(nx-1), eta=j/(ny-1))
// convention and the row-major j*nx+i layout.
TEST(QuadMesh, ReproducesUnitSquare) {
  const int degree = 1;
  const int nx = 4;
  const int ny = 4;
  aether::mesh::Quad quad_mesh(degree, nx, ny, kUnitSquare);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      const auto& node = quad_mesh.nodes()[j * nx + i];
      EXPECT_NEAR(node.x(), static_cast<aether::Real>(i) / (nx - 1), 1e-6);
      EXPECT_NEAR(node.y(), static_cast<aether::Real>(j) / (ny - 1), 1e-6);
    }
  }
}

// Independent (hand-computed) checks of the bilinear map on the diamond:
// the four reference-square corners must land exactly on c0..c3, edge midpoints
// on edge midpoints, and the center on the centroid. Uses a 3x3 grid so the
// interior node sits at (xi,eta)=(0.5,0.5).
TEST(QuadMesh, MapsCornersAndInterior) {
  const int degree = 1;
  const int m = 3;  // 3x3
  aether::mesh::Quad quad_mesh(degree, m, m, kDiamond);

  const auto& nodes = quad_mesh.nodes();

  // Corners: index 0 -> c0, m-1 -> c1, m*m-1 -> c2, (m-1)*m -> c3.
  EXPECT_NEAR(nodes[0].x(), 2.0, 1e-12);
  EXPECT_NEAR(nodes[0].y(), 0.0, 1e-12);

  EXPECT_NEAR(nodes[m - 1].x(), 4.0, 1e-12);
  EXPECT_NEAR(nodes[m - 1].y(), 2.0, 1e-12);

  EXPECT_NEAR(nodes[m * m - 1].x(), 2.0, 1e-12);
  EXPECT_NEAR(nodes[m * m - 1].y(), 4.0, 1e-12);

  EXPECT_NEAR(nodes[(m - 1) * m].x(), 0.0, 1e-12);
  EXPECT_NEAR(nodes[(m - 1) * m].y(), 2.0, 1e-12);

  // Bottom-edge midpoint (i=1,j=0): 0.5*(c0+c1) = (3,1).
  EXPECT_NEAR(nodes[1].x(), 3.0, 1e-12);
  EXPECT_NEAR(nodes[1].y(), 1.0, 1e-12);

  // Center (i=1,j=1): all weights 0.25 -> centroid (2,2).
  EXPECT_NEAR(nodes[m * 1 + 1].x(), 2.0, 1e-12);
  EXPECT_NEAR(nodes[m * 1 + 1].y(), 2.0, 1e-12);
}

// Connectivity is purely topological, so the expected indices match the
// Rectangle scheme regardless of corner geometry.
TEST(QuadMesh, GeneratesCorrectTriangles) {
  const int degree = 1;
  const int m = 4;  // 4x4
  aether::mesh::Quad quad_mesh(degree, m, m, kDiamond);

  const auto& triangles = quad_mesh.triangles();

  const std::size_t expected_count = 2 * (m - 1) * (m - 1);
  ASSERT_EQ(triangles.size(), expected_count);  // 18

  // Per-triangle sanity: every index in range, no degenerate triangles.
  for (const auto& tri : triangles) {
    for (int v = 0; v < 3; ++v) {
      EXPECT_GE(tri[v].get(), 0);
      EXPECT_LT(tri[v].get(), m * m);
    }
    EXPECT_NE(tri[0].get(), tri[1].get());
    EXPECT_NE(tri[1].get(), tri[2].get());
    EXPECT_NE(tri[0].get(), tri[2].get());
  }

  // Exact connectivity: regenerate with the same scheme and compare in order.
  std::size_t k = 0;
  for (int j = 0; j < m - 1; ++j) {
    for (int i = 0; i < m - 1; ++i) {
      const int n0 = j * m + i;
      const int n1 = n0 + 1;
      const int n2 = n0 + m;
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

// Boundary detection is topological: exactly the index ring of the structured
// grid, independent of geometry.
TEST(QuadMesh, BoundaryNodesFormRing) {
  const int degree = 1;
  const int m = 4;  // 4x4
  aether::mesh::Quad quad_mesh(degree, m, m, kDiamond);

  // Expected ring: i or j on a boundary index.
  std::set<int> expected;
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < m; ++i)
      if (i == 0 || i == m - 1 || j == 0 || j == m - 1) expected.insert(j * m + i);

  // Perimeter of an m x m grid: 2*m + 2*m - 4 (corners not double-counted).
  ASSERT_EQ(expected.size(), static_cast<std::size_t>(2 * m + 2 * m - 4));  // 12

  const auto boundary = quad_mesh.boundary_nodes();

  std::set<int> actual;
  for (const auto& node : boundary) {
    EXPECT_GE(node.get(), 0);
    EXPECT_LT(node.get(), m * m);
    actual.insert(node.get());
  }

  // No duplicates, and the set matches the expected ring exactly.
  EXPECT_EQ(actual.size(), boundary.size());
  EXPECT_EQ(actual, expected);

  // Interior nodes must be excluded.
  for (int j = 1; j < m - 1; ++j)
    for (int i = 1; i < m - 1; ++i) EXPECT_EQ(actual.count(j * m + i), 0u);
}

// Requires the CCW orientation guard in the Quad constructor (throws
// std::invalid_argument on non-positive signed area). Remove this test if you
// have not added that guard. The corners below are the degenerate case from
// earlier: c1 == c3, signed area 0.
/*
TEST(QuadMesh, RejectsDegenerateCorners) {
  const std::array<aether::Vec2, 4> degenerate{{{2.0, 0.0}, {0.0, 2.0}, {2.0, 4.0}, {0.0, 2.0}}};
  EXPECT_THROW(aether::mesh::Quad(1, 3, 3, degenerate), std::invalid_argument);
}
*/

namespace {
const std::array<aether::Vec2, 4> kSquare4{{{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}}};

// Corner cut [2,4]x[0,2]: square -> L-shape (reentrant corner at (2,2)).
const std::array<aether::Vec2, 4> kCornerHole{{{2.0, 0.0}, {4.0, 0.0}, {4.0, 2.0}, {2.0, 2.0}}};

// Centred hole [1,3]x[1,3]: square -> annulus (second boundary loop).
const std::array<aether::Vec2, 4> kCenterHole{{{1.0, 1.0}, {3.0, 1.0}, {3.0, 3.0}, {1.0, 3.0}}};

// Hole fully outside the domain: must be a no-op.
const std::array<aether::Vec2, 4> kFarHole{
    {{10.0, 10.0}, {12.0, 10.0}, {12.0, 12.0}, {10.0, 12.0}}};  // NOLINT
}  // namespace

// An explicit empty hole list must reproduce the plain grid (default-arg path).
TEST(QuadMeshHoles, EmptyHolesMatchesPlainGrid) {
  const int degree = 1, m = 5;
  aether::mesh::Quad plain(degree, m, m, kSquare4);
  aether::mesh::Quad empty(degree, m, m, kSquare4, {});

  EXPECT_EQ(plain.num_nodes(), empty.num_nodes());
  EXPECT_EQ(plain.triangles().size(), empty.triangles().size());
  EXPECT_EQ(plain.boundary_nodes().size(), empty.boundary_nodes().size());
}

// A hole disjoint from the domain removes nothing.
TEST(QuadMeshHoles, HoleOutsideDomainIsNoOp) {
  const int degree = 1, m = 5;
  aether::mesh::Quad quad(degree, m, m, kSquare4, {kFarHole});

  EXPECT_EQ(quad.num_nodes(), m * m);
  EXPECT_EQ(quad.triangles().size(), static_cast<std::size_t>(2 * (m - 1) * (m - 1)));
  EXPECT_EQ(quad.boundary_nodes().size(), static_cast<std::size_t>(2 * m + 2 * m - 4));
}

// Corner hole drops the bottom-right quadrant: 4 of 16 cells -> 8 of 32
// triangles gone; orphaned nodes (3,0),(4,0),(3,1),(4,1) -> 25-4=21 nodes.
TEST(QuadMeshHoles, CornerHoleRemovesQuadrant) {
  const int degree = 1, m = 5;
  aether::mesh::Quad quad(degree, m, m, kSquare4, {kCornerHole});

  EXPECT_EQ(quad.triangles().size(), static_cast<std::size_t>(24));
  EXPECT_EQ(quad.num_nodes(), 21);
  EXPECT_EQ(quad.nodes().size(), static_cast<std::size_t>(21));

  // No surviving triangle has its centroid inside the open hole (2,4)x(0,2).
  for (const auto& tri : quad.triangles()) {
    const auto& a = quad.nodes()[tri[0].get()];
    const auto& b = quad.nodes()[tri[1].get()];
    const auto& c = quad.nodes()[tri[2].get()];
    const aether::Real cx = (a.x() + b.x() + c.x()) / 3;
    const aether::Real cy = (a.y() + b.y() + c.y()) / 3;
    EXPECT_FALSE(cx > 2.0 && cx < 4.0 && cy > 0.0 && cy < 2.0);
  }
}

// The cut must turn up as new boundary: the reentrant corner and the two cut
// edges were interior on the full square. L-shape outline has 16 boundary nodes.
TEST(QuadMeshHoles, CornerHoleAddsReentrantBoundary) {
  const int degree = 1, m = 5;
  aether::mesh::Quad quad(degree, m, m, kSquare4, {kCornerHole});

  std::set<std::pair<int, int>> bd;
  for (const auto& n : quad.boundary_nodes()) {
    const auto& p = quad.nodes()[n.get()];
    bd.insert({static_cast<int>(std::lround(p.x())), static_cast<int>(std::lround(p.y()))});
  }

  EXPECT_EQ(quad.boundary_nodes().size(), static_cast<std::size_t>(16));
  EXPECT_EQ(bd.size(), static_cast<std::size_t>(16));  // no duplicates

  EXPECT_TRUE(bd.count({2, 2}));   // reentrant corner
  EXPECT_TRUE(bd.count({2, 1}));   // vertical cut x=2, interior on full square
  EXPECT_TRUE(bd.count({3, 2}));   // horizontal cut y=2, interior on full square
  EXPECT_FALSE(bd.count({1, 1}));  // genuinely interior, must stay off boundary
}

// Compaction must leave every triangle/boundary index in range of the smaller
// node array (guards against indices that still point into the pre-cut grid).
TEST(QuadMeshHoles, IndicesValidAfterCompaction) {
  const int degree = 1, m = 5;
  aether::mesh::Quad quad(degree, m, m, kSquare4, {kCornerHole});

  const int n = quad.num_nodes();
  for (const auto& tri : quad.triangles())
    for (int v = 0; v < 3; ++v) {
      EXPECT_GE(tri[v].get(), 0);
      EXPECT_LT(tri[v].get(), n);
    }
  for (const auto& bn : quad.boundary_nodes()) {
    EXPECT_GE(bn.get(), 0);
    EXPECT_LT(bn.get(), n);
  }
}

// Interior hole: 4 central cells -> 8 triangles gone; only the centre node
// (2,2) is orphaned, so 25-1=24 nodes. On this coarse grid the result is a
// one-cell-thick ring, so all 24 nodes lie on a boundary (16 outer + 8 inner).
TEST(QuadMeshHoles, InteriorHoleCreatesAnnulus) {
  const int degree = 1, m = 5;
  aether::mesh::Quad quad(degree, m, m, kSquare4, {kCenterHole});

  EXPECT_EQ(quad.triangles().size(), static_cast<std::size_t>(24));
  EXPECT_EQ(quad.num_nodes(), 24);

  bool has_center = false;
  for (const auto& p : quad.nodes())
    if (std::lround(p.x()) == 2 && std::lround(p.y()) == 2) has_center = true;
  EXPECT_FALSE(has_center);  // centre node removed

  std::set<std::pair<int, int>> bd;
  for (const auto& n : quad.boundary_nodes()) {
    const auto& p = quad.nodes()[n.get()];
    bd.insert({static_cast<int>(std::lround(p.x())), static_cast<int>(std::lround(p.y()))});
  }
  EXPECT_EQ(quad.boundary_nodes().size(), static_cast<std::size_t>(24));

  const std::array<std::pair<int, int>, 8> inner_ring{
      {{1, 1}, {2, 1}, {3, 1}, {3, 2}, {3, 3}, {2, 3}, {1, 3}, {1, 2}}};
  for (const auto& rc : inner_ring) EXPECT_TRUE(bd.count(rc));
}

// Two disjoint holes must both apply (exercises the loop over the hole list).
TEST(QuadMeshHoles, MultipleDisjointHoles) {
  const int degree = 1, m = 5;
  const std::array<aether::Vec2, 4> hole_bl{{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
  const std::array<aether::Vec2, 4> hole_tr{{{3.0, 3.0}, {4.0, 3.0}, {4.0, 4.0}, {3.0, 4.0}}};
  aether::mesh::Quad quad(degree, m, m, kSquare4, {hole_bl, hole_tr});

  // One corner cell at each opposite corner: 4 triangles, 2 orphaned nodes.
  EXPECT_EQ(quad.triangles().size(), static_cast<std::size_t>(28));
  EXPECT_EQ(quad.num_nodes(), 23);

  bool bl = false, tr = false;
  for (const auto& p : quad.nodes()) {
    if (std::lround(p.x()) == 0 && std::lround(p.y()) == 0) bl = true;
    if (std::lround(p.x()) == 4 && std::lround(p.y()) == 4) tr = true;
  }
  EXPECT_FALSE(bl);
  EXPECT_FALSE(tr);
}
