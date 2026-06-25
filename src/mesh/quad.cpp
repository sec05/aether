// Copyright 2026 Spencer Evans-Cole

#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <utility>
#include <set>

#include "aether/mesh/quad.hpp"

namespace aether::mesh {

namespace {
// Bilinear map: reference [0,1]^2 -> physical quad.
// (xi,eta)=(0,0)->c0, (1,0)->c1, (1,1)->c2, (0,1)->c3.
Vec2 bilinear(const std::array<Vec2, 4>& c, Real xi, Real eta) {
  const Real w0 = (1 - xi) * (1 - eta);
  const Real w1 = xi * (1 - eta);
  const Real w2 = xi * eta;
  const Real w3 = (1 - xi) * eta;
  return Vec2(w0 * c[0].x() + w1 * c[1].x() + w2 * c[2].x() + w3 * c[3].x(),
              w0 * c[0].y() + w1 * c[1].y() + w2 * c[2].y() + w3 * c[3].y());
}

bool point_in_quad(const std::array<Vec2, 4>& q, const Vec2& p) {
  bool inside = false;
  for (int i = 0, j = 3; i < 4; j = i++) {
    if (((q[i].y() > p.y()) != (q[j].y() > p.y())) &&
        (p.x() < (q[j].x() - q[i].x()) * (p.y() - q[i].y()) / (q[j].y() - q[i].y()) + q[i].x())) {
      inside = !inside;
    }
  }
  return inside;
}
}  // namespace

Quad::Quad(int degree, int nx, int ny, std::array<Vec2, 4> corners,
           std::vector<std::array<Vec2, 4>> holes)
    : Mesh(degree, nx * ny), corners_(corners), nx_(nx), ny_(ny) {  // NOLINT

  const Real area = 0.5 * std::abs(corners[0].x() * (corners[1].y() - corners[3].y()) +
                                   corners[1].x() * (corners[2].y() - corners[0].y()) +
                                   corners[2].x() * (corners[3].y() - corners[1].y()) +
                                   corners[3].x() * (corners[0].y() - corners[2].y()));
  assert(area > 1e-12 && "Quad corners must be in winding order, not row-major");
  // 1. Full structured node grid (physical coords via bilinear map).
  std::vector<Vec2> all_nodes;
  all_nodes.reserve(nx * ny);
  for (int j = 0; j < ny; ++j) {
    const Real eta = static_cast<Real>(j) / (ny - 1);
    for (int i = 0; i < nx; ++i) {
      const Real xi = static_cast<Real>(i) / (nx - 1);
      all_nodes.push_back(bilinear(corners_, xi, eta));
    }
  }

  // 2. Triangles, dropping any whose centroid sits inside a hole.
  std::vector<std::array<Index, 3>> kept;
  kept.reserve(2 * (nx - 1) * (ny - 1));
  for (Index j = 0; j < ny - 1; ++j) {
    for (Index i = 0; i < nx - 1; ++i) {
      const Index n0 = j * nx + i, n1 = n0 + 1, n2 = n0 + nx, n3 = n2 + 1;
      for (const auto& t : {std::array<Index, 3>{n0, n1, n2}, std::array<Index, 3>{n1, n3, n2}}) {
        const Vec2 c((all_nodes[t[0]].x() + all_nodes[t[1]].x() + all_nodes[t[2]].x()) / 3,
                     (all_nodes[t[0]].y() + all_nodes[t[1]].y() + all_nodes[t[2]].y()) / 3);
        bool in_hole = false;
        for (const auto& h : holes)
          if (point_in_quad(h, c)) {
            in_hole = true;
            break;
          }
        if (!in_hole) kept.push_back(t);
      }
    }
  }

  // 3. Compact: drop unreferenced nodes but KEEP original row-major order, so
  //    the no-hole case is identical to the plain grid (remap == identity).
  std::vector<bool> used(all_nodes.size(), false);
  for (const auto& t : kept)
    for (Index v : t) used[v] = true;

  std::vector<Index> remap(all_nodes.size(), 0);
  nodes_.clear();
  nodes_.reserve(all_nodes.size());
  for (Index old = 0; old < static_cast<Index>(all_nodes.size()); ++old) {
    if (used[old]) {
      remap[old] = static_cast<Index>(nodes_.size());
      nodes_.push_back(all_nodes[old]);
    }
  }

  triangles_.clear();
  triangles_.reserve(kept.size());
  for (const auto& t : kept)
    triangles_.push_back({NodeIndex{remap[t[0]]}, NodeIndex{remap[t[1]]}, NodeIndex{remap[t[2]]}});
  num_nodes_ = static_cast<int>(nodes_.size());
  num_elements_ = static_cast<int>(triangles_.size());
  compute_boundary();
}

std::vector<NodeIndex> Quad::boundary_nodes() const { return boundary_nodes_; }

void Quad::compute_boundary() {
  struct EdgeRec {
    int count;
    NodeIndex a, b;
  };  // a,b = first-seen (triangle CCW) order
  std::map<std::pair<NodeIndex, NodeIndex>, EdgeRec> edges;
  auto sorted = [](NodeIndex x, NodeIndex y) {
    return x < y ? std::make_pair(x, y) : std::make_pair(y, x);
  };

  for (const auto& tri : triangles_)
    for (int e = 0; e < 3; ++e) {
      NodeIndex a = tri[e], b = tri[(e + 1) % 3];
      auto it = edges.find(sorted(a, b));
      if (it == edges.end())
        edges.insert({sorted(a, b), {1, a, b}});
      else
        it->second.count++;
    }

  boundary_facets_.clear();
  std::set<NodeIndex> bnodes;
  for (const auto& [k, rec] : edges) {
    if (rec.count != 1) continue;  // interior edge
    boundary_facets_.push_back({{rec.a, rec.b}, classify_edge(rec.a, rec.b)});
    bnodes.insert(rec.a);
    bnodes.insert(rec.b);
  }
  boundary_nodes_.assign(bnodes.begin(), bnodes.end());  // keep this for Dirichlet
}

int Quad::classify_edge(NodeIndex a, NodeIndex b) const {
  const Vec2 mid = 0.5 * (nodes_[a.get()] + nodes_[b.get()]);    // adapt to your node accessor
  const double tol = 1e-9 * (corners_[2] - corners_[0]).norm();  // ~bbox diagonal
  for (int s = 0; s < 4; ++s) {
    const Vec2 p = corners_[s], q = corners_[(s + 1) % 4];
    const Vec2 d = q - p, r = mid - p;
    const double cross = d.x() * r.y() - d.y() * r.x();
    const double t = r.dot(d) / d.dot(d);
    if (std::abs(cross) < tol * d.norm() && t > -tol && t < 1.0 + tol)
      return s + 1;  // outer sides -> 1,2,3,4
  }
  return 5;  // not on any outer side => hole rim
}
std::array<NodeIndex, 3> Quad::element_node_indices(int element_index) const {
  return triangles_[element_index];
}

std::vector<Vec2> Quad::element_nodes(int element_index) const {
  const std::array<NodeIndex, 3>& tri = triangles_[element_index];
  return {nodes_[tri[0].get()], nodes_[tri[1].get()], nodes_[tri[2].get()]};
}
}  // namespace aether::mesh
