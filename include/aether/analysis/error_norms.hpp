// Copyright 2026 Spencer Evans-Cole
#pragma once
#include <Eigen/Dense>
#include "../core/types.hpp"
#include "../elements/reference_elements.hpp"
#include "../mesh/mesh.hpp"

namespace aether::analysis {

struct ErrorNorms {
  Real l2 = Real(0);
  Real linf = Real(0);
};

// ‖u_h − u‖ over the whole mesh; u_exact sampled in physical space at time t.
ErrorNorms compute_errors(const mesh::Mesh& mesh, const elements::ReferenceElement& ref,
                          const Eigen::VectorX<Real>& u_h, const BoundaryFn& u_exact,
                          Real t = Real(0));

}  // namespace aether::analysis
