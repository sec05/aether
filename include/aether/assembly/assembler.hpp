// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "../elements/reference_elements.hpp"
#include "../mesh/mesh.hpp"

namespace aether::assembly {
class Assembler {
 public:
  Assembler(const mesh::Mesh& mesh, const elements::ReferenceElement& ref_element)
      : mesh_(mesh), ref_element_(ref_element) {
    // Initialize the global stiffness matrix with the appropriate size
    int num_nodes = mesh_.num_nodes();
    global_stiffness_matrix_.resize(num_nodes, num_nodes);
  }
  ~Assembler() = default;
  void assemble();
  const Eigen::SparseMatrix<Real>* stiffness_matrix() const { return &global_stiffness_matrix_; }
  const Eigen::VectorXd* rhs() const { return &rhs_; }
  // assembler.hpp — add to the public section
  void apply_boundary_conditions(const std::vector<BoundaryCondition>& conditions,
                                 Real t = Real(0));
  // assembler.hpp — public
  void assemble_load(const BoundaryFn& f, Real t = Real(0));

 private:
  void dirichlet_bc(const std::vector<std::pair<NodeIndex, Real>>& bcs);
  void scatter(int i, const Mat3& local_matrix, std::vector<Eigen::Triplet<Real>>* triplets) const;
  Mat3 local_stiffness_matrix(int element_index) const;
  const mesh::Mesh& mesh_;
  const elements::ReferenceElement& ref_element_;
  Eigen::SparseMatrix<Real> global_stiffness_matrix_;
  Eigen::VectorXd rhs_;
};
}  // namespace aether::assembly
