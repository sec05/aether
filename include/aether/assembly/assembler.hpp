// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include "../elements/reference_elements.hpp"
#include "../mesh/mesh.hpp"

namespace aether::assembly {

class Assembler {
 public:
  Assembler(const mesh::Mesh& mesh, const elements::ReferenceElement& ref_element)
      : mesh_(mesh), ref_element_(ref_element) {}

  // Holds references + large matrices; an implicit copy would silently deep-copy
  // the matrices. Opinionated — drop these two lines if you actually need copies.
  Assembler(const Assembler&) = delete;
  Assembler& operator=(const Assembler&) = delete;

  void assemble_stiffness();
  void assemble_mass();
  void assemble_load(const BoundaryFn& f, Real t = Real(0));

  void apply_boundary_conditions(const std::vector<BoundaryCondition>& conditions,
                                 Real t = Real(0));
  void dirichlet_bc(const std::vector<std::pair<NodeIndex, Real>>& bcs);

  const Eigen::SparseMatrix<Real>& stiffness_matrix() const { return global_stiffness_matrix_; }
  const Eigen::SparseMatrix<Real>& mass_matrix() const { return global_mass_matrix_; }
  const Eigen::VectorX<Real>& rhs() const { return rhs_; }

 private:
  void scatter(int element_index, const Mat3& local,
               std::vector<Eigen::Triplet<Real>>* triplets) const;
  Mat3 local_stiffness_matrix(int element_index) const;
  Mat3 local_mass_matrix(int element_index) const;

  const mesh::Mesh& mesh_;
  const elements::ReferenceElement& ref_element_;
  Eigen::SparseMatrix<Real> global_stiffness_matrix_;
  Eigen::SparseMatrix<Real> global_mass_matrix_;
  Eigen::VectorX<Real> rhs_;
};

}  // namespace aether::assembly
