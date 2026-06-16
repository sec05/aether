// Copyright 2026 Spencer Evans-Cole
#include <iostream>

#include <spdlog/spdlog.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "aether/mesh/square.hpp"
#include "aether/version.hpp"
#include "aether/assembly/assembler.hpp"
#include "aether/elements/p1_element.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';
  aether::mesh::Square square_mesh(1, 16);
  spdlog::info("Created square mesh with {} nodes", square_mesh.num_nodes());
  spdlog::info("Num mesh nodes: {}", square_mesh.num_nodes());
  spdlog::info("Num mesh elements: {}", square_mesh.num_elements());
  aether::elements::P1Element ref_element;
  aether::assembly::Assembler assembler(square_mesh, ref_element);
  assembler.assemble();
  spdlog::info("Assembled stiffness matrix with {} nonzeros",
               assembler.stiffness_matrix()->nonZeros());
  // print stiffness matrix for debugging
  const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
  std::ostringstream oss;
  oss << Eigen::MatrixXd(*assembler.stiffness_matrix()).format(fmt);
  spdlog::info("Stiffness matrix:\n{}", oss.str());
  const Eigen::IOFormat rhs_fmt(4, 0, ", ", "\n", "[", "]");
  std::ostringstream rhs_oss;
  rhs_oss << assembler.rhs()->format(rhs_fmt);
  spdlog::info("RHS vector:\n{}", rhs_oss.str());
  return 0;
}
