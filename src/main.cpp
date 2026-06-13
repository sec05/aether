// Copyright 2026 Spencer Evans-Cole
#include <iostream>

#include <spdlog/spdlog.h>

#include "aether/mesh/square.hpp"
#include "aether/version.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';
  aether::mesh::Square square_mesh(1, 16);
  spdlog::info("Created square mesh with {} nodes", square_mesh.num_nodes());
  spdlog::info("Mesh nodes: {}", square_mesh.nodes().size());
  for (const auto& node : square_mesh.nodes()) {
    spdlog::info("Node at ({:.2f}, {:.2f})", node.x(), node.y());
  }
  return 0;
}
