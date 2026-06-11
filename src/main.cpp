#include <iostream>

#include <spdlog/spdlog.h>

#include "aether/version.hpp"

int main() {
  spdlog::set_pattern("[%^%l%$] %v");
  spdlog::info("Starting {}", aether::kProjectName);
  std::cout << aether::kProjectName << ' ' << aether::version_string() << '\n';
  return 0;
}
