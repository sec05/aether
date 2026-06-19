// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <array>
#include <string>
#include <vector>

#include "../core/geometry.hpp"

namespace aether::input {
class Input {
 public:
  // Parses the input file and returns an Input populated with the parsed data.
  // Throws std::runtime_error if the file cannot be opened or is malformed.
  static Input parse(const std::string& filename);

  const std::string& name() const { return name_; }
  const std::array<Vec2, 4>& corners() const { return corners_; }
  const std::string& output_filename() const { return output_filename_; }
  int nx() const { return nx_; }
  int ny() const { return ny_; }
  const std::vector<std::array<Vec2, 4>>& holes() const { return holes_; }

 private:
  std::string name_;
  std::array<Vec2, 4> corners_;
  std::string output_filename_;
  int nx_ = 0;
  int ny_ = 0;
  std::vector<std::array<Vec2, 4>> holes_;
};

}  // namespace aether::input
