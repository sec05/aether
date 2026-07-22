// Copyright 2026 Spencer Evans-Cole

#include "aether/input/input.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace aether::input {
namespace {

std::string trim(const std::string& s) {
  const auto begin = s.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(begin, end - begin + 1);
}

std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

// Parses "(0,0), (4,0), (4,4), (0,4)" into four corners.
// Strips the structural punctuation and reads exactly eight numbers.
std::array<Vec2, 4> parse_corners(const std::string& value, int line_no) {
  std::string buf = value;
  for (char& c : buf)
    if (c == '(' || c == ')' || c == ',') c = ' ';

  std::istringstream iss(buf);
  std::array<Real, 8> v{};
  for (int i = 0; i < 8; ++i) {
    if (!(iss >> v[i]))
      throw std::runtime_error("Corners: expected 4 (x,y) pairs on line " +
                               std::to_string(line_no));
  }

  Real extra;
  if (iss >> extra)  // a 5th coordinate / stray token
    throw std::runtime_error("Corners: too many values on line " + std::to_string(line_no));

  return {{Vec2(v[0], v[1]), Vec2(v[2], v[3]), Vec2(v[4], v[5]), Vec2(v[6], v[7])}};
}

}  // namespace

std::vector<std::array<Vec2, 4>> parse_holes(const std::string& value, int line_no) {
  std::vector<std::array<Vec2, 4>> holes;
  std::istringstream iss(value);
  std::string hole_str;
  while (std::getline(iss, hole_str, ';')) {  // split on semicolons
    hole_str = trim(hole_str);
    if (!hole_str.empty()) {
      holes.push_back(parse_corners(hole_str, line_no));
    }
  }
  return holes;
}

Input Input::parse(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) throw std::runtime_error("Could not open input file: " + filename);

  Input input;
  bool has_name = false;
  bool has_corners = false;
  bool has_output = false;

  std::string line;
  int line_no = 0;
  while (std::getline(infile, line)) {
    ++line_no;
    const std::string trimmed = trim(line);
    if (trimmed.empty() || trimmed[0] == '#') continue;  // blanks and comments

    const auto colon = trimmed.find(':');
    if (colon == std::string::npos)
      throw std::runtime_error("Malformed line " + std::to_string(line_no) +
                               " (expected 'key: value'): " + trimmed);

    const std::string key = to_lower(trim(trimmed.substr(0, colon)));
    const std::string value = trim(trimmed.substr(colon + 1));

    if (key == "name") {
      input.name_ = value;
      has_name = true;
    } else if (key == "corners") {
      input.corners_ = parse_corners(value, line_no);
      has_corners = true;
    } else if (key == "holes") {
      input.holes_ = parse_holes(value, line_no);
    } else if (key == "output") {
      input.output_filename_ = value;
      has_output = true;
    } else if (key == "nx") {
      try {
        input.nx_ = std::stoi(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid Nx value on line " + std::to_string(line_no) + ": " +
                                 e.what());
      }
    } else if (key == "ny") {
      try {
        input.ny_ = std::stoi(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid Ny value on line " + std::to_string(line_no) + ": " +
                                 e.what());
      }
    } else if (key == "time step") {
      try {
        input.time_step_ = std::stod(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid time step value on line " + std::to_string(line_no) +
                                 ": " + e.what());
      }
    } else if (key == "num steps") {
      try {
        input.num_steps_ = std::stoi(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid number of steps value on line " +
                                 std::to_string(line_no) + ": " + e.what());
      }
    } else if (key == "theta") {
      try {
        input.theta_ = std::stod(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid theta value on line " + std::to_string(line_no) + ": " +
                                 e.what());
      }
    } else if (key == "alpha") {
      try {
        input.alpha_ = std::stod(value);
      } catch (const std::exception& e) {
        throw std::runtime_error("Invalid alpha value on line " + std::to_string(line_no) + ": " +
                                 e.what());
      }
    } else {
      throw std::runtime_error("Unknown key '" + key + "' on line " + std::to_string(line_no));
    }
  }

  if (!has_name) throw std::runtime_error("Input file missing required 'Name:' field");
  if (!has_corners) throw std::runtime_error("Input file missing required 'Corners:' field");
  if (!has_output) throw std::runtime_error("Input file missing required 'output:' field");

  return input;
}
}  // namespace aether::input
