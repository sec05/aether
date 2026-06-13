// Copyright 2026 Spencer Evans-Cole
#pragma once

#include <string>

namespace aether {

inline constexpr const char *kProjectName = "aether";
inline constexpr int kVersionMajor = 0;
inline constexpr int kVersionMinor = 1;
inline constexpr int kVersionPatch = 0;

inline std::string version_string() {
  return std::to_string(kVersionMajor) + "." + std::to_string(kVersionMinor) + "." +
         std::to_string(kVersionPatch);
}

}  // namespace aether
