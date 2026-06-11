#pragma once

#include <cstdint>
#include <variant>

namespace aether {

using Index = std::int32_t;

template <typename T, typename E> using Expected = std::variant<T, E>;

template <typename Scalar = double> struct Numeric {
  using ScalarType = Scalar;
};

} // namespace aether
