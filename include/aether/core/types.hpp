#pragma once

#include <algorithm>
#include <cmath>
#include <compare>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>

#if __has_include(<expected>) && defined(__cpp_lib_expected)
  #include <expected>
  #define AETHER_HAS_STD_EXPECTED 1
#else
  #define AETHER_HAS_STD_EXPECTED 0
  // Not on C++23: alias Expected to a backport (e.g. tl::expected) instead.
#endif

namespace aether {

// ---- Scalar ---------------------------------------------------------------

using Real = double;

inline constexpr Real kEpsilon = std::numeric_limits<Real>::epsilon();

inline bool nearly_equal(Real a, Real b,
                         Real rel = 1e-9, Real abs_tol = 1e-12) noexcept {
  return std::abs(a - b) <= std::max(abs_tol, rel * std::max(std::abs(a), std::abs(b)));
}

// ---- Indices --------------------------------------------------------------

using Index = std::int32_t;
inline constexpr Index kInvalidIndex = -1;

// Distinct type per Tag, so node/cell/face/dof indices can't be silently
// swapped.
template <typename Tag>
class TypedIndex {
public:
  using value_type = Index;

  constexpr TypedIndex() noexcept = default;
  explicit constexpr TypedIndex(Index v) noexcept : value_(v) {}

  constexpr Index get() const noexcept { return value_; }
  explicit constexpr operator Index() const noexcept { return value_; }
  constexpr bool valid() const noexcept { return value_ >= 0; }

  friend constexpr auto operator<=>(TypedIndex, TypedIndex) = default;

  constexpr TypedIndex& operator++() noexcept { ++value_; return *this; }
  constexpr TypedIndex operator++(int) noexcept { auto t = *this; ++value_; return t; }

private:
  Index value_ = kInvalidIndex;
};

struct NodeTag; struct CellTag; struct FaceTag; struct EdgeTag; struct DofTag;

using NodeIndex = TypedIndex<NodeTag>;
using CellIndex = TypedIndex<CellTag>;
using FaceIndex = TypedIndex<FaceTag>;
using EdgeIndex = TypedIndex<EdgeTag>;
using DofIndex  = TypedIndex<DofTag>;

// ---- Error handling -------------------------------------------------------

enum class ErrorCode {
  kOk = 0,
  kInvalidArgument,
  kOutOfRange,
  kSingularMatrix,
  kNotConverged,
  kDimensionMismatch,
  kMeshInvalid,
  kNotImplemented,
};

struct Error {
  ErrorCode code = ErrorCode::kInvalidArgument;
  std::string message;
};

#if AETHER_HAS_STD_EXPECTED
template <typename T, typename E = Error>
using Expected = std::expected<T, E>;

inline std::unexpected<Error> fail(ErrorCode code, std::string msg = {}) {
  return std::unexpected<Error>{Error{code, std::move(msg)}};
}
#endif

} // namespace aether

// Hash so typed indices work as keys in unordered_map / unordered_set.
template <typename Tag>
struct std::hash<aether::TypedIndex<Tag>> {
  std::size_t operator()(aether::TypedIndex<Tag> i) const noexcept {
    return std::hash<aether::Index>{}(i.get());
  }
};