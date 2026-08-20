// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PROFILER_H
#define PROFILER_H

/// @file Profiler.h
/// @brief Optional profiling: phase timings, operation counts by operand shape,
/// and container allocation totals.
///
/// Compiled out entirely unless configured with -DENABLE_PROFILING=ON.
///
/// Usage:
/// @code
///   SVMP_PROFILE_PHASE("newton_assembly");                  // times this scope
///   SVMP_PROFILE_OP(mat_mul_return, rows, inner, cols);     // counts one call
/// @endcode
///
/// Two different mechanisms on purpose. Phases wrap stage functions and fire
/// O(1e4) times a run, so a string-keyed map costs nothing next to the work
/// they measure. Operations fire O(1e8) times from inside element and Gauss
/// point loops, where a string hash would cost several times the 3x3 product
/// being counted -- and would then distort the phase timings gathered from the
/// same run. Operation counting is therefore a fixed table plus a one-slot
/// memo: a compare and an increment on the hot path.
///
/// Never place SVMP_PROFILE_PHASE inside an element or Gauss point loop; two
/// clock reads would swamp the body.
///
/// @note The guards are spelled #ifdef. Written as a plain #if on an undefined
/// macro the preprocessor evaluates them as 0 and silently drops everything,
/// which is exactly how the older Array_gather_stats counters came to be dead
/// code that nobody noticed for years.

// ---------------------------------------------------------------------------
// Declared unconditionally: these emit no code, and naming them in the
// disabled-build macros below is what keeps call sites type-checked even when
// profiling is off. Without that, a renamed variable or a dropped operation
// compiles clean for everyone and only breaks for whoever turns the flag on.
// ---------------------------------------------------------------------------
namespace profiler {

/// @brief The mat_fun helpers that allocate, and are therefore worth counting.
enum class Op : unsigned {
  mat_det, mat_dev, mat_dyad_prod, mat_id, mat_inv,
  mat_mul_vec, mat_mul_return, mat_mul_inplace,
  mat_symm, mat_symm_prod,
  ten_asym_prod12, ten_dyad_prod, ten_mddot, ten_symm_prod,
  transpose,
  COUNT
};

}  // namespace profiler

#ifdef ENABLE_PROFILING
#define svmp_profiling_enabled
#endif

#ifdef svmp_profiling_enabled

#include <cstdint>
#include <string>
#include <unordered_map>

#include "Timer.h"

class cmType;
class CmMod;
class SimulationLogger;

namespace profiler {

/// @brief Largest operand dimension representable in a shape key.
/// The value above it is reserved to mark saturation, so a clamped dimension is
/// reported as such rather than silently aliasing onto a real shape.
constexpr unsigned DIM_MAX = 1022;
constexpr unsigned DIM_SATURATED = 1023;
constexpr unsigned DIM_NA = 0;

/// @brief Pack an operand shape into three 10-bit fields.
///
/// Convention: M is rows of the result, K the contracted inner dimension (-1
/// where there is none, e.g. an inverse), N the columns of the result.
inline unsigned shape_key(int m, int k, int n)
{
  auto field = [](int v) -> unsigned {
    if (v < 0) return DIM_NA;
    return (static_cast<unsigned>(v) > DIM_MAX) ? DIM_SATURATED : static_cast<unsigned>(v);
  };
  return (field(m) << 20) | (field(k) << 10) | field(n);
}

void unpack_shape(unsigned key, int& m, int& k, int& n);

/// @brief Name of an operation, for reporting.
const char* op_name(Op op);

/// @brief Counts for one operation, bucketed by operand shape.
///
/// Plain arrays, deliberately: the profiler must not allocate through Array,
/// Vector or std::unordered_map on the hot path, or it would perturb the very
/// allocation totals it exists to report.
struct OpTable {
  static constexpr int CAP = 32;

  unsigned      keys[CAP];
  std::uint64_t calls[CAP];
  int           last;             ///< memo slot; one shape dominates any loop
  int           n_used;
  std::uint64_t overflow_calls;   ///< counted but unbucketed, if CAP is hit

  void slow_add(unsigned key);
};

struct PhaseRecord {
  double        seconds{0.0};
  std::uint64_t calls{0};
  int           order{0};         ///< first-seen, so the report keeps source order
};

struct Registry {
  OpTable ops[static_cast<unsigned>(Op::COUNT)];
  std::unordered_map<std::string, PhaseRecord> phases;
  int phase_order{0};

  void add_phase(const char* name, double seconds);

  /// @brief Reduce across ranks and emit the report. Collective: every rank
  /// must call it. Only the master rank writes.
  void report(const CmMod& cm_mod, const cmType& cm, SimulationLogger& logger,
              const std::string& output_dir);
};

/// @brief Namespace-scope, not a function-local static: a static would add a
/// guard-variable load and branch to every counted call.
extern Registry registry;

/// @brief Hot path. Two loads, a compare, an increment on the common case.
inline void count_op(Op op, unsigned key)
{
  OpTable& table = registry.ops[static_cast<unsigned>(op)];
  if (table.keys[table.last] == key) {
    ++table.calls[table.last];
    return;
  }
  table.slow_add(key);
}

/// @brief Times its enclosing scope, including on early return or unwind.
class ScopedPhase {
  public:
    explicit ScopedPhase(const char* name) : name_(name) { timer_.set_time(); }
    ~ScopedPhase() { registry.add_phase(name_, timer_.get_elapsed_time()); }

    ScopedPhase(const ScopedPhase&) = delete;
    ScopedPhase& operator=(const ScopedPhase&) = delete;

  private:
    const char* name_;
    Timer timer_;
};

}  // namespace profiler

#if defined(_OPENMP)
#pragma message("svmp profiler counters are not thread safe; add atomics before parallelising instrumented regions.")
#endif

#define SVMP_PROFILE_CONCAT_(a, b) a##b
#define SVMP_PROFILE_CONCAT(a, b) SVMP_PROFILE_CONCAT_(a, b)

#define SVMP_PROFILE_PHASE(name) \
  ::profiler::ScopedPhase SVMP_PROFILE_CONCAT(svmp_phase_, __LINE__)(name)

#define SVMP_PROFILE_OP(op, m, k, n) \
  ::profiler::count_op(::profiler::Op::op, ::profiler::shape_key((m), (k), (n)))

/// @brief Record a duration that something else already measured, e.g. a
/// Teuchos::Time that the linear-algebra backend keeps for its own reporting.
/// Avoids adding a second clock around the same region.
#define SVMP_PROFILE_ADD(name, seconds) \
  ::profiler::registry.add_phase(name, (seconds))

#else  // svmp_profiling_enabled

// sizeof operands are unevaluated, so these emit no instructions at any
// optimisation level while still requiring every name and type to resolve.
#define SVMP_PROFILE_PHASE(name) ((void) sizeof(name))
#define SVMP_PROFILE_ADD(name, seconds) ((void) sizeof(name), (void) sizeof(seconds))
#define SVMP_PROFILE_OP(op, m, k, n) \
  ((void) sizeof(::profiler::Op::op), (void) sizeof((m) + (k) + (n)))

#endif  // svmp_profiling_enabled

#endif  // PROFILER_H
