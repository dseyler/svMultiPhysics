// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PROFILER_H
#define PROFILER_H

/// @file Profiler.h
/// @brief Optional profiling: phase timings and container allocation totals.
///
/// Compiled out entirely unless configured with -DENABLE_PROFILING=ON.
///
/// Usage:
/// @code
///   SVMP_PROFILE_PHASE("newton_assembly");   // times this scope
/// @endcode
///
/// Phases wrap stage functions and fire O(1e4) times a run, so a string-keyed
/// map costs nothing next to the work they measure.
///
/// Never place SVMP_PROFILE_PHASE inside an element or Gauss point loop; two
/// clock reads would swamp the body, and a string hash per iteration would
/// distort the very timings being gathered.
///
/// @note Guard on #ifdef, never #if. These macros are defined without a value,
/// and #if on an undefined macro evaluates to 0, silently dropping the guarded
/// code with no diagnostic.

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

struct PhaseRecord {
  double        seconds{0.0};
  std::uint64_t calls{0};
  int           order{0};         ///< first-seen, so the report keeps source order
};

struct Registry {
  std::unordered_map<std::string, PhaseRecord> phases;
  int phase_order{0};

  void add_phase(const char* name, double seconds);

  /// @brief Reduce across ranks and emit the report. Collective: every rank
  /// must call it. Only the master rank writes.
  void report(const CmMod& cm_mod, const cmType& cm, SimulationLogger& logger,
              const std::string& output_dir);
};

extern Registry registry;

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
#pragma message("svmp profiler phase records are not thread safe; add atomics before parallelising instrumented regions.")
#endif

#define SVMP_PROFILE_CONCAT_(a, b) a##b
#define SVMP_PROFILE_CONCAT(a, b) SVMP_PROFILE_CONCAT_(a, b)

#define SVMP_PROFILE_PHASE(name) \
  ::profiler::ScopedPhase SVMP_PROFILE_CONCAT(svmp_phase_, __LINE__)(name)

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

#endif  // svmp_profiling_enabled

#endif  // PROFILER_H
