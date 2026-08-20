// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Profiler.h"

#ifdef svmp_profiling_enabled

#include <algorithm>
#include <cstdio>
#include <set>
#include <vector>

#include "Array.h"
#include "Array3.h"
#include "CmMod.h"
#include "SimulationLogger.h"
#include "Vector.h"

namespace profiler {

namespace {
const char* const kOpNames[] = {
  "mat_det", "mat_dev", "mat_dyad_prod", "mat_id", "mat_inv",
  "mat_mul_vec", "mat_mul_return", "mat_mul_inplace",
  "mat_symm", "mat_symm_prod",
  "ten_asym_prod12", "ten_dyad_prod", "ten_mddot", "ten_symm_prod",
  "transpose",
};
static_assert(sizeof(kOpNames) / sizeof(kOpNames[0]) ==
              static_cast<unsigned>(Op::COUNT),
              "kOpNames must stay in step with enum class Op");
}  // namespace

const char* op_name(Op op) { return kOpNames[static_cast<unsigned>(op)]; }

/// @brief Zero-initialised at load time; keys[] of 0 never matches a real key
/// because a valid shape always has a non-zero result-row field.
Registry registry{};

void OpTable::slow_add(unsigned key)
{
  for (int i = 0; i < n_used; i++) {
    if (keys[i] == key) { calls[i] += 1; last = i; return; }
  }
  if (n_used == CAP) { overflow_calls += 1; return; }
  keys[n_used] = key;
  calls[n_used] = 1;
  last = n_used;
  n_used += 1;
}

void unpack_shape(unsigned key, int& m, int& k, int& n)
{
  auto field = [](unsigned v) -> int {
    if (v == DIM_NA) return -1;
    if (v == DIM_SATURATED) return -2;      // reported as ">1022"
    return static_cast<int>(v);
  };
  m = field((key >> 20) & 0x3FF);
  k = field((key >> 10) & 0x3FF);
  n = field(key & 0x3FF);
}

void Registry::add_phase(const char* name, double seconds)
{
  auto& record = phases[name];
  if (record.calls == 0) { record.order = phase_order++; }
  record.seconds += seconds;
  record.calls += 1;
}

namespace {


/// @brief Pack a set of labels into a '\0'-separated buffer for MPI transport.
std::vector<char> pack(const std::vector<std::string>& labels)
{
  std::vector<char> buffer;
  for (const auto& label : labels) {
    buffer.insert(buffer.end(), label.begin(), label.end());
    buffer.push_back('\0');
  }
  return buffer;
}

std::vector<std::string> unpack(const std::vector<char>& buffer)
{
  std::vector<std::string> labels;
  std::string current;
  for (char c : buffer) {
    if (c == '\0') { labels.push_back(current); current.clear(); }
    else { current.push_back(c); }
  }
  return labels;
}

/// @brief Union a rank-local label set across all ranks.
///
/// Ranks are not assumed to hold the same labels: a partition may contain
/// element types another does not, giving it shape keys the others never see.
/// Every rank ends up with the same ordered vector.
std::vector<std::string> union_labels(const std::vector<std::string>& local, MPI_Comm comm)
{
  int nprocs = 1;
  MPI_Comm_size(comm, &nprocs);
  if (nprocs == 1) { return local; }

  auto send = pack(local);
  int send_size = static_cast<int>(send.size());

  std::vector<int> sizes(nprocs);
  MPI_Allgather(&send_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);

  std::vector<int> displs(nprocs, 0);
  int total = 0;
  for (int i = 0; i < nprocs; i++) { displs[i] = total; total += sizes[i]; }

  std::vector<char> recv(total);
  MPI_Allgatherv(send.data(), send_size, MPI_CHAR,
                 recv.data(), sizes.data(), displs.data(), MPI_CHAR, comm);

  auto all = unpack(recv);
  std::set<std::string> unique(all.begin(), all.end());
  return std::vector<std::string>(unique.begin(), unique.end());
}

/// @brief Reduced statistics for one quantity across ranks.
struct Reduced {
  double sum{0.0}, min{0.0}, max{0.0};
  double mean(int nprocs) const { return sum / nprocs; }
  double imbalance(int nprocs) const
  {
    const double m = mean(nprocs);
    return (m > 0.0) ? max / m : 1.0;
  }
};

Reduced reduce_value(double local, MPI_Comm comm)
{
  Reduced r;
  MPI_Allreduce(&local, &r.sum, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(&local, &r.min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&local, &r.max, 1, MPI_DOUBLE, MPI_MAX, comm);
  return r;
}

}  // namespace

void Registry::report(const CmMod& cm_mod, const cmType& cm, SimulationLogger& logger,
                      const std::string& output_dir)
{
  const MPI_Comm comm = cm.com();
  const int nprocs = cm.np();
  const bool master = cm.mas(cm_mod);
  const std::string prefix = output_dir.empty() ? std::string("") : output_dir + "/";

  // ---------------------------------------------------------------- phases
  std::vector<std::string> phase_names;
  for (const auto& entry : phases) { phase_names.push_back(entry.first); }
  phase_names = union_labels(phase_names, comm);

  // Report in first-seen order on the master rank, which follows the order the
  // phases occur in the time loop; ranks that never entered a phase sort last.
  std::stable_sort(phase_names.begin(), phase_names.end(),
                   [this](const std::string& a, const std::string& b) {
                     auto ia = phases.find(a), ib = phases.find(b);
                     int oa = (ia == phases.end()) ? 1 << 30 : ia->second.order;
                     int ob = (ib == phases.end()) ? 1 << 30 : ib->second.order;
                     return oa < ob;
                   });

  struct PhaseOut { std::string name; Reduced sec; long calls; };
  std::vector<PhaseOut> phase_out;
  double total_time = 0.0;

  for (const auto& name : phase_names) {
    auto it = phases.find(name);
    const double local = (it == phases.end()) ? 0.0 : it->second.seconds;
    const long local_calls = (it == phases.end()) ? 0 : it->second.calls;

    Reduced sec = reduce_value(local, comm);
    long calls = 0;
    MPI_Allreduce(&local_calls, &calls, 1, MPI_LONG, MPI_SUM, comm);

    phase_out.push_back({name, sec, calls});
    total_time = std::max(total_time, sec.max);
  }

  // ---------------------------------------------------------------- ops
  //
  // Shape keys are data dependent: a partition holding a different element type
  // sees shapes the others never do, so the key sets are unioned rather than
  // assumed identical. Every rank derives the same canonical ordering by
  // sorting the concatenation, so no further communication is needed.
  std::vector<std::uint64_t> local_keys;
  for (unsigned op = 0; op < static_cast<unsigned>(Op::COUNT); op++) {
    const OpTable& table = ops[op];
    for (int i = 0; i < table.n_used; i++) {
      local_keys.push_back((static_cast<std::uint64_t>(op) << 32) | table.keys[i]);
    }
  }

  std::vector<std::uint64_t> union_keys = local_keys;
  if (nprocs > 1) {
    int local_n = static_cast<int>(local_keys.size());
    std::vector<int> counts(nprocs);
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

    std::vector<int> displs(nprocs, 0);
    int total = 0;
    for (int i = 0; i < nprocs; i++) { displs[i] = total; total += counts[i]; }

    union_keys.assign(total, 0);
    MPI_Allgatherv(local_keys.data(), local_n, MPI_UINT64_T,
                   union_keys.data(), counts.data(), displs.data(), MPI_UINT64_T, comm);
  }
  std::sort(union_keys.begin(), union_keys.end());
  union_keys.erase(std::unique(union_keys.begin(), union_keys.end()), union_keys.end());

  struct OpOut { std::string name; int m, k, n; double calls; double imbalance; };
  std::vector<OpOut> op_out;

  for (std::uint64_t gkey : union_keys) {
    const unsigned op = static_cast<unsigned>(gkey >> 32);
    const unsigned key = static_cast<unsigned>(gkey & 0xFFFFFFFFu);

    double local = 0.0;
    const OpTable& table = ops[op];
    for (int i = 0; i < table.n_used; i++) {
      if (table.keys[i] == key) { local = static_cast<double>(table.calls[i]); break; }
    }

    // Counts go through as double: exact to 2^53, well past any call count, and
    // it sidesteps the int/double-only dispatch in cmType::reduce.
    double sum = local, max = local;
    if (nprocs > 1) {
      MPI_Allreduce(&local, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(&local, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
    }

    int m, k, n;
    unpack_shape(key, m, k, n);
    const double mean = sum / nprocs;
    op_out.push_back({op_name(static_cast<Op>(op)), m, k, n, sum,
                      (mean > 0.0) ? max / mean : 1.0});
  }

  std::sort(op_out.begin(), op_out.end(),
            [](const OpOut& a, const OpOut& b) { return a.calls > b.calls; });

  if (!master) { return; }

  // ---------------------------------------------------------------- report
  // M = rows of result, K = contracted dimension ("-" when there is none),
  // N = cols of result. ">1022" marks a dimension that saturated the key.
  auto dim = [](int v) {
    if (v == -1) return std::string("-");
    if (v == -2) return std::string(">1022");
    return std::to_string(v);
  };
  auto shape_string = [&dim](int m, int k, int n) {
    return dim(m) + " " + dim(k) + " " + dim(n);
  };

  char line[256];
  const std::string sep(78, '=');
  const std::string dash(78, '-');

  logger << std::endl << sep << std::endl;
  logger << " Profiling summary -- " << nprocs << " rank(s)" << std::endl;
  logger << sep << std::endl;

  logger << " Phase timings" << std::endl << dash << std::endl;
  snprintf(line, sizeof(line), "  %-22s %10s %8s %10s %10s %8s %10s",
           "phase", "mean(s)", "%", "min(s)", "max(s)", "imbal", "calls");
  logger << line << std::endl;
  for (const auto& p : phase_out) {
    snprintf(line, sizeof(line), "  %-22s %10.3f %7.1f%% %10.3f %10.3f %8.2f %10ld",
             p.name.c_str(), p.sec.mean(nprocs),
             (total_time > 0.0) ? 100.0 * p.sec.mean(nprocs) / total_time : 0.0,
             p.sec.min, p.sec.max, p.sec.imbalance(nprocs), p.calls);
    logger << line << std::endl;
  }

  logger << std::endl << " Operation counts by operand shape" << std::endl << dash << std::endl;
  snprintf(line, sizeof(line), "  %-24s %-16s %16s %8s", "operation", "shape", "calls", "imbal");
  logger << line << std::endl;
  for (const auto& o : op_out) {
    snprintf(line, sizeof(line), "  %-24s %-16s %16.0f %8.2f",
             o.name.c_str(), shape_string(o.m, o.k, o.n).c_str(), o.calls, o.imbalance);
    logger << line << std::endl;
  }

  // ---------------------------------------------------------------- CSV
  // The project has no CSV precedent, but the analysis scripts consume it and
  // the table above stays the human-facing output.
  if (FILE* f = fopen((prefix + "profile_phases.csv").c_str(), "w")) {
    fprintf(f, "phase,mean_s,pct,min_s,max_s,imbalance,calls,nranks\n");
    for (const auto& p : phase_out) {
      fprintf(f, "%s,%.9g,%.9g,%.9g,%.9g,%.9g,%ld,%d\n", p.name.c_str(),
              p.sec.mean(nprocs),
              (total_time > 0.0) ? 100.0 * p.sec.mean(nprocs) / total_time : 0.0,
              p.sec.min, p.sec.max, p.sec.imbalance(nprocs), p.calls, nprocs);
    }
    fclose(f);
  }

  if (FILE* f = fopen((prefix + "profile_ops.csv").c_str(), "w")) {
    fprintf(f, "operation,m,k,n,calls,imbalance,nranks\n");
    for (const auto& o : op_out) {
      fprintf(f, "%s,%d,%d,%d,%.0f,%.9g,%d\n",
              o.name.c_str(), o.m, o.k, o.n, o.calls, o.imbalance, nprocs);
    }
    fclose(f);
  }

  // ---------------------------------------------------------------- memory
  // Read straight from the container statics rather than through
  // Array<T>::memory()/stats(), which are only defined for some instantiations.
  struct MemRow { const char* name; long allocated; double in_use_kb; double returned_kb; };
  const MemRow mem[] = {
    {"Array<double>",  Array<double>::num_allocated,
     Array<double>::memory_in_use / 1024.0,  Array<double>::memory_returned / 1024.0},
    {"Array<int>",     Array<int>::num_allocated,
     Array<int>::memory_in_use / 1024.0,     Array<int>::memory_returned / 1024.0},
    {"Array3<double>", Array3<double>::num_allocated,
     Array3<double>::memory_in_use / 1024.0, Array3<double>::memory_returned / 1024.0},
    {"Vector<double>", Vector<double>::num_allocated,
     Vector<double>::memory_in_use / 1024.0, Vector<double>::memory_returned / 1024.0},
    {"Vector<int>",    Vector<int>::num_allocated,
     Vector<int>::memory_in_use / 1024.0,    Vector<int>::memory_returned / 1024.0},
  };

  logger << std::endl << " Container allocations (master rank)" << std::endl << dash << std::endl;
  snprintf(line, sizeof(line), "  %-18s %16s %16s %16s",
           "container", "allocations", "in use (KB)", "returned (KB)");
  logger << line << std::endl;
  for (const auto& m : mem) {
    snprintf(line, sizeof(line), "  %-18s %16ld %16.1f %16.1f",
             m.name, m.allocated, m.in_use_kb, m.returned_kb);
    logger << line << std::endl;
  }

  if (FILE* f = fopen((prefix + "profile_memory.csv").c_str(), "w")) {
    fprintf(f, "container,allocations,in_use_kb,returned_kb,nranks\n");
    for (const auto& m : mem) {
      fprintf(f, "%s,%ld,%.9g,%.9g,%d\n",
              m.name, m.allocated, m.in_use_kb, m.returned_kb, nprocs);
    }
    fclose(f);
  }

  logger << dash << std::endl;
  logger << " Profiling CSV -> profile_phases.csv, profile_ops.csv, profile_memory.csv"
         << std::endl << sep << std::endl;
}

}  // namespace profiler

#endif  // svmp_profiling_enabled
