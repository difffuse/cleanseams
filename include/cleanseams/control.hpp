// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/types.hpp"
#include <atomic>
#include <chrono>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
namespace cleanseams {
using MonotonicClock=std::chrono::steady_clock;
// The time source a deadline reads (MonotonicClock::now unless injected).
using ClockSource=std::function<MonotonicClock::time_point()>;
// Source control.expired/deadline_at: an optional monotonic instant. Empty never expires.
class CLEANSEAMS_CORE_EXPORT Deadline {
public:
    Deadline()=default;
    explicit Deadline(MonotonicClock::time_point at):at_(at) {}
    // A deadline read from an injected clock, shared by its copies (deterministic tests of
    // expiry partway through a stage; the source's control.expired reads time.monotonic).
    Deadline(MonotonicClock::time_point at,std::shared_ptr<const ClockSource> clock):at_(at),clock_(std::move(clock)) {}
    // Finite, nonnegative seconds after `start` (now by default); nullopt means none.
    static Deadline after(std::optional<double> seconds,std::optional<MonotonicClock::time_point> start=std::nullopt);
    static Deadline expired_now() {return Deadline(MonotonicClock::time_point::min());}
    // Expires at the (checks + 1)-th expired() call on it or its copies: a clock that
    // advances one tick per reading, from tick 0, against an instant `checks` ticks away.
    static Deadline after_checks(std::uint64_t checks);
    bool expired() const {return at_ && (clock_?(*clock_)():MonotonicClock::now())>=*at_;}
    const std::optional<MonotonicClock::time_point>& at() const {return at_;}
private:
    std::optional<MonotonicClock::time_point> at_;
    std::shared_ptr<const ClockSource> clock_;
};
// Cooperative cancellation shared between a session owner and its workers.
class CLEANSEAMS_CORE_EXPORT CancellationToken {
public:
    CancellationToken():flag_(std::make_shared<std::atomic<bool>>(false)) {}
    void request() const {flag_->store(true,std::memory_order_relaxed);}
    bool requested() const {return flag_->load(std::memory_order_relaxed);}
private:
    std::shared_ptr<std::atomic<bool>> flag_;
};
// Source Control: check() throws a `cancelled` ContractError when cancellation was
// requested, otherwise reports whether the deadline has passed.
struct CLEANSEAMS_CORE_EXPORT Control {
    CancellationToken cancel;
    Deadline deadline;
    bool check() const;
};
// Throws ContractError(cancelled) when `token` was requested.
CLEANSEAMS_CORE_EXPORT void throw_if_cancelled(const CancellationToken& token);
// ------------------------------------------------------------ native scheduling (T24)
// Bounded work batches and native evaluation threads. Every thread count gives the same
// results: only independent charts of one immutable snapshot are computed in parallel,
// and their results are committed in the original order. threads=1 is the reference.
struct SchedulingOptions {
    std::uint64_t threads=1;                 // native evaluation threads, the calling thread included
    std::uint64_t batch_max_charts=16;       // charts per evaluation or host-solve batch
    std::uint64_t batch_max_triangles=50000; // triangles per batch; a larger chart is a batch of its own
};
inline constexpr std::uint64_t max_native_threads=64;
// GB performance: the product default thread count (CLI without --threads, a session's
// automatic worker_count 0, the Blender extension): the machine's hardware threads, at
// most default_native_thread_cap, at least 1. Every count gives the same results; the
// library-level SchedulingOptions keep the reference serial default of one thread.
inline constexpr std::uint64_t default_native_thread_cap=8;
CLEANSEAMS_CORE_EXPORT std::uint64_t default_native_threads();
// Checks the bounds (1..64 threads, positive batch limits) and returns the options.
CLEANSEAMS_CORE_EXPORT SchedulingOptions resolve_scheduling(SchedulingOptions);
// A bounded set of native worker threads owned by one run. Workers run only native
// tasks on owned data; they never call the host. Nested solver parallelism is off:
// a task that asks for parallel work runs it inline, SuperLU is the serial build with
// its internal serial BLAS, and Eigen runs single-threaded.
class CLEANSEAMS_CORE_EXPORT WorkerPool {
public:
    explicit WorkerPool(std::uint64_t threads);
    ~WorkerPool(); // stops and joins every worker
    WorkerPool(const WorkerPool&)=delete;
    WorkerPool& operator=(const WorkerPool&)=delete;
    std::uint64_t threads() const;
    // Runs task(i) once for every i in [0, n) on the calling thread and the workers, and
    // returns when every started task has finished. Tasks write their own result slots;
    // the caller commits them in index order, so completion order never matters. The
    // failure of the lowest failing index is rethrown. Once `cancel` is requested, tasks
    // that have not started are skipped (the caller then discards the whole call).
    void for_each(std::size_t n,const std::function<void(std::size_t)>& task,const CancellationToken* cancel=nullptr);
    // GB performance: like for_each, but also from inside a worker task: the calling thread
    // runs every index no idle worker claims, idle workers help, and the call returns when
    // every index finished. The failure of the lowest failing index is rethrown.
    void for_each_nested(std::size_t n,const std::function<void(std::size_t)>& task);
    struct Stats {
        std::uint64_t parallel_calls=0,inline_calls=0,tasks=0,skipped=0,max_width=0;
        std::uint64_t nested_solver_threads=1; // largest solver/Eigen thread count seen inside a task
    };
    Stats stats() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
// Makes `pool` the calling thread's active pool for the scope (nesting restores).
class CLEANSEAMS_CORE_EXPORT WorkerScope {
public:
    explicit WorkerScope(WorkerPool* pool);
    ~WorkerScope();
    WorkerScope(const WorkerScope&)=delete;
    WorkerScope& operator=(const WorkerScope&)=delete;
private:
    WorkerPool* previous_=nullptr;
};
CLEANSEAMS_CORE_EXPORT WorkerPool* active_workers();
// True inside a worker task (on a worker or on the calling thread while it helps).
CLEANSEAMS_CORE_EXPORT bool in_worker_task();
// task(i) for every i in [0, n): on the active pool, or inline in index order when there
// is none, when n < 2 or inside a worker task (no nested threads).
CLEANSEAMS_CORE_EXPORT void parallel_for(std::size_t n,const std::function<void(std::size_t)>& task,const CancellationToken* cancel=nullptr);
// GB performance: task(i) for every i in [0, n) with idle workers of the active pool helping,
// also inside a worker task (the caller runs whatever no idle worker takes); inline in index
// order without a pool. For a few independent pieces of one computation whose results do not
// depend on the thread that computes them (the two factorizations of a chart solve).
CLEANSEAMS_CORE_EXPORT void parallel_for_nested(std::size_t n,const std::function<void(std::size_t)>& task);
// Process-wide gauges for leak checks: worker threads alive and pools alive.
struct WorkerCounters {std::uint64_t live_threads=0,live_pools=0;};
CLEANSEAMS_CORE_EXPORT WorkerCounters worker_counters();
}
