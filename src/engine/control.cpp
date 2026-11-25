// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/control.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
namespace cleanseams {
Deadline Deadline::after(std::optional<double> seconds,std::optional<MonotonicClock::time_point> start) {
    if(!seconds)return Deadline();
    if(!std::isfinite(*seconds) || *seconds<0)
        throw ContractError(ErrorCode::invalid_option,"deadline_s","deadline_s must be finite and non-negative, or None");
    const auto base=start.value_or(MonotonicClock::now());
    const auto limit=std::chrono::duration<double>(MonotonicClock::time_point::max()-base).count();
    if(*seconds>=limit)return Deadline(MonotonicClock::time_point::max());
    return Deadline(base+std::chrono::duration_cast<MonotonicClock::duration>(std::chrono::duration<double>(*seconds)));
}
Deadline Deadline::after_checks(std::uint64_t checks) {
    auto ticks=std::make_shared<std::atomic<std::int64_t>>(0);
    auto clock=std::make_shared<const ClockSource>([ticks] {return MonotonicClock::time_point(MonotonicClock::duration(ticks->fetch_add(1)));});
    return Deadline(MonotonicClock::time_point(MonotonicClock::duration(static_cast<std::int64_t>(checks))),std::move(clock));
}
bool Control::check() const {
    throw_if_cancelled(cancel);
    return deadline.expired();
}
void throw_if_cancelled(const CancellationToken& token) {
    if(token.requested())throw ContractError(ErrorCode::cancelled,"control","operation cancelled");
}
std::uint64_t default_native_threads() {
    const std::uint64_t hardware=std::thread::hardware_concurrency();
    return std::clamp<std::uint64_t>(hardware,1,default_native_thread_cap);
}
SchedulingOptions resolve_scheduling(SchedulingOptions o) {
    if(o.threads<1 || o.threads>max_native_threads)
        throw ContractError(ErrorCode::invalid_option,"threads","native threads must be between 1 and "+std::to_string(max_native_threads));
    if(o.batch_max_charts<1)throw ContractError(ErrorCode::invalid_option,"batch_max_charts","batch_max_charts must be positive");
    if(o.batch_max_triangles<1)throw ContractError(ErrorCode::invalid_option,"batch_max_triangles","batch_max_triangles must be positive");
    return o;
}
namespace {
std::atomic<std::uint64_t> live_threads{0},live_pools{0};
thread_local WorkerPool* active_pool=nullptr;
thread_local int task_depth=0; // > 0 while this thread runs a worker task
// Solver threads a task would use: SuperLU is the serial build (internal serial BLAS);
// Eigen is not built with OpenMP, so nbThreads() is 1. Measured, not assumed.
std::uint64_t solver_threads() {return static_cast<std::uint64_t>(std::max(1,Eigen::nbThreads()));}
void single_threaded_solvers() {
    static std::once_flag once;
    std::call_once(once,[] {Eigen::initParallel();Eigen::setNbThreads(1);});
}
}
struct WorkerPool::Impl {
    struct Job {
        std::size_t n=0;
        const std::function<void(std::size_t)>* task=nullptr;
        const CancellationToken* cancel=nullptr;
        TopologyCache* cache=nullptr; // the caller's active cache, shared by the workers
        std::uint64_t memory_budget=default_memory_budget; // the caller's run allowance (Gate A F3)
        std::atomic<std::size_t> next{0},done{0};
        std::size_t active=0; // workers inside this job (guarded by the pool mutex)
        std::vector<std::exception_ptr> errors;
        std::atomic<std::uint64_t> skipped{0},nested{1};
    };
    WorkerPool* owner=nullptr;
    std::uint64_t threads=1;
    std::vector<std::thread> workers;
    mutable std::mutex mutex;
    std::condition_variable wake,finished;
    // Open jobs, innermost last (GB performance: a task may open a nested job that idle
    // workers help with). Guarded by the pool mutex.
    std::vector<Job*> jobs;
    bool stop=false;
    std::mutex caller; // one top-level for_each at a time; a second caller runs inline
    Stats stats;
    // Claims indices until none is left. Called by workers and by the caller.
    static void work(Job& j) {
        CacheScope scope(j.cache);MemoryBudgetScope budget(j.memory_budget);
        ++task_depth;
        for(;;) {
            const auto i=j.next.fetch_add(1);
            if(i>=j.n)break;
            if(j.cancel && j.cancel->requested())++j.skipped;
            else {
                auto seen=j.nested.load();const auto now=solver_threads();
                while(now>seen && !j.nested.compare_exchange_weak(seen,now)) {}
                try {(*j.task)(i);} catch(...) {j.errors[i]=std::current_exception();}
            }
            j.done.fetch_add(1);
        }
        --task_depth;
    }
    // The innermost open job with an unclaimed index (the pool mutex held).
    Job* claimable() const {
        for(auto it=jobs.rbegin();it!=jobs.rend();++it)if((*it)->next.load()<(*it)->n)return *it;
        return nullptr;
    }
    void run_worker() {
        ++live_threads;
        active_pool=owner; // a task may open a nested job on this pool
        for(;;) {
            Job* j=nullptr;
            {
                std::unique_lock lock(mutex);
                wake.wait(lock,[&] {return stop || (j=claimable())!=nullptr;});
                if(stop)break;
                ++j->active;
            }
            work(*j);
            {
                std::lock_guard lock(mutex);
                if(--j->active==0)finished.notify_all();
            }
        }
        --live_threads;
    }
    // Runs the job on the calling thread with the help of idle workers and returns once
    // every task finished and every worker left it (then the job may be destroyed).
    void run(Job& j) {
        {std::lock_guard lock(mutex);jobs.push_back(&j);}
        wake.notify_all();
        work(j);
        std::unique_lock lock(mutex);
        finished.wait(lock,[&] {return j.done.load()==j.n && j.active==0;});
        jobs.erase(std::find(jobs.begin(),jobs.end(),&j));
        ++stats.parallel_calls;stats.tasks+=j.n-j.skipped.load();stats.skipped+=j.skipped.load();
        stats.max_width=std::max<std::uint64_t>(stats.max_width,std::min<std::uint64_t>(j.n,threads));
        stats.nested_solver_threads=std::max(stats.nested_solver_threads,j.nested.load());
    }
};
WorkerPool::WorkerPool(std::uint64_t threads):impl_(std::make_unique<Impl>()) {
    if(threads<1 || threads>max_native_threads)
        throw ContractError(ErrorCode::invalid_option,"threads","native threads must be between 1 and "+std::to_string(max_native_threads));
    single_threaded_solvers();
    impl_->threads=threads;impl_->owner=this;
    ++live_pools;
    try {
        for(std::uint64_t k=1;k<threads;++k)impl_->workers.emplace_back([impl=impl_.get()] {impl->run_worker();});
    } catch(...) {
        {std::lock_guard lock(impl_->mutex);impl_->stop=true;}
        impl_->wake.notify_all();
        for(auto& w:impl_->workers)w.join();
        --live_pools;throw;
    }
}
WorkerPool::~WorkerPool() {
    {std::lock_guard lock(impl_->mutex);impl_->stop=true;}
    impl_->wake.notify_all();
    for(auto& w:impl_->workers)if(w.joinable())w.join();
    --live_pools;
}
std::uint64_t WorkerPool::threads() const {return impl_->threads;}
WorkerPool::Stats WorkerPool::stats() const {std::lock_guard lock(impl_->mutex);return impl_->stats;}
void WorkerPool::for_each(std::size_t n,const std::function<void(std::size_t)>& task,const CancellationToken* cancel) {
    if(n==0)return;
    auto& s=*impl_;
    // Nested requests (inside a task) run inline before the owner lock is touched: the
    // calling thread may already hold it. A second concurrent caller also runs inline.
    std::unique_lock<std::mutex> owner;
    if(n>=2 && !s.workers.empty() && task_depth==0)owner=std::unique_lock(s.caller,std::try_to_lock);
    if(!owner.owns_lock()) {
        // Inline, in index order: exactly the serial schedule.
        {std::lock_guard lock(s.mutex);++s.stats.inline_calls;s.stats.max_width=std::max<std::uint64_t>(s.stats.max_width,1);}
        std::uint64_t skipped=0,ran=0;
        for(std::size_t i=0;i<n;++i) {
            if(cancel && cancel->requested()) {++skipped;continue;}
            ++ran;task(i);
        }
        std::lock_guard lock(s.mutex);s.stats.tasks+=ran;s.stats.skipped+=skipped;
        return;
    }
    Impl::Job j;j.n=n;j.task=&task;j.cancel=cancel;j.cache=active_cache();j.memory_budget=active_memory_budget();j.errors.resize(n);
    s.run(j);
    for(auto& error:j.errors)if(error)std::rethrow_exception(error);
}
void WorkerPool::for_each_nested(std::size_t n,const std::function<void(std::size_t)>& task) {
    if(n==0)return;
    auto& s=*impl_;
    if(n<2 || s.workers.empty()) {for(std::size_t i=0;i<n;++i)task(i);return;}
    Impl::Job j;j.n=n;j.task=&task;j.cache=active_cache();j.memory_budget=active_memory_budget();j.errors.resize(n);
    s.run(j);
    for(auto& error:j.errors)if(error)std::rethrow_exception(error);
}
WorkerScope::WorkerScope(WorkerPool* pool):previous_(active_pool) {active_pool=pool;}
WorkerScope::~WorkerScope() {active_pool=previous_;}
WorkerPool* active_workers() {return active_pool;}
bool in_worker_task() {return task_depth>0;}
void parallel_for(std::size_t n,const std::function<void(std::size_t)>& task,const CancellationToken* cancel) {
    if(active_pool && task_depth==0 && n>1) {active_pool->for_each(n,task,cancel);return;}
    for(std::size_t i=0;i<n;++i) {
        if(cancel && cancel->requested())continue;
        task(i);
    }
}
void parallel_for_nested(std::size_t n,const std::function<void(std::size_t)>& task) {
    if(active_pool && n>1) {active_pool->for_each_nested(n,task);return;}
    for(std::size_t i=0;i<n;++i)task(i);
}
WorkerCounters worker_counters() {return {live_threads.load(),live_pools.load()};}
}
