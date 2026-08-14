// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <chrono>
#include <future>

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

namespace {
inline void myCheckInterruptFn(void * /*dummy*/) { R_CheckUserInterrupt(); }
} // anonymous namespace

// Check for interrupts and throw the sentinel exception if one is pending
inline bool hasUserInterrupt() { return (R_ToplevelExec(myCheckInterruptFn, NULL) == FALSE); }

// Wrap a function call such that we will check for R user interrupts
// The function itself is run in a background thread, and the main thread checks for an
// R user interrupt every 100ms.
// The function must take a pointer to as std::atomic<bool> as the last argument, and exit
// early if the bool gets set to true
// NOTE: It is EXTREMELY IMPORTANT that no R objects are created/destroyed inside the spawned
// thread, which includes destructors that mess with R's GC protection.
template <class F, class... Args>
std::invoke_result_t<F, Args..., std::atomic<bool> *>
run_with_R_interrupt_check(F &&f, Args &&...args) {
    std::atomic<bool> interrupt(false);
    auto job =
        std::async(std::launch::async, std::forward<F>(f), std::forward<Args>(args)..., &interrupt);
    while (job.wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
        if (hasUserInterrupt()) {
            interrupt = true;
        }
    }
    if (interrupt) {
        throw Rcpp::internal::InterruptedException();
    }
    return job.get();
}

// Wrap and parallelize a function call while checking for R interrupts
// Arguments:
//  Function f:
//     parameters are: (int) task_start, (int) task_end, and pointer to std::atomic<bool>
//     Will run tasks with IDs [task_start, task_end).
//     Return value if any is ignored.
//     Should terminate early if the std::atomic<bool> becomes true
//  int n_tasks:
//     Number of tasks to run
//  int n_threads:
//     Number of threads to run with
// Parallelization assigns a contigous range of tasks to each thread
template <class F>
void run_parallel_with_R_interrupt_check(F &&f, size_t n_tasks, size_t n_threads) {
    std::vector<std::future<void>> task_vec;
    std::atomic<bool> interrupt = false;
    size_t idx = 0;
    n_threads = std::max<size_t>(1, n_threads);

    for (size_t i = 0; i < n_threads; i++) {
        size_t items = (n_tasks - idx) / (n_threads - i);
        task_vec.push_back(std::async(std::launch::async, std::forward<F>(f), idx, idx+items, &interrupt));
        idx += items;
    }
    // Wait for work to finish while checking for interrupts
    for (size_t i = 0; i < n_threads; i++) {
        while(task_vec[i].wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
            if (hasUserInterrupt()) {
                interrupt = true;
            }
        }
    }
    if (interrupt) throw Rcpp::internal::InterruptedException();
}