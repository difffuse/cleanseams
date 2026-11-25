// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/control.hpp"
#include "cleanseams/export.hpp"
#include <chrono>
#include <string>
#include <string_view>
#include <vector>
namespace cleanseams {
struct ProcessResult {
    int exit_code=-1;          // exit status, or 128 + signal
    std::string out,err;       // captured, each bounded by `max_output` bytes
    bool timed_out=false,out_truncated=false,err_truncated=false;
    bool cancelled=false;
};
struct ProcessOptions {
    std::chrono::milliseconds timeout{120000};
    std::size_t max_output=16u<<20;
    CancellationToken cancel;
};
// Runs argv directly (UTF-8 on Windows), never through a shell. Captures bounded
// stdout/stderr while writing stdin. Owns ordinary descendants (POSIX process group /
// Windows job); terminates them on direct-child exit, timeout or cancellation.
// An already requested cancellation prevents launch. Embedded NUL arguments are invalid.
CLEANSEAMS_CORE_EXPORT ProcessResult run_process(const std::vector<std::string>& argv,std::string_view input={},const ProcessOptions& options={});
}
