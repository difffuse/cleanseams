// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/types.hpp"
#include <stdexcept>
#include <array>
#include <string_view>
#include <string>
namespace cleanseams {
enum class ErrorCode { invalid_option, invalid_input, unsupported_version, constraint_conflict,
    unsupported_topology, numerical_failure, stale_state, host_failure, cancelled, memory_limit, io_error };
inline constexpr std::array<std::string_view,11> error_codes{"invalid_option","invalid_input","unsupported_version","constraint_conflict","unsupported_topology","numerical_failure","stale_state","host_failure","cancelled","memory_limit","io_error"};
inline std::string_view error_code_name(ErrorCode code) {return error_codes.at(static_cast<std::size_t>(code));}
class CLEANSEAMS_CORE_EXPORT ContractError : public std::invalid_argument {
public:
    ErrorCode code;
    std::string field;
    ContractError(ErrorCode c, std::string path, std::string message)
      : std::invalid_argument(path + ": " + message), code(c), field(std::move(path)) {}
};
}
