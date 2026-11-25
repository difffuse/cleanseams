// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/export.hpp"
#include <nlohmann/json.hpp>
// Source core/verdict.py: evidence-based layout verdicts, comparisons and the claim
// aggregate, over the same dictionaries (JSON objects) with Python semantics: None and
// NaN are distinct (both "missing"), booleans are not measurements, integers stay
// integers, and reason strings use Python repr/str formatting. Inputs the source
// rejects with a TypeError/AttributeError/IndexError raise ContractError here.
// JSON objects are unordered here; a container-valued measurement is printed with
// sorted keys where Python would keep insertion order.
namespace cleanseams::verdicts {
CLEANSEAMS_CORE_EXPORT bool finite(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT nlohmann::json measurements(const nlohmann::json& run);
CLEANSEAMS_CORE_EXPORT bool rejected(const nlohmann::json& run);
CLEANSEAMS_CORE_EXPORT nlohmann::json verdict(const nlohmann::json& run,const nlohmann::json& meta=nullptr);
CLEANSEAMS_CORE_EXPORT nlohmann::json compare(const nlohmann::json& a,const nlohmann::json& b,const nlohmann::json& meta=nullptr);
CLEANSEAMS_CORE_EXPORT nlohmann::json aggregate(const nlohmann::json& class_runs);
// Python repr()/str() of a JSON value (floats via the shortest round-trip repr).
CLEANSEAMS_CORE_EXPORT std::string python_repr(const nlohmann::json&);
CLEANSEAMS_CORE_EXPORT std::string python_str(const nlohmann::json&);
}
