// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Boundary checks of the CPython module (T25): every array, dictionary field and JSON
// string from Python is validated here, while the GIL is held, before any native work.
// Nothing returned by these helpers refers to Python memory after the call that made it.
#include "cleanseams/mesh.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <nlohmann/json.hpp>
#include <cstdint>
#include <cstring>
#include <span>
#include <string>
#include <vector>
namespace cs_bindings {
namespace py=pybind11;
// Exact dtype (no implicit conversion), rank, C-contiguous storage and element alignment.
void require_array(const py::array& array,const py::dtype& dtype,std::size_t alignment,int rank,const char* name);
template<class T> void require_array(const py::array& array,int rank,const char* name) {
    require_array(array,py::dtype::of<T>(),alignof(T),rank,name);
}
// Second dimension of a rank-2 array (rows of `width` values).
void require_width(const py::array& array,py::ssize_t width,const char* name);
template<class T> std::span<const T> elements(const py::array& array) {
    return {static_cast<const T*>(array.data()),static_cast<std::size_t>(array.size())};
}
template<class T> std::vector<T> copy_elements(const py::array& array) {const auto e=elements<T>(array);return {e.begin(),e.end()};}
// dict[key] as a NumPy array (ValueError otherwise, KeyError never escapes).
py::array dict_array(const py::dict& dict,const char* key);
// An (n,2) float64 map copied into owned storage; `finite` refuses NaN and infinities
// (pins). Incoming and host-solved UV maps may hold them: the engine decides (Gate A F6/F9).
std::vector<cleanseams::Vec2> copy_uv(const py::array& array,const char* name,bool finite);
// Typed scalar fields of a host dictionary: ValueError naming the field on any mismatch.
std::uint64_t dict_uint(const py::dict& dict,const char* key,const char* owner);
bool dict_bool(const py::dict& dict,const char* key,const char* owner);
std::string dict_str(const py::dict& dict,const char* key,const char* owner);
// A sequence of finite real numbers (list, tuple, mathutils vector, NumPy array); bool
// and str are refused.
std::vector<double> numbers(const py::handle& value,const char* name);
// JSON text from Python through the bounded native parser (ContractError on bad input).
nlohmann::json parse_json(const std::string& text,const char* field);
template<class T> py::array_t<T> owned_array(const std::vector<T>& values) {
    py::array_t<T> out(static_cast<py::ssize_t>(values.size()));
    if(!values.empty())std::memcpy(out.mutable_data(),values.data(),values.size()*sizeof(T));
    return out;
}
template<class T,std::size_t N> py::array_t<T> owned_array(const std::vector<std::array<T,N>>& values) {
    static_assert(sizeof(std::array<T,N>)==N*sizeof(T));
    py::array_t<T> out(std::vector<py::ssize_t>{static_cast<py::ssize_t>(values.size()),static_cast<py::ssize_t>(N)});
    if(!values.empty())std::memcpy(out.mutable_data(),values.data(),values.size()*N*sizeof(T));
    return out;
}
}
