// SPDX-License-Identifier: GPL-3.0-or-later
#include "arrays.hpp"
#include "cleanseams/io.hpp"
#include <cmath>
namespace cs_bindings {
void require_array(const py::array& array,const py::dtype& dtype,std::size_t alignment,int rank,const char* name) {
    if(!array.dtype().is(dtype) || array.ndim()!=rank || !(array.flags() & py::array::c_style) ||
       reinterpret_cast<std::uintptr_t>(array.data())%alignment!=0)
        throw py::value_error(std::string(name)+": exact dtype, rank, alignment and C-contiguous storage required");
}
void require_width(const py::array& array,py::ssize_t width,const char* name) {
    if(array.ndim()!=2 || array.shape(1)!=width)
        throw py::value_error(std::string(name)+" must have shape (n,"+std::to_string(width)+")");
}
py::array dict_array(const py::dict& dict,const char* key) {
    if(!dict.contains(key) || !py::isinstance<py::array>(dict[key]))throw py::value_error(std::string(key)+": NumPy array required");
    return py::cast<py::array>(dict[key]);
}
std::vector<cleanseams::Vec2> copy_uv(const py::array& array,const char* name,bool finite) {
    require_array<double>(array,2,name);require_width(array,2,name);
    const auto e=elements<double>(array);std::vector<cleanseams::Vec2> out;out.reserve(e.size()/2);
    for(std::size_t i=0;i<e.size();i+=2) {
        if(finite && (!std::isfinite(e[i]) || !std::isfinite(e[i+1])))throw py::value_error(std::string(name)+" must be finite");
        out.push_back({e[i],e[i+1]});
    }
    return out;
}
namespace {
py::object field(const py::dict& dict,const char* key,const char* owner) {
    if(!dict.contains(key))throw py::value_error(std::string(owner)+"."+key+" is required");
    return dict[key];
}
}
std::uint64_t dict_uint(const py::dict& dict,const char* key,const char* owner) {
    const auto value=field(dict,key,owner);
    // bool is an int subclass in Python; a flag is never an identity.
    if(!py::isinstance<py::int_>(value) || py::isinstance<py::bool_>(value))
        throw py::value_error(std::string(owner)+"."+key+" must be an integer");
    try {return py::cast<std::uint64_t>(value);}
    catch(const py::cast_error&) {throw py::value_error(std::string(owner)+"."+key+" must be a non-negative 64-bit integer");}
}
bool dict_bool(const py::dict& dict,const char* key,const char* owner) {
    const auto value=field(dict,key,owner);
    if(!py::isinstance<py::bool_>(value))throw py::value_error(std::string(owner)+"."+key+" must be a bool");
    return py::cast<bool>(value);
}
std::string dict_str(const py::dict& dict,const char* key,const char* owner) {
    const auto value=field(dict,key,owner);
    if(!py::isinstance<py::str>(value))throw py::value_error(std::string(owner)+"."+key+" must be a string");
    return py::cast<std::string>(value);
}
std::vector<double> numbers(const py::handle& value,const char* name) {
    const auto error=std::string(name)+" must be a sequence of finite numbers";
    if(py::isinstance<py::str>(value) || py::isinstance<py::bytes>(value) || !py::isinstance<py::sequence>(value))throw py::value_error(error);
    std::vector<double> out;
    for(const auto& item:py::reinterpret_borrow<py::sequence>(value)) {
        if(py::isinstance<py::bool_>(item) || py::isinstance<py::str>(item) || py::isinstance<py::bytes>(item))throw py::value_error(error);
        double x=0;
        try {x=py::cast<double>(py::float_(item));} catch(const py::error_already_set&) {throw py::value_error(error);}
        if(!std::isfinite(x))throw py::value_error(error);
        out.push_back(x);
    }
    return out;
}
nlohmann::json parse_json(const std::string& text,const char* field) {return cleanseams::parse_untrusted_json(text,field);}
}
