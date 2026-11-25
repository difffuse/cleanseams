// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/errors.hpp"
#include <nlohmann/json.hpp>
#include <limits>
#include <optional>
#include <type_traits>
namespace cleanseams::detail {
using Json=nlohmann::json;
template<class T> struct optional_type : std::false_type {};
template<class T> struct optional_type<std::optional<T>> : std::true_type { using value_type=T; };
template<class T> concept Reflected=requires(T& t) {visit_fields(t,[](auto,auto&){});};
template<class T> concept Sequence=requires(T& t){t.begin();t.end();typename T::value_type;} && !std::is_same_v<T,std::string> && !std::is_same_v<T,Json>;
template<class T> Json encode(const T& value) {
    if constexpr(optional_type<T>::value) return value ? encode(*value) : Json(nullptr);
    else if constexpr(Reflected<T>) { Json result=Json::object(); visit_fields(value,[&](auto name,const auto& v){result[name]=encode(v);}); return result; }
    else if constexpr(Sequence<T>) { Json result=Json::array();for(const auto& v:value) result.push_back(encode(v));return result; }
    else return value;
}
template<class T> T decode(const Json& value, const std::string& path) {
    if constexpr(optional_type<T>::value) {if(value.is_null()) return std::nullopt;return decode<typename T::value_type>(value,path);}
    else if constexpr(Reflected<T>) {T result;visit_fields(result,[&](auto name,auto& v){using V=std::remove_cvref_t<decltype(v)>;v=decode<V>(value.at(name),path+"."+name);});return result;}
    else if constexpr(Sequence<T>) {
        T result{};
        if constexpr(requires{result.resize(value.size());}) result.resize(value.size());
        if(result.size()!=value.size()) throw ContractError(ErrorCode::invalid_input,path,"wrong array length");
        for(std::size_t i=0;i<value.size();++i) result[i]=decode<typename T::value_type>(value[i],path+"["+std::to_string(i)+"]");
        return result;
    } else if constexpr(std::is_integral_v<T> && !std::is_same_v<T,bool>) {
        if(value.is_number_unsigned()) {
            if(value.get<std::uint64_t>()>static_cast<std::uint64_t>(std::numeric_limits<T>::max()))
                throw ContractError(ErrorCode::invalid_input,path,"integer overflow");
        } else if constexpr(std::is_unsigned_v<T>) {
            if(value.get<std::int64_t>()<0) throw ContractError(ErrorCode::invalid_input,path,"negative unsigned value");
        }
        return value.get<T>();
    } else return value.get<T>();
}
template<class T> Json schema_for(const T& value) {
    Json schema;
    if constexpr(optional_type<T>::value) {schema=schema_for(value.value_or(typename T::value_type{}));schema["type"]=Json::array({schema["type"],"null"});}
    else if constexpr(Reflected<T>) {
        schema={{"type","object"},{"additionalProperties",false},{"properties",Json::object()},{"required",Json::array()}};
        visit_fields(value,[&](auto name,const auto& v){schema["properties"][name]=schema_for(v);schema["required"].push_back(name);});
    } else if constexpr(Sequence<T>) {
        schema={{"type","array"},{"items",schema_for(typename T::value_type{})}};
        if constexpr(!requires(T t){t.resize(1);}) {schema["minItems"]=value.size();schema["maxItems"]=value.size();}
    } else if constexpr(std::is_same_v<T,Json>) schema=Json::object();
    else if constexpr(std::is_same_v<T,bool>) schema={{"type","boolean"}};
    else if constexpr(std::is_integral_v<T>) schema={{"type","integer"},{"minimum",std::numeric_limits<T>::min()},{"maximum",std::numeric_limits<T>::max()}};
    else if constexpr(std::is_floating_point_v<T>) schema={{"type","number"}};
    else schema={{"type","string"}};
    schema["default"]=encode(value);
    return schema;
}
void validate_schema(const Json&,const Json&,const std::string&,ErrorCode=ErrorCode::invalid_option);
Json merge_defaults(const Json&,const Json&);
}
