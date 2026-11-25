// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/detail/codec.hpp"
#include <cmath>
namespace cleanseams::detail {
namespace {
// GB performance: the validation below reads sub-schemas by reference and builds a field's
// path only when it is reported (a result holds hundreds of thousands of array elements).
// Checks, their order and every error's code, path and reason are the ones of the former
// copying implementation.
struct FieldPath {
    const FieldPath* parent=nullptr;
    const std::string* root=nullptr; // the caller's path (root segment only)
    const std::string* key=nullptr;  // object member
    std::size_t index=0;             // array element when key and root are null
    std::string str() const {
        if(root)return *root;
        const auto base=parent->str();
        if(key)return base.empty()?*key:base+"."+*key;
        return base+"["+std::to_string(index)+"]";
    }
};
const Json& empty_object() {static const Json empty=Json::object();return empty;}
const Json* member(const Json& schema,const std::string& key) {
    if(!schema.is_object())return nullptr;
    const auto found=schema.find(key);
    return found==schema.end()?nullptr:&*found;
}
void validate(const Json& value,const Json& schema,const FieldPath& path,ErrorCode code) {
    auto fail=[&](const char* reason){throw ContractError(code,path.str(),reason);};
    if(value.is_number_float() && !std::isfinite(value.get<double>())) fail("must be finite");
    if(!schema.is_object())(void)schema.value("type",Json("")); // non-object schema: the former value() error
    const Json* type_spec=member(schema,"type");
    if(value.is_null() && type_spec && type_spec->is_array()) return;
    static const std::string untyped;
    const std::string type=!type_spec?untyped:type_spec->is_array()?(*type_spec)[0].get<std::string>():type_spec->get<std::string>();
    if((type=="object" && !value.is_object()) || (type=="array" && !value.is_array()) ||
       (type=="string" && !value.is_string()) || (type=="boolean" && !value.is_boolean()) ||
       (type=="integer" && !value.is_number_integer()) || (type=="number" && !value.is_number())) fail("wrong type");
    if(const Json* e=member(schema,"enum");e && std::find(e->begin(),e->end(),value)==e->end()) {
        const auto p=path.str();
        if(p=="version" || p.ends_with(".version")) throw ContractError(ErrorCode::unsupported_version,p,"unsupported schema version");
        fail("unknown enum value");
    }
    if(value.is_number()) {
        if(const Json* low=member(schema,"minimum");low && value<*low) fail("below minimum");
        if(const Json* high=member(schema,"maximum");high && value>*high) fail("above maximum");
    }
    if(value.is_object()) {
        const Json* found=member(schema,"properties");const Json& props=found?*found:empty_object();
        for(auto it=value.begin();it!=value.end();++it) {
            const FieldPath child{&path,nullptr,&it.key(),0};
            const Json* sub=member(props,it.key());
            if(!sub && !schema.value("additionalProperties",true))throw ContractError(code,child.str(),"unknown field");
            validate(it.value(),sub?*sub:empty_object(),child,code);
        }
        if(const Json* required=member(schema,"required"))
            for(const auto& key:*required)
                if(!value.contains(key.get<std::string>())) {
                    const auto name=key.get<std::string>();
                    throw ContractError(code,FieldPath{&path,nullptr,&name,0}.str(),"missing required field");
                }
    }
    if(value.is_array()) {
        if(const Json* least=member(schema,"minItems");least && value.size()<least->get<std::size_t>()) fail("too few array elements");
        if(const Json* most=member(schema,"maxItems");most && value.size()>most->get<std::size_t>()) fail("too many array elements");
        const Json* found=member(schema,"items");const Json& items=found?*found:empty_object();
        for(std::size_t i=0;i<value.size();++i) validate(value[i],items,FieldPath{&path,nullptr,nullptr,i},code);
    }
}
}
void validate_schema(const Json& value,const Json& schema,const std::string& path,ErrorCode code) {
    validate(value,schema,FieldPath{nullptr,&path,nullptr,0},code);
}
Json merge_defaults(const Json& defaults,const Json& updates) {
    if(!defaults.is_object() || !updates.is_object()) return updates;
    auto result=defaults;
    for(auto it=updates.begin();it!=updates.end();++it)
        result[it.key()]=defaults.contains(it.key())?merge_defaults(defaults[it.key()],it.value()):it.value();
    return result;
}
}
