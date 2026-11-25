// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/io.py: read_obj and write_obj.
#include "cleanseams/io.hpp"
#include "cleanseams/errors.hpp"
#include "detail.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
namespace cleanseams {
namespace {
[[noreturn]] void bad(const std::string& message) {throw ContractError(ErrorCode::invalid_input,"obj",message);}
// Length of the Python str.isspace() character at text[i] (0 when it is not one): the
// ASCII separators \t\n\v\f\r, \x1c-\x1f and space, plus the UTF-8 encodings of U+0085,
// U+00A0, U+1680, U+2000-U+200A, U+2028, U+2029, U+202F, U+205F and U+3000. The source
// decodes the file as UTF-8 text and splits records with str.split().
std::size_t python_space(std::string_view text,std::size_t i) {
    const auto at=[&](std::size_t k) {return i+k<text.size()?static_cast<unsigned char>(text[i+k]):0u;};
    const auto c=at(0);
    if((c>=0x09 && c<=0x0d) || (c>=0x1c && c<=0x20))return 1;
    if(c==0xc2 && (at(1)==0x85 || at(1)==0xa0))return 2;
    if(c==0xe1 && at(1)==0x9a && at(2)==0x80)return 3;
    if(c==0xe2 && at(1)==0x80 && ((at(2)>=0x80 && at(2)<=0x8a) || at(2)==0xa8 || at(2)==0xa9 || at(2)==0xaf))return 3;
    if(c==0xe2 && at(1)==0x81 && at(2)==0x9f)return 3;
    if(c==0xe3 && at(1)==0x80 && at(2)==0x80)return 3;
    return 0;
}
// Python str.split() with no separator.
std::vector<std::string_view> split(std::string_view text) {
    std::vector<std::string_view> out;std::size_t i=0;
    while(i<text.size()) {
        while(i<text.size())if(const auto n=python_space(text,i))i+=n;else break;
        std::size_t j=i;while(j<text.size() && !python_space(text,j))++j;
        if(j>i)out.push_back(text.substr(i,j-i));
        i=j;
    }
    return out;
}
double python_float(std::string_view token) {
    const auto value=output_detail::python_float_text(token);
    if(!value)bad("could not convert string to float: '"+std::string(token)+"'");
    return *value;
}
// Python int(token) in base 10 for a face index: an optional sign, then ASCII digits with
// single underscores between them (PEP 515); no fractions, exponents or second sign.
// Returns the value, or nullopt with the canonical decimal text when it exceeds int64
// (Python keeps it exactly; such an index is outside any mesh).
struct PythonInt {std::optional<std::int64_t> value;std::string text;};
PythonInt python_int(std::string_view token) {
    std::string_view body=token;bool negative=false;
    if(!body.empty() && (body[0]=='+' || body[0]=='-')) {negative=body[0]=='-';body.remove_prefix(1);}
    std::string digits;bool ok=!body.empty();
    for(std::size_t i=0;ok && i<body.size();++i) {
        if(std::isdigit(static_cast<unsigned char>(body[i])))digits.push_back(body[i]);
        else ok=body[i]=='_' && i>0 && i+1<body.size() && std::isdigit(static_cast<unsigned char>(body[i-1])) &&
                std::isdigit(static_cast<unsigned char>(body[i+1]));
    }
    if(!ok)bad("invalid literal for int() with base 10: '"+std::string(token)+"'");
    const auto first=digits.find_first_not_of('0');
    const std::string canonical=first==std::string::npos?"0":(negative?"-":"")+digits.substr(first);
    std::uint64_t magnitude=0;const auto [ptr,ec]=std::from_chars(digits.data(),digits.data()+digits.size(),magnitude);
    const auto limit=static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if(ec!=std::errc() || ptr!=digits.data()+digits.size() || magnitude>limit)return {std::nullopt,canonical};
    const auto value=static_cast<std::int64_t>(magnitude);
    return {negative?-value:value,canonical};
}
MeshIndex resolve(std::string_view token,std::size_t count) {
    const auto parsed=python_int(token);
    auto outside=[&] {bad("OBJ index "+parsed.text+" is outside the "+std::to_string(count)+" defined elements");};
    if(!parsed.value)outside();
    const auto i=*parsed.value;
    if(i==0)bad("OBJ indices are one-based; zero is invalid");
    const auto result=i>0?i-1:static_cast<std::int64_t>(count)+i;
    if(!(result>=0 && result<static_cast<std::int64_t>(count)))outside();
    return static_cast<MeshIndex>(result);
}
std::string g17(double x) {
    if(std::isnan(x))return "nan";
    if(std::isinf(x))return x>0?"inf":"-inf";
    char buffer[64];std::snprintf(buffer,sizeof buffer,"%.17g",x);return buffer;
}
}
ObjMesh parse_obj(std::string_view text,std::uint64_t memory_budget) {
    if(text.size()>memory_budget)throw ContractError(ErrorCode::memory_limit,"obj","OBJ text exceeds the memory budget");
    ObjMesh out;std::vector<Vec2> vt;std::vector<std::optional<std::vector<MeshIndex>>> ft;std::uint64_t used=0;
    auto charge=[&](std::uint64_t bytes) {used+=bytes;if(used>memory_budget)throw ContractError(ErrorCode::memory_limit,"obj","OBJ exceeds the memory budget");};
    // The source reads the file in text mode: universal newlines end a line at "\n",
    // "\r\n" or a lone "\r".
    std::size_t start=0;
    while(start<=text.size()) {
        auto end=text.find_first_of("\r\n",start);if(end==std::string_view::npos)end=text.size();
        auto line=text.substr(start,end-start);
        start=end+(end+1<text.size() && text[end]=='\r' && text[end+1]=='\n'?2:1);
        if(const auto hash=line.find('#');hash!=std::string_view::npos)line=line.substr(0,hash);
        const auto tokens=split(line);
        if(tokens.empty()) {if(end==text.size())break;continue;}
        if(tokens[0]=="v") {
            if(tokens.size()<4)bad("vertex records need three coordinates");
            out.positions.push_back({python_float(tokens[1]),python_float(tokens[2]),python_float(tokens[3])});charge(24);
        } else if(tokens[0]=="vt") {
            if(tokens.size()<3)bad("texture records need two coordinates");
            vt.push_back({python_float(tokens[1]),python_float(tokens[2])});charge(16);
        } else if(tokens[0]=="f") {
            std::vector<MeshIndex> vs,ts;
            for(std::size_t k=1;k<tokens.size();++k) {
                const auto tok=tokens[k];const auto slash=tok.find('/');
                vs.push_back(resolve(tok.substr(0,slash),out.positions.size()));
                if(slash!=std::string_view::npos) {
                    auto rest=tok.substr(slash+1);const auto second=rest.find('/');const auto part=rest.substr(0,second);
                    if(!part.empty())ts.push_back(resolve(part,vt.size()));
                }
            }
            charge(8*(vs.size()+ts.size())+16);
            ft.push_back(ts.size()==vs.size()?std::optional(ts):std::nullopt);
            out.faces.push_back(std::move(vs));
        }
        if(end==text.size())break;
    }
    if(!vt.empty() && std::all_of(ft.begin(),ft.end(),[](const auto& t) {return t.has_value();})) {
        std::vector<Vec2> uv;for(const auto& t:ft)for(auto i:*t)uv.push_back(vt[i]);
        out.uv=std::move(uv);
    }
    return out;
}
ObjMesh read_obj(const std::filesystem::path& path,std::uint64_t memory_budget) {
    return parse_obj(read_text_file(path,memory_budget,"obj"),memory_budget);
}
std::string format_obj(std::span<const Vec3> positions,const std::vector<std::vector<MeshIndex>>& faces,std::optional<std::span<const Vec2>> uv,
                       std::span<const EdgeEndpoints> seam_lines) {
    std::size_t corners=0;for(const auto& f:faces)corners+=f.size();
    if(uv && uv->size()!=corners)throw ContractError(ErrorCode::invalid_input,"obj.uv","UV count must equal corner count");
    std::ostringstream out;
    for(const auto& p:positions)out<<"v "<<g17(p[0])<<' '<<g17(p[1])<<' '<<g17(p[2])<<'\n';
    if(uv)for(const auto& t:*uv)out<<"vt "<<g17(t[0])<<' '<<g17(t[1])<<'\n';
    std::size_t c=0;
    for(const auto& f:faces) {
        out<<'f';
        for(std::size_t i=0;i<f.size();++i) {
            out<<' '<<f[i]+1;
            if(uv)out<<'/'<<c+i+1;
        }
        out<<'\n';if(uv)c+=f.size();
    }
    for(const auto& e:seam_lines)out<<"l "<<e[0]+1<<' '<<e[1]+1<<'\n';
    return out.str();
}
void write_obj(const std::filesystem::path& path,std::span<const Vec3> positions,const std::vector<std::vector<MeshIndex>>& faces,
               std::optional<std::span<const Vec2>> uv,std::span<const EdgeEndpoints> seam_lines) {
    const auto text=format_obj(positions,faces,uv,seam_lines);
    std::ofstream stream(path,std::ios::binary);
    if(!stream)throw ContractError(ErrorCode::io_error,"obj","cannot open the output file");
    stream<<text;if(!stream)throw ContractError(ErrorCode::io_error,"obj","cannot write the output file");
}
}
