// SPDX-License-Identifier: GPL-3.0-or-later
// Case files (plan section 25.2): inline JSON or JSON plus checked binary arrays.
#include "cleanseams/io.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <set>
namespace cleanseams {
namespace {
using Json=nlohmann::json;
namespace fs=std::filesystem;
[[noreturn]] void fail(const std::string& field,const std::string& message,ErrorCode code=ErrorCode::invalid_input) {throw ContractError(code,field,message);}
// ----------------------------------------------------------------- SHA-256 (FIPS 180-4)
struct Sha256 {
    std::array<std::uint32_t,8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::array<std::uint8_t,64> block{};std::size_t used=0;std::uint64_t length=0;
    static constexpr std::array<std::uint32_t,64> k{
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
        0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
        0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
        0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    void compress() {
        std::array<std::uint32_t,64> w{};
        for(std::size_t i=0;i<16;++i)w[i]=std::uint32_t(block[4*i])<<24 | std::uint32_t(block[4*i+1])<<16 | std::uint32_t(block[4*i+2])<<8 | block[4*i+3];
        for(std::size_t i=16;i<64;++i) {
            const auto s0=std::rotr(w[i-15],7)^std::rotr(w[i-15],18)^(w[i-15]>>3),s1=std::rotr(w[i-2],17)^std::rotr(w[i-2],19)^(w[i-2]>>10);
            w[i]=w[i-16]+s0+w[i-7]+s1;
        }
        auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
        for(std::size_t i=0;i<64;++i) {
            const auto t1=hh+(std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25))+((e&f)^(~e&g))+k[i]+w[i];
            const auto t2=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
            hh=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=hh;
    }
    void update(std::span<const std::uint8_t> bytes) {
        for(auto byte:bytes) {block[used++]=byte;if(used==64) {compress();used=0;}}
        length+=bytes.size();
    }
    std::string hex() {
        const auto bits=length*8;std::uint8_t one=0x80;update({&one,1});
        std::uint8_t zero=0;while(used!=56)update({&zero,1});
        for(int i=7;i>=0;--i) {const auto byte=static_cast<std::uint8_t>(bits>>(8*i));update({&byte,1});}
        static constexpr char digits[]="0123456789abcdef";std::string out;
        for(auto word:h)for(int i=3;i>=0;--i) {const auto byte=(word>>(8*i))&0xff;out+=digits[byte>>4];out+=digits[byte&15];}
        return out;
    }
};
// ----------------------------------------------------------------- binary arrays
struct Dtype {const char* name;std::size_t itemsize;const char* ext;};
constexpr std::array<Dtype,4> dtypes{{{"float64-le",8,"f64le"},{"int64-le",8,"i64le"},{"uint64-le",8,"u64le"},{"uint8",1,"u8"}}};
const Dtype& dtype_of(const std::string& name,const std::string& field) {
    for(const auto& d:dtypes)if(name==d.name)return d;
    fail(field+".dtype","unsupported dtype");
}
struct Budget {
    std::uint64_t left;
    void take(std::uint64_t bytes,const std::string& field) {
        if(bytes>left)fail(field,"exceeds the memory budget",ErrorCode::memory_limit);
        left-=bytes;
    }
};
// A non-negative JSON integer (parsed files give unsigned, constructed values signed).
bool counting(const Json& j) {return j.is_number_unsigned() || (j.is_number_integer() && j.get<std::int64_t>()>=0);}
// A checked descriptor: returns the raw little-endian bytes and the shape.
std::pair<std::vector<std::uint8_t>,std::vector<std::uint64_t>> read_descriptor(const Json& d,const fs::path& base,const std::string& field,
                                                                                 const char* expected,Budget& budget) {
    static const std::set<std::string> keys{"file","dtype","shape","byte_offset","byte_length","sha256"};
    for(auto& [key,value]:d.items())if(!keys.contains(key))fail(field+"."+key,"unknown descriptor field");
    for(const auto& key:keys)if(!d.contains(key))fail(field+"."+key,"missing descriptor field");
    if(!d.at("file").is_string() || !d.at("dtype").is_string() || !d.at("sha256").is_string() || !d.at("shape").is_array() ||
       !counting(d.at("byte_offset")) || !counting(d.at("byte_length")))fail(field,"descriptor field has the wrong type");
    const auto file=d.at("file").get<std::string>();
    const auto& type=dtype_of(d.at("dtype").get<std::string>(),field);
    if(std::string(type.name)!=expected)fail(field+".dtype","this field requires "+std::string(expected));
    const fs::path relative(file);
    if(file.empty() || relative.is_absolute() || file.find('\\')!=std::string::npos || file.find(':')!=std::string::npos)fail(field+".file","array path must be relative");
    // NUL or other control bytes would be cut or reinterpreted by the operating system.
    if(std::any_of(file.begin(),file.end(),[](char c) {return static_cast<unsigned char>(c)<0x20;}))fail(field+".file","array path must not contain control characters");
    for(const auto& part:relative)if(part==".." || part==".")fail(field+".file","array path must not contain . or ..");
    std::uint64_t count=1;
    for(const auto& n:d.at("shape")) {
        if(!counting(n))fail(field+".shape","shape entries must be non-negative integers");
        const auto v=n.get<std::uint64_t>();
        if(v && count>std::numeric_limits<std::uint64_t>::max()/v)fail(field+".shape","shape overflows",ErrorCode::memory_limit);
        count*=v;
    }
    if(count>std::numeric_limits<std::uint64_t>::max()/type.itemsize)fail(field+".shape","shape overflows",ErrorCode::memory_limit);
    const auto length=d.at("byte_length").get<std::uint64_t>(),offset=d.at("byte_offset").get<std::uint64_t>();
    if(count*type.itemsize!=length)fail(field+".byte_length","byte length does not match shape and dtype");
    if(offset>std::numeric_limits<std::uint64_t>::max()-length)fail(field+".byte_offset","offset overflows");
    std::error_code ec;
    // Resolution failures (overlong names, permissions) are input errors, not exceptions.
    const auto root=fs::weakly_canonical(base,ec);if(ec)fail(field+".file","case directory cannot be resolved",ErrorCode::io_error);
    const auto target=fs::weakly_canonical(base/relative,ec);if(ec)fail(field+".file","array path cannot be resolved",ErrorCode::io_error);
    const auto rel=target.lexically_relative(root);
    if(rel.empty() || *rel.begin()=="..")fail(field+".file","array path escapes the case directory");
    if(!fs::is_regular_file(target,ec))fail(field+".file","array file does not exist",ErrorCode::io_error);
    const auto size=fs::file_size(target,ec);
    if(ec || offset+length>size)fail(field+".file","array file is truncated");
    budget.take(length,field);
    std::vector<std::uint8_t> bytes(length);
    std::ifstream stream(target,std::ios::binary);stream.seekg(static_cast<std::streamoff>(offset));
    if(length && !stream.read(reinterpret_cast<char*>(bytes.data()),static_cast<std::streamsize>(length)))fail(field+".file","cannot read the array",ErrorCode::io_error);
    Sha256 hash;hash.update(bytes);
    if(hash.hex()!=d.at("sha256").get<std::string>())fail(field+".sha256","array hash mismatch");
    std::vector<std::uint64_t> shape;for(const auto& n:d.at("shape"))shape.push_back(n.get<std::uint64_t>());
    return {std::move(bytes),std::move(shape)};
}
template<class T> T load_le(const std::uint8_t* p) {
    std::array<std::uint8_t,sizeof(T)> raw{};std::memcpy(raw.data(),p,sizeof(T));
    if constexpr(std::endian::native==std::endian::big)std::reverse(raw.begin(),raw.end());
    T value;std::memcpy(&value,raw.data(),sizeof(T));return value;
}
bool is_descriptor(const Json& j) {return j.is_object() && j.contains("file");}
// Flat numbers of an inline array (nested lists) or a descriptor, with the row width.
std::vector<double> floats(const Json& j,const fs::path& base,const std::string& field,std::size_t width,Budget& budget,bool finite=true) {
    std::vector<double> out;
    if(is_descriptor(j)) {
        auto [bytes,shape]=read_descriptor(j,base,field,"float64-le",budget);
        if(width>1 ? (shape.size()!=2 || shape[1]!=width) : shape.size()!=1)fail(field+".shape","unexpected array shape");
        for(std::size_t i=0;i<bytes.size();i+=8)out.push_back(load_le<double>(&bytes[i]));
    } else {
        if(!j.is_array())fail(field,"must be an array");
        budget.take(8*j.size()*std::max<std::size_t>(width,1),field);
        for(const auto& row:j) {
            if(width>1) {
                if(!row.is_array() || row.size()!=width)fail(field,"rows must have "+std::to_string(width)+" numbers");
                for(const auto& x:row) {if(!x.is_number())fail(field,"entries must be numbers");out.push_back(x.get<double>());}
            } else {if(!row.is_number())fail(field,"entries must be numbers");out.push_back(row.get<double>());}
        }
    }
    if(finite)for(double x:out)if(!std::isfinite(x))fail(field,"values must be finite");
    return out;
}
std::vector<std::int64_t> integers(const Json& j,const fs::path& base,const std::string& field,std::size_t width,const char* dtype,Budget& budget) {
    std::vector<std::int64_t> out;
    if(is_descriptor(j)) {
        auto [bytes,shape]=read_descriptor(j,base,field,dtype,budget);
        if(width>1 ? (shape.size()!=2 || shape[1]!=width) : shape.size()!=1)fail(field+".shape","unexpected array shape");
        const std::string type=dtype;
        if(type=="uint8")for(auto b:bytes)out.push_back(b);
        else for(std::size_t i=0;i<bytes.size();i+=8) {
            if(type=="uint64-le") {
                const auto v=load_le<std::uint64_t>(&bytes[i]);
                if(v>static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))fail(field,"value exceeds the index range");
                out.push_back(static_cast<std::int64_t>(v));
            } else out.push_back(load_le<std::int64_t>(&bytes[i]));
        }
    } else {
        if(!j.is_array())fail(field,"must be an array");
        budget.take(8*j.size()*std::max<std::size_t>(width,1),field);
        for(const auto& row:j) {
            auto one=[&](const Json& x) {if(!x.is_number_integer())fail(field,"entries must be integers");out.push_back(x.get<std::int64_t>());};
            if(width>1) {if(!row.is_array() || row.size()!=width)fail(field,"rows must have "+std::to_string(width)+" integers");for(const auto& x:row)one(x);}
            else one(row);
        }
    }
    return out;
}
Json descriptor_for(std::span<const std::uint8_t> bytes,const char* dtype,std::vector<std::uint64_t> shape,const fs::path& directory,const std::string& name) {
    const auto& type=dtype_of(dtype,name);
    const auto relative=fs::path("arrays")/(name+"."+type.ext);
    fs::create_directories(directory/"arrays");
    const auto target=directory/relative;
    std::error_code ec;bool same=false;
    if(fs::is_regular_file(target,ec) && fs::file_size(target,ec)==bytes.size()) {
        std::vector<char> existing(bytes.size());std::ifstream in(target,std::ios::binary);
        same=!bytes.size() || (in.read(existing.data(),static_cast<std::streamsize>(existing.size())) &&
             std::memcmp(existing.data(),bytes.data(),bytes.size())==0);
    }
    if(!same) {
        std::ofstream out(target,std::ios::binary|std::ios::trunc);
        if(!out.write(reinterpret_cast<const char*>(bytes.data()),static_cast<std::streamsize>(bytes.size())))fail(name,"cannot write the array",ErrorCode::io_error);
    }
    Sha256 hash;hash.update(bytes);
    return {{"file",relative.generic_string()},{"dtype",dtype},{"shape",shape},{"byte_offset",0},{"byte_length",bytes.size()},{"sha256",hash.hex()}};
}
template<class T> void append_le(std::vector<std::uint8_t>& out,T value) {
    std::array<std::uint8_t,sizeof(T)> raw{};std::memcpy(raw.data(),&value,sizeof(T));
    if constexpr(std::endian::native==std::endian::big)std::reverse(raw.begin(),raw.end());
    out.insert(out.end(),raw.begin(),raw.end());
}
}
std::string sha256_hex(std::span<const std::uint8_t> bytes) {Sha256 h;h.update(bytes);return h.hex();}
Json parse_untrusted_json(std::string_view text,const std::string& field,std::size_t max_depth) {
    // Structural nesting outside string literals; the parser itself is iterative.
    std::size_t depth=0;bool in_string=false,escaped=false;
    for(const char c:text) {
        if(in_string) {
            if(escaped)escaped=false;else if(c=='\\')escaped=true;else if(c=='"')in_string=false;
            continue;
        }
        if(c=='"')in_string=true;
        else if(c=='[' || c=='{') {if(++depth>max_depth)fail(field,"JSON nesting exceeds "+std::to_string(max_depth)+" levels");}
        else if((c==']' || c=='}') && depth>0)--depth;
    }
    try {return Json::parse(text);}
    catch(const Json::exception& error) {fail(field,std::string("invalid JSON: ")+error.what());}
}
std::string read_text_file(const fs::path& path,std::uint64_t limit,const char* field) {
    std::error_code ec;
    if(!fs::is_regular_file(path,ec))fail(field,"file does not exist: "+path.string(),ErrorCode::io_error);
    const auto size=fs::file_size(path,ec);
    if(ec)fail(field,"cannot read the file size",ErrorCode::io_error);
    if(size>limit)fail(field,"file exceeds the memory budget",ErrorCode::memory_limit);
    std::string text(size,'\0');std::ifstream stream(path,std::ios::binary);
    if(size && !stream.read(text.data(),static_cast<std::streamsize>(size)))fail(field,"cannot read the file",ErrorCode::io_error);
    return text;
}
std::vector<std::vector<MeshIndex>> case_faces(const MeshInput& mesh) {
    std::vector<std::vector<MeshIndex>> faces;
    for(std::size_t f=0;f+1<mesh.face_offsets.size();++f) {
        std::vector<MeshIndex> face;
        for(auto c=mesh.face_offsets[f];c<mesh.face_offsets[f+1];++c)face.push_back(static_cast<MeshIndex>(mesh.corner_vertices[c]));
        faces.push_back(std::move(face));
    }
    return faces;
}
CaseData case_from_json(const Json& j,const fs::path& base,const CaseLimits& limits) {
    if(!j.is_object())fail("case","a case must be a JSON object");
    Budget budget{limits.memory_budget};
    if(j.contains("schema") && j.at("schema")!="cleanseams.case")fail("schema","not a cleanseams.case document",ErrorCode::unsupported_version);
    CaseData c;
    if(j.contains("version")) {
        if(!j.at("version").is_number_integer() || j.at("version").get<std::int64_t>()!=1)fail("version","unsupported case version",ErrorCode::unsupported_version);
    }
    if(j.contains("coordinate_space")) {
        if(j.at("coordinate_space")!="object")fail("coordinate_space","only object space is supported");
    }
    for(const auto* key:{"positions","face_offsets","corner_vertices"})if(!j.contains(key))fail(key,"missing required field");
    const auto p=floats(j.at("positions"),base,"positions",3,budget);
    for(std::size_t i=0;i<p.size();i+=3)c.mesh.positions.push_back({p[i],p[i+1],p[i+2]});
    const auto offsets=integers(j.at("face_offsets"),base,"face_offsets",1,"uint64-le",budget);
    c.mesh.face_offsets.clear();
    for(auto o:offsets) {if(o<0)fail("face_offsets","offsets must be non-negative");c.mesh.face_offsets.push_back(static_cast<Offset>(o));}
    c.mesh.corner_vertices=integers(j.at("corner_vertices"),base,"corner_vertices",1,"int64-le",budget);
    const auto nv=c.mesh.positions.size(),nc=c.mesh.corner_vertices.size();
    const auto& fo=c.mesh.face_offsets;
    if(fo.empty() || fo.front()!=0 || fo.back()!=nc || !std::is_sorted(fo.begin(),fo.end()))fail("face_offsets","offsets must span corners monotonically from zero");
    for(auto v:c.mesh.corner_vertices)if(v<0 || static_cast<std::uint64_t>(v)>=nv)fail("corner_vertices","index out of range");
    auto pairs=[&](const char* key,std::vector<EdgeEndpoints>& out) {
        if(!j.contains(key) || j.at(key).is_null())return;
        const auto flat=integers(j.at(key),base,key,2,"int64-le",budget);
        for(std::size_t i=0;i<flat.size();i+=2) {
            for(std::size_t k=0;k<2;++k)if(flat[i+k]<0 || static_cast<std::uint64_t>(flat[i+k])>=nv)fail(key,"edge endpoint out of range");
            out.push_back({flat[i],flat[i+1]});
        }
    };
    pairs("seam_edges",c.seam_edges);pairs("locked_edges",c.locked_edges);
    if(j.contains("anchors") && !j.at("anchors").is_null()) {
        c.anchors=integers(j.at("anchors"),base,"anchors",1,"int64-le",budget);
        for(auto a:c.anchors)if(a<0 || static_cast<std::uint64_t>(a)>=nv)fail("anchors","index out of range");
    }
    if(j.contains("pins") && !j.at("pins").is_null()) {
        if(!j.at("pins").is_array())fail("pins","must be an array of {corner, uv} records");
        std::set<TransferIndex> seen;
        for(const auto& pin:j.at("pins")) {
            if(!pin.is_object() || !pin.contains("corner") || !pin.contains("uv") || !pin.at("corner").is_number_integer() || !pin.at("uv").is_array() ||
               pin.at("uv").size()!=2 || !pin.at("uv")[0].is_number() || !pin.at("uv")[1].is_number())fail("pins","pins are {corner, uv: [u, v]} records");
            CornerPin cp;cp.corner=pin.at("corner").get<TransferIndex>();cp.uv={pin.at("uv")[0].get<double>(),pin.at("uv")[1].get<double>()};
            if(cp.corner<0 || static_cast<std::uint64_t>(cp.corner)>=nc)fail("pins","pinned corner out of range");
            if(!std::isfinite(cp.uv[0]) || !std::isfinite(cp.uv[1]))fail("pins","pin coordinates must be finite");
            if(!seen.insert(cp.corner).second)fail("pins","duplicate corner pin");
            c.pins.push_back(cp);
        }
    }
    if(j.contains("uv") && !j.at("uv").is_null()) {
        // Non-finite UVs (binary arrays only; JSON has no NaN) mean no baseline (Gate A F6).
        const auto uv=floats(j.at("uv"),base,"uv",2,budget,false);
        if(uv.size()!=2*nc)fail("uv","UV count must equal corner count");
        std::vector<Vec2> rows;for(std::size_t i=0;i<uv.size();i+=2)rows.push_back({uv[i],uv[i+1]});
        c.uv=std::move(rows);
    }
    if(j.contains("operation")) {
        if(!j.at("operation").is_string())fail("operation","operation must be a string");
        c.operation=j.at("operation").get<std::string>();
    }
    if(c.operation!="pipeline" && c.operation!="postprocess_layout")fail("operation","unknown operation");
    if(c.operation=="postprocess_layout" && !c.uv)fail("uv","postprocess_layout requires supplied UVs");
    if(j.contains("selected_faces")) {
        const auto& selected=j.at("selected_faces");
        if(!selected.is_array() || (!selected.empty() && selected.size()!=c.mesh.face_offsets.size()-1))fail("selected_faces","mask must match original faces");
        for(const auto& value:selected) {
            if(!value.is_number_integer() || (value!=0 && value!=1))fail("selected_faces","mask must contain only zero or one");
            c.selected_faces.push_back(value.get<std::uint8_t>());
        }
    }
    static const std::set<std::string> known{"operation","selected_faces","schema","version","coordinate_space","positions","face_offsets","corner_vertices","seam_edges",
                                             "locked_edges","anchors","pins","uv"};
    for(auto& [key,value]:j.items())if(!known.contains(key))c.extra[key]=value;
    return c;
}
CaseData read_case(const fs::path& path,const CaseLimits& limits) {
    const auto text=read_text_file(path,limits.memory_budget,"case");
    const auto j=parse_untrusted_json(text,"case");
    CaseLimits rest=limits;rest.memory_budget-=std::min<std::uint64_t>(rest.memory_budget,text.size());
    return case_from_json(j,path.parent_path().empty()?fs::path("."):path.parent_path(),rest);
}
Json case_json(const CaseData& c,const fs::path& path,std::size_t inline_limit) {
    const auto directory=path.parent_path().empty()?fs::path("."):path.parent_path();
    const auto stem=path.stem().string();
    Json j=c.extra;
    if(c.operation!="pipeline")j["operation"]=c.operation;
    if(!c.selected_faces.empty())j["selected_faces"]=c.selected_faces;
    j["schema"]="cleanseams.case";j["version"]=c.version;j["coordinate_space"]=c.coordinate_space;
    auto emit=[&](const char* key,std::size_t scalars,auto inline_value,auto bytes,const char* dtype,std::vector<std::uint64_t> shape) {
        if(scalars>=inline_limit) {auto raw=bytes();j[key]=descriptor_for(raw,dtype,std::move(shape),directory,stem+"-"+key);}
        else j[key]=inline_value();
    };
    const auto nv=c.mesh.positions.size(),nc=c.mesh.corner_vertices.size();
    emit("positions",3*nv,[&] {return Json(c.mesh.positions);},[&] {std::vector<std::uint8_t> b;for(const auto& p:c.mesh.positions)for(double x:p)append_le(b,x);return b;},
         "float64-le",{nv,3});
    emit("face_offsets",c.mesh.face_offsets.size(),[&] {return Json(c.mesh.face_offsets);},
         [&] {std::vector<std::uint8_t> b;for(auto o:c.mesh.face_offsets)append_le(b,static_cast<std::uint64_t>(o));return b;},"uint64-le",{c.mesh.face_offsets.size()});
    emit("corner_vertices",nc,[&] {return Json(c.mesh.corner_vertices);},
         [&] {std::vector<std::uint8_t> b;for(auto v:c.mesh.corner_vertices)append_le(b,static_cast<std::int64_t>(v));return b;},"int64-le",{nc});
    for(const auto& [key,edges]:{std::pair{"seam_edges",&c.seam_edges},std::pair{"locked_edges",&c.locked_edges}})
        emit(key,2*edges->size(),[&] {return Json(*edges);},
             [&] {std::vector<std::uint8_t> b;for(const auto& e:*edges)for(auto v:e)append_le(b,static_cast<std::int64_t>(v));return b;},"int64-le",{edges->size(),2});
    j["anchors"]=c.anchors;
    Json pins=Json::array();for(const auto& p:c.pins)pins.push_back({{"corner",p.corner},{"uv",{p.uv[0],p.uv[1]}}});j["pins"]=pins;
    // JSON has no NaN or infinity: a map holding one is always written as a binary array.
    const bool nonfinite=c.uv && std::any_of(c.uv->begin(),c.uv->end(),[](const Vec2& p) {return !std::isfinite(p[0]) || !std::isfinite(p[1]);});
    if(!c.uv)j["uv"]=nullptr;
    else emit("uv",nonfinite?std::max<std::size_t>(inline_limit,2*c.uv->size()):2*c.uv->size(),[&] {return Json(*c.uv);},[&] {std::vector<std::uint8_t> b;for(const auto& t:*c.uv)for(double x:t)append_le(b,x);return b;},
              "float64-le",{c.uv->size(),2});
    return j;
}
void write_case(const fs::path& path,const CaseData& c,std::size_t inline_limit) {
    const auto j=case_json(c,path,inline_limit);
    std::ofstream stream(path,std::ios::binary|std::ios::trunc);
    if(!stream)fail("case","cannot open the output file",ErrorCode::io_error);
    stream<<j.dump(2)<<'\n';
    if(!stream)fail("case","cannot write the output file",ErrorCode::io_error);
}
}
