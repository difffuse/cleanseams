// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include <set>
#include <tuple>
namespace cleanseams {
using namespace quality_detail;
namespace {
struct Segment {MeshIndex corner,next;Vec2 a,b;};
struct Bounds {std::int64_t x0,y0,x1,y1;std::uint64_t nx,ny;};
struct Entry {std::int64_t x,y;MeshIndex segment;};
using Pair=std::array<MeshIndex,2>;
struct BoundaryStopped {};
template<class T,class Less>
bool checked_sort(std::vector<T>& values,Less less,const std::function<bool()>& stopped) {
    try {
        std::sort(values.begin(),values.end(),[&](const auto& a,const auto& b) {
            if(stopped())throw BoundaryStopped{};
            return less(a,b);
        });
    } catch(const BoundaryStopped&) {return false;}
    return !stopped();
}
void memory_require(bool value,const char* field) {
    if(!value) throw ContractError(ErrorCode::memory_limit,field,"crossing storage exceeds memory allowance");
}
bool proper(Vec2 p,Vec2 q,Vec2 r,Vec2 s,double eps) {
    return orient(r,s,p)*orient(r,s,q)<-eps && orient(p,q,r)*orient(p,q,s)<-eps;
}
bool nonconforming_contact(Vec2 a,Vec2 b,Vec2 c,Vec2 d) {
    for(std::size_t k=0;k<2;++k)
        if(std::max(a[k],b[k])<std::min(c[k],d[k]) || std::max(c[k],d[k])<std::min(a[k],b[k]))return false;
    const int ac=orientation_sign(a,b,c),ad=orientation_sign(a,b,d);
    const int ca=orientation_sign(c,d,a),cb=orientation_sign(c,d,b);
    if(ac*ad<0 && ca*cb<0)return true;
    auto interior=[](Vec2 p,Vec2 q,Vec2 r,int orientation) {
        return orientation==0 && r!=p && r!=q &&
            r[0]>=std::min(p[0],q[0]) && r[0]<=std::max(p[0],q[0]) &&
            r[1]>=std::min(p[1],q[1]) && r[1]<=std::max(p[1],q[1]);
    };
    return interior(a,b,c,ac) || interior(a,b,d,ad) || interior(c,d,a,ca) || interior(c,d,b,cb);
}
// A chart has its own bucket keys; processing charts separately prevents
// inter-chart pairing and bounds peak temporary storage by one boundary.
bool candidates(const std::vector<Segment>& sides,double cell,std::uint64_t workspace,
                std::vector<Pair>& pairs,std::uint64_t& emitted_entries,const std::function<bool()>& stopped) {
    if(workspace==0 || !std::isfinite(cell)) return false;
    const auto entry_limit=std::min<std::uint64_t>(std::max<std::uint64_t>(20*sides.size(),100000),workspace/(2*sizeof(Entry)));
    if(entry_limit<sides.size()) return false;
    std::vector<Bounds> bounds(sides.size());std::uint64_t total=0;
    for(;;) {
        bool fit=true;total=0;
        for(std::size_t i=0;i<sides.size();++i) {
            if(stopped())return false;
            const auto& s=sides[i];auto& b=bounds[i];
            auto coordinate=[&](double value,std::int64_t& out) {
                const double x=std::floor(value/cell);
                if(!std::isfinite(x) || x < -0x1p63 || x >= 0x1p63) return false;
                out=static_cast<std::int64_t>(x);return true;
            };
            if(!coordinate(std::min(s.a[0],s.b[0]),b.x0) || !coordinate(std::max(s.a[0],s.b[0]),b.x1) ||
               !coordinate(std::min(s.a[1],s.b[1]),b.y0) || !coordinate(std::max(s.a[1],s.b[1]),b.y1)) {fit=false;break;}
            // Unsigned subtraction is exact here even across zero; check before
            // adding one so the full signed-coordinate range cannot wrap.
            const auto dx=static_cast<std::uint64_t>(b.x1)-static_cast<std::uint64_t>(b.x0);
            const auto dy=static_cast<std::uint64_t>(b.y1)-static_cast<std::uint64_t>(b.y0);
            if(dx>=entry_limit || dy>=entry_limit) {fit=false;break;}
            b.nx=dx+1;b.ny=dy+1;
            if(b.nx>(entry_limit-total)/b.ny) {fit=false;break;}
            total+=b.nx*b.ny;
        }
        if(fit) break;
        if(cell>std::numeric_limits<double>::max()/2) return false;
        cell*=2;
    }
    std::vector<Entry> entries;entries.reserve(static_cast<std::size_t>(total));
    for(std::size_t i=0;i<bounds.size();++i) {
        const auto& b=bounds[i];
        for(std::uint64_t y=0;y<b.ny;++y) for(std::uint64_t x=0;x<b.nx;++x) {
            if(stopped())return false;
            entries.push_back({b.x0+static_cast<std::int64_t>(x),b.y0+static_cast<std::int64_t>(y),static_cast<MeshIndex>(i)});
        }
    }
    emitted_entries+=total;
    if(!checked_sort(entries,[](const auto& a,const auto& b){return std::tie(a.x,a.y,a.segment)<std::tie(b.x,b.y,b.segment);},stopped))return false;
    const auto capacity=(workspace-total*sizeof(Entry))/sizeof(Pair);
    std::uint64_t pair_count=0;
    for(std::size_t start=0;start<entries.size();) {
        if(stopped())return false;
        auto end=start+1;while(end<entries.size() && entries[end].x==entries[start].x && entries[end].y==entries[start].y) {if(stopped())return false;++end;}
        const auto n=static_cast<std::uint64_t>(end-start);
        if(n && n-1>2*(capacity-pair_count)/n) return false;
        pair_count+=n*(n-1)/2;start=end;
    }
    pairs.reserve(static_cast<std::size_t>(pair_count));
    for(std::size_t start=0;start<entries.size();) {
        if(stopped())return false;
        auto end=start+1;while(end<entries.size() && entries[end].x==entries[start].x && entries[end].y==entries[start].y) {if(stopped())return false;++end;}
        for(auto a=start;a<end;++a) for(auto b=a+1;b<end;++b) {
            if(stopped())return false;
            pairs.push_back({std::min(entries[a].segment,entries[b].segment),std::max(entries[a].segment,entries[b].segment)});
        }
        start=end;
    }
    if(!checked_sort(pairs,std::less<Pair>{},stopped))return false;
    try {pairs.erase(std::unique(pairs.begin(),pairs.end(),[&](const auto& a,const auto& b){if(stopped())throw BoundaryStopped{};return a==b;}),pairs.end());}
    catch(const BoundaryStopped&) {return false;}
    return true;
}
}
bool segments_have_nonconforming_contact(Vec2 a,Vec2 b,Vec2 c,Vec2 d) {
    require(finite(a) && finite(b) && finite(c) && finite(d),"quality.contacts","contacts need finite points");
    return nonconforming_contact(a,b,c,d);
}
bool segments_cross(Vec2 p,Vec2 q,Vec2 r,Vec2 s,double eps) {
    require(std::isfinite(eps) && eps>=0,"quality.crossing_epsilon","epsilon must be finite and nonnegative");
    return proper(p,q,r,s,eps);
}
BoundaryCrossingReport boundary_crossings(const CutSurface& cut,CornerUvView uv,
    const BoundaryCrossingOptions& options,std::optional<std::span<const MeshIndex>> charts) {
    const auto& mesh=cut.mesh();const auto& ct=cut.tables();validate_view(mesh,uv);
    require(!options.cell || (std::isfinite(*options.cell) && *options.cell>0),"quality.crossing_cell","cell must be finite and positive");
    const auto budget=std::min<std::uint64_t>(options.memory_budget,std::numeric_limits<std::size_t>::max());
    const auto nf=mesh.input().face_offsets.size()-1;
    const auto selection_size=charts?charts->size():ct.n_charts;
    memory_require(selection_size<=budget/32 && options.ignored_faces.size()<=budget/32,"quality.crossings.inputs");
    for(auto f:options.ignored_faces) require(f>=0 && static_cast<std::uint64_t>(f)<nf,"quality.crossings.ignore","face out of range");
    if(charts) for(auto ch:*charts) require(ch<ct.n_charts,"quality.crossings.charts","chart out of range");
    BoundaryCrossingReport result;
    std::uint64_t polls=0;
    auto stopped=[&] {
        if(result.complete && options.stopped && (polls++ & 1023)==0 && options.stopped())result.complete=false;
        return !result.complete;
    };
    if(stopped())return result;
    if(charts) result.selected_charts.assign(charts->begin(),charts->end());
    else {result.selected_charts.resize(ct.n_charts);std::iota(result.selected_charts.begin(),result.selected_charts.end(),0);}
    if(!checked_sort(result.selected_charts,std::less<MeshIndex>{},stopped))return result;
    require(std::adjacent_find(result.selected_charts.begin(),result.selected_charts.end())==result.selected_charts.end(),"quality.crossings.charts","duplicate chart IDs");
    auto ignored=options.ignored_faces;if(!checked_sort(ignored,std::less<TransferIndex>{},stopped))return result;ignored.erase(std::unique(ignored.begin(),ignored.end()),ignored.end());
    std::uint64_t ns=0;
    for(auto ch:result.selected_charts) for(auto li:ct.chart_loops[ch]) {if(stopped())return result;ns+=ct.loops[li].size();}
    // Conservative live allocation allowance: segments, bounds, touched faces,
    // selected charts, count arrays, policy copies and per-chart traversal data.
    std::uint64_t base=0;
    for(auto [count,stride]:std::array<std::pair<std::uint64_t,std::uint64_t>,4>{{{ns,256},{ct.n_charts,32},{selection_size,8},{options.ignored_faces.size(),16}}}) {
        memory_require(count<=(budget-base)/stride,"quality.crossings.memory_budget");base+=count*stride;
    }
    const auto workspace=std::min(options.workspace_bytes,(budget-base)/2);
    const auto retained_limit=(budget-base-workspace)/sizeof(BoundaryCrossingPair);
    result.counts.assign(ct.n_charts,0);result.testable.assign(ct.n_charts,0);
    std::set<MeshIndex> touched_faces;
    auto append=[&](BoundaryCrossingPair pair) {
        if(!options.retain_pairs) return;
        if(result.pairs.size()==result.pairs.capacity()) {
            memory_require(result.pairs.size()<retained_limit,"quality.crossings.pairs");
            const auto next=std::min<std::uint64_t>(retained_limit,std::max<std::uint64_t>(8,2*result.pairs.capacity()));
            // Include old and new buffers while vector reserve reallocates.
            memory_require(next+result.pairs.capacity()<=retained_limit,"quality.crossings.pairs");
            result.pairs.reserve(static_cast<std::size_t>(next));
        }
        result.pairs.push_back(pair);
    };
    for(auto ch:result.selected_charts) {
        std::vector<Segment> sides;std::size_t raw_count=0;bool all_finite=true;
        for(auto li:ct.chart_loops[ch]) raw_count+=ct.loops[li].size();
        if(raw_count<4) continue;
        sides.reserve(raw_count);
        for(auto li:ct.chart_loops[ch]) for(auto c:ct.loops[li]) {
            if(stopped())return result;
            const auto face=mesh.tables().corner_face[c];
            if(std::binary_search(ignored.begin(),ignored.end(),static_cast<TransferIndex>(face))) continue;
            const auto next=mesh.corner_next(c);const auto a=point(uv,c),b=point(uv,next);
            sides.push_back({c,next,a,b});all_finite=all_finite && finite(a) && finite(b);
        }
        if(!all_finite || sides.size()<2) continue;
        result.testable[ch]=sides.size()>=4;
        auto consider=[&](std::size_t i,std::size_t j) {
            ++result.candidate_pairs;const auto& a=sides[i];const auto& b=sides[j];bool vertex_share=false;
            if(options.nonconforming_contacts) {
                if(nonconforming_contact(a.a,a.b,b.a,b.b)) {
                    ++result.counts[ch];append({ch,std::min(a.corner,b.corner),std::max(a.corner,b.corner)});
                    touched_faces.insert(mesh.tables().corner_face[a.corner]);touched_faces.insert(mesh.tables().corner_face[b.corner]);
                }
                return;
            }
            for(auto x:{a.corner,a.next}) for(auto y:{b.corner,b.next}) {
                if(ct.sector_of_corner[x]==ct.sector_of_corner[y]) return;
                vertex_share=vertex_share || mesh.input().corner_vertices[x]==mesh.input().corner_vertices[y];
            }
            vertex_share=vertex_share && ct.chart_faces[ch].size()==1;
            if(!proper(a.a,a.b,b.a,b.b,1e-12)) return;
            if(result.testable[ch] && !(options.single_face_slit && vertex_share)) {
                ++result.counts[ch];append({ch,std::min(a.corner,b.corner),std::max(a.corner,b.corner)});
            }
            // Source crossing_faces always enables the keyhole rule, and checks
            // the four-side minimum before removing ignored sides.
            if(!vertex_share) {touched_faces.insert(mesh.tables().corner_face[a.corner]);touched_faces.insert(mesh.tables().corner_face[b.corner]);}
        };
        double cell=options.cell.value_or(0);
        if(!options.cell) {
            for(const auto& s:sides) {if(stopped())return result;const auto x=s.b[0]-s.a[0],y=s.b[1]-s.a[1];cell+=std::sqrt(x*x+y*y);}
            cell=std::max(cell/static_cast<double>(sides.size())*2,1e-9);
        }
        std::vector<Pair> pairs;
        if(candidates(sides,cell,workspace,pairs,result.bucket_entries,stopped)) {
            for(auto p:pairs) {if(stopped())return result;consider(p[0],p[1]);}
        } else {
            if(stopped())return result;
            result.used_exhaustive=true;
            for(std::size_t i=0;i<sides.size();++i) for(std::size_t j=i+1;j<sides.size();++j) {if(stopped())return result;consider(i,j);}
        }
    }
    result.crossing_faces.assign(touched_faces.begin(),touched_faces.end());
    if(!checked_sort(result.pairs,std::less<BoundaryCrossingPair>{},stopped))return result;
    if(options.stopped && options.stopped())result.complete=false;
    return result;
}
}
