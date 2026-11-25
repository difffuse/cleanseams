// SPDX-License-Identifier: GPL-3.0-or-later
#include "detail.hpp"
#include <array>
#include <optional>
namespace cleanseams {
using namespace quality_detail;
nlohmann::json ChartValidity::to_json() const {
    auto valid_json=nlohmann::json::array(),measured_json=nlohmann::json::array();
    for(auto x:valid) valid_json.push_back(x!=0);for(auto x:measured) measured_json.push_back(x!=0);
    return {{"flipped",flipped},{"nonfinite",nonfinite},{"collapsed",collapsed},{"crossings",crossings},
        {"valid",valid_json},{"measured",measured_json}};
}
QualityEvaluation::QualityEvaluation(std::shared_ptr<const QualityGeometry> geometry,CornerUvView uv,
    const QualityPolicy& policy,std::optional<std::span<const MeshIndex>> charts):snapshot_(geometry,uv,policy,charts),whole_snapshot_(!charts) {
    const auto& cut=geometry->cut();const auto& ct=cut.tables();const auto n=ct.n_charts;
    BoundaryCrossingOptions metric_options;metric_options.memory_budget=geometry->memory_budget();metric_options.ignored_faces=snapshot_.policy().exempt_faces;
    // The source _metrics passes exempt=False straight to boundary_crossings.
    // Its scalar face-mask conversion interprets False as face ID zero. Preserve
    // this raw diagnostic edge case; chart_validity uses the actual empty mask.
    if(!policy.exempt_input_degenerate && cut.mesh().input().face_offsets.size()>1) metric_options.ignored_faces={0};
    auto canonicalize=[](auto& ids){std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());};
    canonicalize(metric_options.ignored_faces);
    metric_crossings_=std::make_shared<const BoundaryCrossingReport>(boundary_crossings(cut,uv,metric_options,charts));
    BoundaryCrossingOptions validity_options;validity_options.memory_budget=geometry->memory_budget();
    for(const auto& p:snapshot_.polygons()) if(p.exempt) validity_options.ignored_faces.push_back(p.face);
    canonicalize(validity_options.ignored_faces);
    if(metric_options.ignored_faces==validity_options.ignored_faces) validation_crossings_=metric_crossings_;
    else validation_crossings_=std::make_shared<const BoundaryCrossingReport>(boundary_crossings(cut,uv,validity_options,charts));
    auto& v=validity_;v.flipped.assign(n,0);v.nonfinite.assign(n,0);v.collapsed.assign(n,0);
    v.crossings=validation_crossings_->counts;v.valid.assign(n,1);v.measured.assign(n,0);
    std::vector<MeshIndex> offenders=validation_crossings_->crossing_faces;
    std::vector<MeshIndex> nonfinite_faces;
    for(const auto& p:snapshot_.polygons()) {
        v.measured[p.chart]=1;
        if(p.exempt) continue;
        v.flipped[p.chart]+=p.raw_flipped;v.nonfinite[p.chart]+=!p.finite;
        if(!p.finite) nonfinite_faces.push_back(p.face);
        if(p.raw_flipped || !p.finite) offenders.push_back(p.face);
    }
    for(const auto& t:snapshot_.triangles()) {
        const bool collapsed=t.area3>1e-14 && std::abs(t.raw_area2)<=1e-12*t.area3 && !t.exempt;
        if(collapsed) offenders.push_back(t.face);
        const bool finite_triangle=std::isfinite(t.raw_area2) && !std::binary_search(nonfinite_faces.begin(),nonfinite_faces.end(),t.face);
        v.collapsed[t.chart]+=collapsed && finite_triangle;
    }
    canonicalize(offenders);
    // Boundary-crossing faces were already filtered by the same exemption mask.
    offending_faces_=cut_detail::rows(n,[&](auto add){for(auto face:offenders) add(ct.chart_of_face[face],face);});
    for(MeshIndex ch=0;ch<n;++ch) v.valid[ch]=v.flipped[ch]==0 && v.nonfinite[ch]==0 && v.collapsed[ch]==0 && v.crossings[ch]==0;
}
bool QualityEvaluation::chart_valid(MeshIndex chart) const {
    require(chart<validity_.measured.size() && validity_.measured[chart],"quality.chart","chart was not measured");return validity_.valid[chart]!=0;
}
std::vector<MeshIndex> QualityEvaluation::offending_faces(MeshIndex chart) const {
    require(chart<validity_.measured.size() && validity_.measured[chart],"quality.chart","chart was not measured");
    const auto row=offending_faces_[chart];return {row.begin(),row.end()};
}
nlohmann::json QualityEvaluation::metrics() const {
    require(whole_snapshot_,"quality.metrics","whole-snapshot metrics require all charts");
    const auto r=snapshot_.report();auto result=r.to_json();
    const auto crossings=std::accumulate(metric_crossings_->counts.begin(),metric_crossings_->counts.end(),std::uint64_t{0});
    result["crossings"]=crossings;result["valid"]=r.nonfinite==0 && r.flipped==0 && crossings==0 && r.collapsed_uv==0;
    const auto& ct=snapshot_.geometry().cut().tables();result["charts"]=ct.n_charts;
    result["disks"]=std::count(ct.chart_is_disk.begin(),ct.chart_is_disk.end(),std::uint8_t{1});return result;
}
// ------------------------------------------------ strict overlap (Phase B, T28)
namespace {
// Knuth/Dekker error-free transformations; std::fma is explicit, never a contraction.
inline void two_sum(double a,double b,double& sum,double& error) {
    sum=a+b;const double bv=sum-a;const double av=sum-bv;error=(a-av)+(b-bv);
}
inline void two_product(double a,double b,double& product,double& error) {product=a*b;error=std::fma(a,b,-product);}
// Shewchuk's orient2d filter bound (3 + 16 eps) eps for the double evaluation below.
constexpr double orientation_bound=3.3306690738754716e-16;
int sign_of(double x) {return (x>0)-(x<0);}
// Exact sign of bx*cy - bx*ay - ax*cy - cx*by + cx*ay + ax*by: each product is split into
// two doubles, the twelve terms are summed as a nonoverlapping expansion (grow-expansion)
// and the most significant nonzero component carries the sign.
int exact_orientation(Vec2 a,Vec2 b,Vec2 c) {
    std::array<double,12> terms{};std::size_t k=0;
    auto product=[&](double x,double y,bool negative) {
        double p=0,e=0;two_product(x,y,p,e);terms[k++]=negative?-p:p;terms[k++]=negative?-e:e;
    };
    product(b[0],c[1],false);product(b[0],a[1],true);product(a[0],c[1],true);
    product(c[0],b[1],true);product(c[0],a[1],false);product(a[0],b[1],false);
    std::array<double,13> expansion{};std::size_t length=0;
    for(const double term:terms) {
        double q=term;
        for(std::size_t i=0;i<length;++i) {double sum=0,error=0;two_sum(q,expansion[i],sum,error);expansion[i]=error;q=sum;}
        expansion[length++]=q;
    }
    for(std::size_t i=length;i-->0;) if(expansion[i]!=0) return sign_of(expansion[i]);
    return 0;
}
double twice_area(const std::array<Vec2,3>& t) {return orient(t[0],t[1],t[2]);}
// The triangle with counter-clockwise vertices; nullopt when degenerate.
std::optional<std::array<Vec2,3>> counter_clockwise(const std::array<Vec2,3>& t,std::uint64_t* fallbacks) {
    bool exact=false;const int s=orientation_sign(t[0],t[1],t[2],&exact);
    if(exact && fallbacks) ++*fallbacks;
    if(s==0) return std::nullopt;
    if(s>0) return t;
    return std::array<Vec2,3>{t[0],t[2],t[1]};
}
// True when an edge line of counter-clockwise `a` has all of `b` in its closed outer half-plane.
bool separated_by_edges(const std::array<Vec2,3>& a,const std::array<Vec2,3>& b,std::uint64_t* fallbacks) {
    for(std::size_t i=0;i<3;++i) {
        const auto p=a[i],q=a[(i+1)%3];bool all_outside=true;
        for(const auto& r:b) {
            bool exact=false;const int s=orientation_sign(p,q,r,&exact);
            if(exact && fallbacks) ++*fallbacks;
            if(s>0) {all_outside=false;break;}
        }
        if(all_outside) return true;
    }
    return false;
}
struct OverlapTriangle {
    std::array<Vec2,3> p{};Vec2 low{},high{};double area=0;
    MeshIndex id=0,face=0,chart=0;
};
// m5: the broad phase is a bounding volume hierarchy over the chart-local UV boxes, traversed
// pairwise: every pair of boxes whose open interiors intersect is emitted exactly once, as
// (lower index, higher index), whatever the chart's shape: O(n log n + pairs). The boxes are
// ordered along a Z-order (Morton) curve of their centres, normalized per axis to the group's
// bounds (ties: the index), and the order is halved recursively down to leaves of at most 8
// boxes; node bounds are the union of their boxes. The uniform square grid it replaces
// coarsened on long strip charts until a column held a large share of the triangles
// (quadratic work: an 80k-triangle strip took seconds where a square took a few hundredths).
struct BoxTree {
    struct Node {Vec2 low{},high{};std::uint32_t begin=0,end=0,left=0,right=0;bool leaf() const {return left==0;}};
    std::vector<Node> nodes;std::vector<std::uint32_t> order;std::uint64_t leaves=0;
};
constexpr std::uint32_t box_leaf_size=8;
bool open_overlap(const UvBox& a,const UvBox& b) {return a.low[0]<b.high[0] && b.low[0]<a.high[0] && a.low[1]<b.high[1] && b.low[1]<a.high[1];}
// The 31 low bits of x on the even bit positions of the result.
std::uint64_t spread_bits(std::uint64_t v) {
    v&=0x7fffffffULL;
    v=(v|(v<<16))&0x0000ffff0000ffffULL;v=(v|(v<<8))&0x00ff00ff00ff00ffULL;v=(v|(v<<4))&0x0f0f0f0f0f0f0f0fULL;
    v=(v|(v<<2))&0x3333333333333333ULL;v=(v|(v<<1))&0x5555555555555555ULL;
    return v;
}
BoxTree build_box_tree(std::span<const UvBox> boxes,std::uint64_t memory_budget,std::uint64_t fixed_bytes) {
    BoxTree tree;const auto n=boxes.size();
    require(n<std::numeric_limits<std::uint32_t>::max()/4,"quality.overlap","too many triangles for the overlap broad phase");
    const auto max_nodes=2*(static_cast<std::uint64_t>(n)/(box_leaf_size/2)+1);
    const auto bytes=fixed_bytes+max_nodes*sizeof(BoxTree::Node)+n*(sizeof(std::uint32_t)+2*sizeof(std::uint64_t));
    if(bytes>memory_budget) throw ContractError(ErrorCode::memory_limit,"quality.overlap","overlap broad phase exceeds the memory budget");
    constexpr double inf=std::numeric_limits<double>::infinity();
    Vec2 clow{inf,inf},chigh{-inf,-inf};
    for(const auto& b:boxes)for(std::size_t d=0;d<2;++d) {const double c=0.5*(b.low[d]+b.high[d]);clow[d]=std::min(clow[d],c);chigh[d]=std::max(chigh[d],c);}
    std::vector<std::pair<std::uint64_t,std::uint32_t>> keyed(n);
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t key=0;
        for(std::size_t d=0;d<2;++d) {
            const double extent=chigh[d]-clow[d],c=0.5*(boxes[i].low[d]+boxes[i].high[d]);
            const double t=extent>0?(c-clow[d])/extent:0.0;
            const auto q=static_cast<std::uint64_t>(std::clamp(t,0.0,1.0)*2147483647.0);
            key|=spread_bits(q)<<d;
        }
        keyed[i]={key,static_cast<std::uint32_t>(i)};
    }
    std::sort(keyed.begin(),keyed.end());
    tree.order.resize(n);for(std::size_t i=0;i<n;++i)tree.order[i]=keyed[i].second;
    keyed.clear();keyed.shrink_to_fit();
    tree.nodes.reserve(static_cast<std::size_t>(max_nodes));tree.nodes.push_back({});tree.nodes[0].end=static_cast<std::uint32_t>(n);
    // Children are created after their parent: the reverse creation order fills leaves first.
    for(std::size_t index=0;index<tree.nodes.size();++index) {
        const auto begin=tree.nodes[index].begin,end=tree.nodes[index].end;
        if(end-begin<=box_leaf_size) {++tree.leaves;continue;}
        const auto mid=begin+(end-begin)/2;
        const auto left=static_cast<std::uint32_t>(tree.nodes.size());
        tree.nodes[index].left=left;tree.nodes[index].right=left+1;
        BoxTree::Node l,r;l.begin=begin;l.end=mid;r.begin=mid;r.end=end;
        tree.nodes.push_back(l);tree.nodes.push_back(r);
    }
    for(std::size_t index=tree.nodes.size();index-->0;) {
        auto& node=tree.nodes[index];
        if(node.leaf()) {
            node.low={inf,inf};node.high={-inf,-inf};
            for(auto k=node.begin;k<node.end;++k) {
                const auto& b=boxes[tree.order[k]];
                for(std::size_t d=0;d<2;++d) {node.low[d]=std::min(node.low[d],b.low[d]);node.high[d]=std::max(node.high[d],b.high[d]);}
            }
        } else {
            const auto& a=tree.nodes[node.left];const auto& b=tree.nodes[node.right];
            for(std::size_t d=0;d<2;++d) {node.low[d]=std::min(a.low[d],b.low[d]);node.high[d]=std::max(a.high[d],b.high[d]);}
        }
    }
    return tree;
}
// Every pair (i < j) of `boxes` with intersecting open interiors, each exactly once; returns
// the number of leaves.
template<class Visit>
std::uint64_t for_each_box_pair(std::span<const UvBox> boxes,std::uint64_t memory_budget,std::uint64_t fixed_bytes,Visit&& visit) {
    if(boxes.size()<2) return boxes.size();
    const auto tree=build_box_tree(boxes,memory_budget,fixed_bytes);
    const auto& nodes=tree.nodes;const auto& order=tree.order;
    auto emit=[&](std::uint32_t x,std::uint32_t y) {
        const auto i=std::min(x,y),j=std::max(x,y);
        if(open_overlap(boxes[i],boxes[j])) visit(static_cast<std::size_t>(i),static_cast<std::size_t>(j));
    };
    // A node paired with itself covers the pairs inside it; two distinct nodes the pairs across.
    std::vector<std::array<std::uint32_t,2>> stack{{0,0}};
    while(!stack.empty()) {
        const auto [a,b]=stack.back();stack.pop_back();
        const auto& A=nodes[a];const auto& B=nodes[b];
        if(a==b) {
            if(A.leaf()) {for(auto i=A.begin;i<A.end;++i)for(auto j=i+1;j<A.end;++j) emit(order[i],order[j]);}
            else {stack.push_back({A.left,A.right});stack.push_back({A.right,A.right});stack.push_back({A.left,A.left});}
            continue;
        }
        // Node boxes contain their boxes, so two boxes with intersecting open interiors lie in nodes
        // whose open interiors intersect too (a.low < b.high <= B.high and A.low <= a.low): the
        // open test prunes touching subtrees and loses no pair.
        if(!(A.low[0]<B.high[0] && B.low[0]<A.high[0] && A.low[1]<B.high[1] && B.low[1]<A.high[1])) continue;
        if(A.leaf() && B.leaf()) {for(auto i=A.begin;i<A.end;++i)for(auto j=B.begin;j<B.end;++j) emit(order[i],order[j]);}
        else if(B.leaf() || (!A.leaf() && A.end-A.begin>=B.end-B.begin)) {stack.push_back({A.right,b});stack.push_back({A.left,b});}
        else {stack.push_back({a,B.right});stack.push_back({a,B.left});}
    }
    return tree.leaves;
}
}
std::vector<std::array<std::size_t,2>> uv_box_pairs(std::span<const UvBox> boxes,std::uint64_t memory_budget) {
    for(const auto& b:boxes) require(finite(b.low) && finite(b.high),"quality.overlap","boxes need finite corners");
    std::vector<std::array<std::size_t,2>> out;
    for_each_box_pair(boxes,memory_budget,boxes.size()*sizeof(UvBox),[&](std::size_t i,std::size_t j) {out.push_back({i,j});});
    std::sort(out.begin(),out.end());
    return out;
}
namespace {
OverlapReport overlaps(const QualityGeometry& geometry,CornerUvView uv,const OverlapOptions& options,
                       std::optional<std::span<const MeshIndex>> charts,bool cross) {
    require(std::isfinite(options.area_tolerance) && options.area_tolerance>=0,"quality.overlap","area tolerance must be finite and nonnegative");
    const auto& cut=geometry.cut();const auto& mesh=cut.mesh();const auto& ct=cut.tables();validate_view(mesh,uv);
    const auto nf=mesh.input().face_offsets.size()-1;
    std::vector<std::uint8_t> ignored(nf,0);
    for(auto f:options.ignored_faces) {require(f>=0 && static_cast<std::uint64_t>(f)<nf,"quality.overlap","ignored face out of range");ignored[static_cast<std::size_t>(f)]=1;}
    if(options.ignore_input_degenerate)for(std::size_t f=0;f<nf;++f)ignored[f]=ignored[f] || mesh.geometry().face_degenerate_kind[f]!=0;
    OverlapReport report;report.cross_chart=cross;report.counts.assign(ct.n_charts,0);report.tested.assign(ct.n_charts,1);
    if(charts) {
        for(auto ch:*charts) require(ch<ct.n_charts,"quality.overlap","chart out of range");
        report.selected_charts.assign(charts->begin(),charts->end());
    } else {report.selected_charts.resize(ct.n_charts);std::iota(report.selected_charts.begin(),report.selected_charts.end(),MeshIndex{0});}
    std::sort(report.selected_charts.begin(),report.selected_charts.end());
    report.selected_charts.erase(std::unique(report.selected_charts.begin(),report.selected_charts.end()),report.selected_charts.end());
    const auto& tri=mesh.triangulation();
    auto collect=[&](MeshIndex ch,std::vector<OverlapTriangle>& out) {
        for(auto id:geometry.chart_triangles()[ch]) {
            const auto face=tri.faces[id];if(ignored[face]) continue;
            ++report.triangles;std::array<Vec2,3> p{};bool good=true;
            for(std::size_t k=0;k<3;++k) {p[k]=point(uv,tri.corners[id][k]);good=good && finite(p[k]);}
            if(!good) {++report.untestable;report.tested[ch]=0;continue;}
            const auto ccw=counter_clockwise(p,&report.exact_fallbacks);
            if(!ccw) {++report.degenerate;continue;}
            OverlapTriangle t;t.p=*ccw;t.id=id;t.face=face;t.chart=ch;t.area=0.5*std::abs(twice_area(t.p));
            t.low={std::min({t.p[0][0],t.p[1][0],t.p[2][0]}),std::min({t.p[0][1],t.p[1][1],t.p[2][1]})};
            t.high={std::max({t.p[0][0],t.p[1][0],t.p[2][0]}),std::max({t.p[0][1],t.p[1][1],t.p[2][1]})};
            out.push_back(t);
        }
    };
    std::vector<MeshIndex> faces;
    // Retained pairs: the max_pairs smallest by triangle IDs, independent of the visiting order.
    auto by_ids=[](const OverlapPair& x,const OverlapPair& y) {return std::pair(x.first_triangle,x.second_triangle)<std::pair(y.first_triangle,y.second_triangle);};
    auto retain=[&](const OverlapPair& pair) {
        report.pairs.push_back(pair);
        if(report.pairs.size()>=2*std::max<std::size_t>(options.max_pairs,1)) {
            std::nth_element(report.pairs.begin(),report.pairs.begin()+static_cast<std::ptrdiff_t>(options.max_pairs),report.pairs.end(),by_ids);
            report.pairs.resize(options.max_pairs);report.pairs_truncated=true;
        }
    };
    auto test_group=[&](const std::vector<OverlapTriangle>& tris) {
        if(tris.size()<2) return;
        std::vector<UvBox> boxes;boxes.reserve(tris.size());for(const auto& t:tris) boxes.push_back({t.low,t.high});
        // Interiors can meet only when the AABB interiors do: every such pair once, lower index first.
        report.cells+=for_each_box_pair(boxes,options.memory_budget,tris.size()*(sizeof(OverlapTriangle)+sizeof(UvBox)),[&](std::size_t i,std::size_t j) {
            const auto& a=tris[i];const auto& b=tris[j];
            if(cross && a.chart==b.chart) return;
            ++report.candidate_pairs;++report.exact_tests;
            if(separated_by_edges(a.p,b.p,&report.exact_fallbacks) || separated_by_edges(b.p,a.p,&report.exact_fallbacks)) return;
            const double area=triangle_intersection_area(a.p,b.p);
            if(!(area>options.area_tolerance*std::min(a.area,b.area))) {++report.contacts;return;}
            ++report.overlaps;++report.counts[a.chart];if(cross)++report.counts[b.chart];
            faces.push_back(a.face);faces.push_back(b.face);
            if(options.retain_pairs) {
                const bool order=a.id<b.id;const auto& f=order?a:b;const auto& s=order?b:a;
                retain({f.id,s.id,f.face,s.face,f.chart,s.chart,area});
            }
        });
        report.cell_entries+=tris.size();
    };
    if(cross) {
        std::vector<OverlapTriangle> all;for(auto ch:report.selected_charts) collect(ch,all);test_group(all);
    } else for(auto ch:report.selected_charts) {std::vector<OverlapTriangle> tris;collect(ch,tris);test_group(tris);}
    std::sort(faces.begin(),faces.end());faces.erase(std::unique(faces.begin(),faces.end()),faces.end());report.overlapping_faces=std::move(faces);
    std::sort(report.pairs.begin(),report.pairs.end(),by_ids);
    if(report.pairs.size()>options.max_pairs) {report.pairs.resize(options.max_pairs);report.pairs_truncated=true;}
    return report;
}
}
int orientation_sign(Vec2 a,Vec2 b,Vec2 c,bool* exact) {
    if(exact) *exact=false;
    require(finite(a) && finite(b) && finite(c),"quality.orientation","orientation needs finite points");
    const double left=(b[0]-a[0])*(c[1]-a[1]),right=(c[0]-a[0])*(b[1]-a[1]),det=left-right;
    const double bound=orientation_bound*(std::abs(left)+std::abs(right));
    if(std::isfinite(det) && std::abs(det)>bound) return sign_of(det);
    if(exact) *exact=true;
    return exact_orientation(a,b,c);
}
bool triangle_interiors_intersect(const std::array<Vec2,3>& first,const std::array<Vec2,3>& second,std::uint64_t* exact_fallbacks) {
    const auto a=counter_clockwise(first,exact_fallbacks),b=counter_clockwise(second,exact_fallbacks);
    if(!a || !b) return false;
    return !separated_by_edges(*a,*b,exact_fallbacks) && !separated_by_edges(*b,*a,exact_fallbacks);
}
double triangle_intersection_area(const std::array<Vec2,3>& first,const std::array<Vec2,3>& second) {
    auto ccw=[](std::array<Vec2,3> t) {if(twice_area(t)<0) std::swap(t[1],t[2]);return t;};
    const auto a=ccw(first),b=ccw(second);
    if(twice_area(a)==0 || twice_area(b)==0) return 0;
    // Sutherland-Hodgman: clip b by the three half-planes of a.
    std::vector<Vec2> polygon(b.begin(),b.end()),next;
    for(std::size_t i=0;i<3 && !polygon.empty();++i) {
        const auto p=a[i],q=a[(i+1)%3];next.clear();
        for(std::size_t k=0;k<polygon.size();++k) {
            const auto s=polygon[k],e=polygon[(k+1)%polygon.size()];
            const double ds=orient(p,q,s),de=orient(p,q,e);
            if(ds>=0) next.push_back(s);
            if((ds>=0)!=(de>=0)) {const double t=ds/(ds-de);next.push_back({s[0]+t*(e[0]-s[0]),s[1]+t*(e[1]-s[1])});}
        }
        polygon.swap(next);
    }
    double area=0;
    for(std::size_t k=0;k<polygon.size();++k) {const auto& s=polygon[k];const auto& e=polygon[(k+1)%polygon.size()];area+=s[0]*e[1]-e[0]*s[1];}
    return std::max(0.0,0.5*area);
}
OverlapReport chart_overlaps(const QualityGeometry& geometry,CornerUvView uv,const OverlapOptions& options,std::optional<std::span<const MeshIndex>> charts) {
    return overlaps(geometry,uv,options,charts,false);
}
OverlapReport cross_chart_overlaps(const QualityGeometry& geometry,CornerUvView uv,const OverlapOptions& options,std::optional<std::span<const MeshIndex>> charts) {
    return overlaps(geometry,uv,options,charts,true);
}
nlohmann::json OverlapReport::to_json() const {
    auto pair_json=nlohmann::json::array();
    for(const auto& p:pairs) pair_json.push_back({{"triangles",{p.first_triangle,p.second_triangle}},{"faces",{p.first_face,p.second_face}},
        {"charts",{p.first_chart,p.second_chart}},{"area",p.area}});
    auto tested_json=nlohmann::json::array();for(auto t:tested) tested_json.push_back(t!=0);
    return {{"policy",policy},{"cross_chart",cross_chart},{"counts",counts},{"tested",tested_json},{"overlapping_faces",overlapping_faces},
        {"pairs",pair_json},{"pairs_truncated",pairs_truncated},{"triangles",triangles},{"degenerate",degenerate},{"untestable",untestable},
        {"cells",cells},{"cell_entries",cell_entries},{"candidate_pairs",candidate_pairs},{"exact_tests",exact_tests},
        {"exact_fallbacks",exact_fallbacks},{"contacts",contacts},{"overlaps",overlaps}};
}
}
