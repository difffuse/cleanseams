// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/seam_metrics.py (all functions).
#include "cleanseams/metrics.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/seam_graph.hpp"
#include "cleanseams/detail/numpy_blas.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <numbers>
#include <numeric>
namespace cleanseams {
namespace {
constexpr double nan_value=std::numeric_limits<double>::quiet_NaN();
void require(bool condition,const char* field,const char* message) {
    if(!condition)throw ContractError(ErrorCode::invalid_input,field,message);
}
void check_seams(const MeshSnapshot& m,std::span<const std::uint8_t> seams,const char* field) {
    require(seams.size()==m.tables().edge_vertices.size(),field,"one flag per mesh edge required");
}
void check_uv(const MeshSnapshot& m,std::span<const Vec2> uv) {
    require(uv.size()==m.input().corner_vertices.size(),"metrics.uv","one UV per original corner required");
}
bool finite(const Vec2& p) {return std::isfinite(p[0]) && std::isfinite(p[1]);}
double sum_where(std::span<const double> values,auto&& keep) {
    std::vector<double> kept;for(std::size_t i=0;i<values.size();++i)if(keep(i))kept.push_back(values[i]);
    return detail::numpy_sum(kept);
}
// np.argmin / np.argmax: the first NaN if any, else the first extreme.
template<class Better> std::size_t arg_extreme(std::span<const double> values,Better better) {
    for(std::size_t i=0;i<values.size();++i)if(std::isnan(values[i]))return i;
    std::size_t best=0;for(std::size_t i=1;i<values.size();++i)if(better(values[i],values[best]))best=i;
    return best;
}
nlohmann::json number(double x) {return nlohmann::json(x);}
}
double seam_length(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    check_seams(m,seams,"metrics.seams");
    return sum_where(m.geometry().edge_length,[&](std::size_t e) {return seams[e]!=0;});
}
double seam_length_norm(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    const double area=detail::numpy_sum(m.geometry().face_area);
    return area>0?seam_length(m,seams)/std::sqrt(area):nan_value;
}
double seam_length_ratio(const MeshSnapshot& m,std::span<const std::uint8_t> seams,std::span<const std::uint8_t> reference) {
    const double ref=seam_length(m,reference);
    return ref>0?seam_length(m,seams)/ref:nan_value;
}
std::vector<std::uint8_t> seams_from_uv(const MeshSnapshot& m,std::span<const Vec2> uv,double eps) {
    check_uv(m,uv);
    const auto& t=m.tables();const auto& vertex=m.input().corner_vertices;
    std::vector<std::uint8_t> seams(t.edge_vertices.size(),0);
    // max(|p - q|) over both coordinates; a NaN coordinate makes the row NaN (never a seam).
    auto apart=[&](MeshIndex p,MeshIndex q) {
        const double d0=std::abs(uv[p][0]-uv[q][0]),d1=std::abs(uv[p][1]-uv[q][1]);
        return !std::isnan(d0) && !std::isnan(d1) && std::max(d0,d1)>eps;
    };
    for(std::size_t e=0;e<seams.size();++e) {
        const auto corners=t.edge_corners[e];if(corners.size()!=2)continue;
        const auto c0=corners[0],c1=corners[1],n0=m.corner_next(c0),n1=m.corner_next(c1);
        const bool consistent=vertex[c0]==vertex[n1];
        const auto pa=consistent?n1:c1,pb=consistent?c1:n1;
        seams[e]=apart(c0,pa) || apart(n0,pb);
    }
    return seams;
}
IslandLabels islands_from_seams(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    check_seams(m,seams,"metrics.seams");
    const auto& t=m.tables();const auto nf=m.input().face_offsets.size()-1;
    std::vector<MeshIndex> parent(nf);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](MeshIndex x) {while(parent[x]!=x)x=parent[x]=parent[parent[x]];return x;};
    for(std::size_t e=0;e<seams.size();++e) {
        if(seams[e] || t.edge_faces[e].size()!=2)continue;
        const auto a=root(t.edge_faces[e][0]),b=root(t.edge_faces[e][1]);
        if(a!=b)parent[std::max(a,b)]=std::min(a,b);
    }
    IslandLabels out;out.labels.resize(nf);std::vector<MeshIndex> label(nf,std::numeric_limits<MeshIndex>::max());
    for(MeshIndex f=0;f<nf;++f) {const auto r=root(f);if(label[r]==std::numeric_limits<MeshIndex>::max())label[r]=out.count++;out.labels[f]=label[r];}
    return out;
}
std::vector<std::int64_t> seam_hop_distance(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    check_seams(m,seams,"metrics.seams");
    const auto& t=m.tables();std::vector<std::int64_t> dist(m.input().positions.size(),-1);
    std::vector<MeshIndex> start;
    for(std::size_t e=0;e<seams.size();++e)if(seams[e])for(auto v:t.edge_vertices[e])start.push_back(v);
    std::sort(start.begin(),start.end());start.erase(std::unique(start.begin(),start.end()),start.end());
    std::deque<MeshIndex> queue;for(auto v:start) {dist[v]=0;queue.push_back(v);}
    while(!queue.empty()) {
        const auto v=queue.front();queue.pop_front();
        for(auto e:t.vert_edges[v]) {const auto w=m.edge_other(e,v);if(dist[w]<0) {dist[w]=dist[v]+1;queue.push_back(w);}}
    }
    return dist;
}
std::vector<double> edge_hops(const MeshSnapshot& m,std::span<const std::uint8_t>,std::span<const std::uint8_t> target) {
    const auto d=seam_hop_distance(m,target);const auto& t=m.tables();
    std::vector<double> hops(t.edge_vertices.size());
    for(std::size_t e=0;e<hops.size();++e) {
        const auto a=d[t.edge_vertices[e][0]],b=d[t.edge_vertices[e][1]];
        hops[e]=a<0 || b<0?static_cast<double>(m.input().positions.size()):static_cast<double>(std::max(a,b));
        if(target[e])hops[e]=0;
    }
    return hops;
}
nlohmann::json SeamAgreement::to_json() const {
    return {{"exact",number(exact)},{"one_ring",number(one_ring)},{"mean_hops",number(mean_hops)},{"exact_of_ours",number(exact_of_ours)},
            {"one_ring_of_ours",number(one_ring_of_ours)},{"mean_hops_of_ours",number(mean_hops_of_ours)},
            {"island_count_ratio",number(island_count_ratio)},{"islands",islands},{"reference_islands",reference_islands}};
}
SeamAgreement seam_agreement(const MeshSnapshot& m,std::span<const std::uint8_t> seams,std::span<const std::uint8_t> reference) {
    check_seams(m,seams,"metrics.seams");check_seams(m,reference,"metrics.reference");
    const auto& L=m.geometry().edge_length;
    auto block=[&](std::span<const std::uint8_t> a,std::span<const std::uint8_t> b) {
        const double la=sum_where(L,[&](std::size_t e) {return a[e]!=0;});
        if(!(la>0))return std::array<double,3>{nan_value,nan_value,nan_value};
        const auto hops=edge_hops(m,a,b);
        std::vector<double> weighted;for(std::size_t e=0;e<L.size();++e)if(a[e])weighted.push_back(hops[e]*L[e]);
        return std::array<double,3>{100.0*sum_where(L,[&](std::size_t e) {return a[e] && b[e];})/la,
                                    100.0*sum_where(L,[&](std::size_t e) {return a[e] && hops[e]<=1;})/la,
                                    detail::numpy_sum(weighted)/la};
    };
    const auto [e,r,h]=block(reference,seams);const auto [e2,r2,h2]=block(seams,reference);
    SeamAgreement out;out.exact=e;out.one_ring=r;out.mean_hops=h;out.exact_of_ours=e2;out.one_ring_of_ours=r2;out.mean_hops_of_ours=h2;
    out.islands=islands_from_seams(m,seams).count;out.reference_islands=islands_from_seams(m,reference).count;
    out.island_count_ratio=out.reference_islands?static_cast<double>(out.islands)/static_cast<double>(out.reference_islands):nan_value;
    return out;
}
std::uint64_t dangling_seam_ends(const MeshSnapshot& m,std::span<const std::uint8_t> seams) {
    check_seams(m,seams,"metrics.seams");
    const auto& t=m.tables();std::vector<std::uint64_t> degree(m.input().positions.size(),0);
    for(std::size_t e=0;e<seams.size();++e)if(seams[e])for(auto v:t.edge_vertices[e])++degree[v];
    std::uint64_t count=0;for(std::size_t v=0;v<degree.size();++v)count+=degree[v]==1 && !t.vert_boundary[v];
    return count;
}
std::vector<Vec2> convex_hull(std::span<const Vec2> points) {
    std::vector<Vec2> P;for(const auto& p:points)if(finite(p))P.push_back(p);
    if(P.size()<3)return P;
    // np.unique(P, axis=0): lexicographically sorted distinct rows.
    std::sort(P.begin(),P.end(),[](const Vec2& a,const Vec2& b) {return a[0]!=b[0]?a[0]<b[0]:a[1]<b[1];});
    P.erase(std::unique(P.begin(),P.end(),[](const Vec2& a,const Vec2& b) {return a[0]==b[0] && a[1]==b[1];}),P.end());
    if(P.size()<3)return P;
    auto cross=[](const Vec2& o,const Vec2& a,const Vec2& b) {return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);};
    std::vector<Vec2> lower,upper;
    for(const auto& p:P) {while(lower.size()>=2 && cross(lower[lower.size()-2],lower.back(),p)<=0)lower.pop_back();lower.push_back(p);}
    for(auto it=P.rbegin();it!=P.rend();++it) {while(upper.size()>=2 && cross(upper[upper.size()-2],upper.back(),*it)<=0)upper.pop_back();upper.push_back(*it);}
    lower.pop_back();upper.pop_back();lower.insert(lower.end(),upper.begin(),upper.end());
    return lower;
}
MinAreaRect min_area_rect(std::span<const Vec2> points) {
    const auto H=convex_hull(points);
    if(H.empty())return {};
    if(H.size()<3) {
        Vec2 lo=H[0],hi=H[0];
        for(const auto& p:H)for(std::size_t k=0;k<2;++k) {lo[k]=std::min(lo[k],p[k]);hi[k]=std::max(hi[k],p[k]);}
        return {hi[0]-lo[0],hi[1]-lo[1],0.0};
    }
    const auto n=H.size();std::vector<Vec2> U,W;
    for(std::size_t i=0;i<n;++i) {
        const Vec2 e{H[(i+1)%n][0]-H[i][0],H[(i+1)%n][1]-H[i][1]};
        const double length=std::sqrt(e[0]*e[0]+e[1]*e[1]);
        if(length>0) {U.push_back({e[0]/length,e[1]/length});W.push_back({-U.back()[1],U.back()[0]});}
    }
    if(U.empty())return {};
    // H @ U.T and H @ W.T through NumPy/Accelerate (square when every hull edge is kept).
    const auto px=detail::accelerate_square_projection(H,U),py=detail::accelerate_square_projection(H,W);
    const auto h=U.size();std::vector<double> width(h),height(h),area(h);
    for(std::size_t j=0;j<h;++j) {
        double x0=px[j],x1=px[j],y0=py[j],y1=py[j];
        for(std::size_t i=1;i<n;++i) {x0=std::min(x0,px[i*h+j]);x1=std::max(x1,px[i*h+j]);y0=std::min(y0,py[i*h+j]);y1=std::max(y1,py[i*h+j]);}
        width[j]=x1-x0;height[j]=y1-y0;area[j]=width[j]*height[j];
    }
    const auto i=arg_extreme(area,std::less<>());
    return {width[i],height[i],std::atan2(U[i][1],U[i][0])};
}
SkylinePack skyline_pack(std::span<const std::array<double,2>> sizes,double width,bool allow_rotate) {
    struct Segment {double x,y,w;};
    auto seg_less=[](const Segment& a,const Segment& b) {
        if(a.x!=b.x)return a.x<b.x;if(a.y!=b.y)return a.y<b.y;return a.w<b.w;
    };
    const auto n=sizes.size();std::vector<std::size_t> order(n);std::iota(order.begin(),order.end(),0);
    auto key=[&](std::size_t i) {
        const auto& s=sizes[i];const double big=s[0]>=s[1]?s[0]:s[1],small=s[0]<=s[1]?s[0]:s[1];
        return std::array<double,2>{-big,-small};
    };
    std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
        const auto ka=key(a),kb=key(b);
        if(ka[0]!=kb[0])return ka[0]<kb[0];if(ka[1]!=kb[1])return ka[1]<kb[1];return a<b;
    });
    std::vector<Segment> skyline{{0.0,0.0,width}};
    SkylinePack out;out.placements.resize(n);
    for(auto i:order) {
        bool have=false;double best_top=0,best_x=0,best_w=0,best_h=0;bool best_rot=false;
        for(int r=0;r<(allow_rotate?2:1);++r) {
            const bool rot=r==1;const double w=rot?sizes[i][1]:sizes[i][0],h=rot?sizes[i][0]:sizes[i][1];
            for(std::size_t j=0;j<skyline.size();++j) {
                const auto [sx,sy,sw]=skyline[j];(void)sw;
                if(sx+w>width+1e-12 && j>0)continue;
                double y=sy;std::size_t k=j;
                while(k<skyline.size() && skyline[k].x<sx+w-1e-12) {y=y>=skyline[k].y?y:skyline[k].y;++k;}
                const double top=y+h;
                if(!have || top<best_top || (top==best_top && sx<best_x)) {have=true;best_top=top;best_x=sx;best_rot=rot;best_w=w;best_h=h;}
            }
        }
        const double top=best_top,x=best_x,w=best_w,h=best_h;
        out.placements[i]={x,top-h,best_rot};
        out.height=out.height>=top?out.height:top;
        std::vector<Segment> next;
        for(const auto& s:skyline) {
            const double ex=s.x+s.w;
            if(ex<=x+1e-12 || s.x>=x+w-1e-12) {next.push_back(s);continue;}
            if(s.x<x)next.push_back({s.x,s.y,x-s.x});
            if(ex>x+w)next.push_back({x+w,s.y,ex-x-w});
        }
        next.push_back({x,top,w});
        std::stable_sort(next.begin(),next.end(),seg_less);
        std::vector<Segment> merged;
        for(const auto& s:next) {
            if(!merged.empty() && std::abs(merged.back().y-s.y)<1e-12 && std::abs(merged.back().x+merged.back().w-s.x)<1e-9)
                merged.back().w=merged.back().w+s.w;
            else merged.push_back(s);
        }
        skyline=std::move(merged);
    }
    return out;
}
std::vector<IslandRect> island_rects(const CutSurface& cs,std::span<const Vec2> uv) {
    const auto& m=cs.mesh();check_uv(m,uv);const auto& t=cs.tables();
    std::vector<IslandRect> out;
    for(MeshIndex ch=0;ch<t.n_charts;++ch) {
        std::vector<Vec2> pts;for(auto f:t.chart_faces[ch])for(auto c:m.face_corners(f))pts.push_back(uv[c]);
        const auto rect=min_area_rect(pts);double area=0.0;
        for(auto f:t.chart_faces[ch]) {
            std::vector<double> x,y,xr,yr;bool ok=true;
            for(auto c:m.face_corners(f)) {x.push_back(uv[c][0]);y.push_back(uv[c][1]);ok=ok && finite(uv[c]);}
            if(!ok)continue;
            xr.assign(x.begin()+1,x.end());xr.push_back(x[0]);yr.assign(y.begin()+1,y.end());yr.push_back(y[0]);
            // np.dot of a strided column with a rolled copy (Accelerate strided ddot).
            area+=0.5*std::abs(detail::accelerate_strided_dot(x,yr)-detail::accelerate_strided_dot(xr,y));
        }
        out.push_back({rect.width,rect.height,area});
    }
    return out;
}
double packing_efficiency(const CutSurface& cs,std::span<const Vec2> uv,double padding,std::int64_t iterations) {
    const auto rects=island_rects(cs,uv);
    if(rects.empty() || std::any_of(rects.begin(),rects.end(),[](const IslandRect& r) {return !std::isfinite(r.width) || !std::isfinite(r.height);}))
        return nan_value;
    double area=0;for(const auto& r:rects)area+=r.area;
    if(!(area>0))return nan_value;
    auto fits=[&](double k) {
        std::vector<std::array<double,2>> inflated;
        for(const auto& r:rects) {
            const std::array<double,2> s{k*r.width+2*padding,k*r.height+2*padding};
            if(std::min(s[0],s[1])>1.0+1e-12)return false;
            inflated.push_back(s);
        }
        return skyline_pack(inflated,1.0).height<=1.0+1e-12;
    };
    double total=0;for(const auto& r:rects)total+=r.width*r.height;
    double hi=total>0?1.0/std::sqrt(total):1.0;
    if(fits(hi))return std::min(1.0,hi*hi*area);
    double lo=0.0;
    for(std::int64_t i=0;i<iterations;++i) {const double mid=0.5*(lo+hi);if(fits(mid))lo=mid;else hi=mid;}
    return lo*lo*area;
}
TriangleDistortion triangle_distortion(const std::shared_ptr<const CutSurface>& cs,std::span<const Vec2> uv) {
    require(bool(cs),"metrics.cut","cut surface owner required");check_uv(cs->mesh(),uv);
    const QualitySnapshot snapshot(quality_geometry(cs),CornerUvView{uv});
    const auto& tris=snapshot.triangles();const auto n=cs->tables().n_charts;
    TriangleDistortion out;
    std::vector<double> d3(n,0.0),duv(n,0.0);std::vector<std::uint8_t> ok(tris.size());
    for(std::size_t i=0;i<tris.size();++i) {
        const auto& t=tris[i];out.chart.push_back(t.chart);out.area3.push_back(t.area3);
        ok[i]=finite(t.singular_values) && std::isfinite(t.area2);
        d3[t.chart]+=ok[i]?t.area3:0.0;duv[t.chart]+=ok[i]?std::abs(t.area2):0.0;
    }
    std::vector<double> dens(n);for(std::size_t ch=0;ch<n;++ch)dens[ch]=d3[ch]>0?duv[ch]/std::max(d3[ch],1e-300):1.0;
    // np.maximum propagates NaN.
    auto maximum=[](double a,double b) {return std::isnan(a) || std::isnan(b)?nan_value:std::max(a,b);};
    for(std::size_t i=0;i<tris.size();++i) {
        const auto& sv=tris[i].singular_values;
        const double aniso=sv[0]/maximum(sv[1],1e-300);
        const double r=sv[0]*sv[1]/maximum(dens[tris[i].chart],1e-300);
        const double area=maximum(r,1.0/maximum(r,1e-300));
        out.anisotropy.push_back(ok[i] && std::isfinite(aniso)?aniso:1e9);
        out.area.push_back(ok[i] && std::isfinite(area)?area:1e9);
    }
    return out;
}
nlohmann::json ChartWorst::to_json() const {
    return {{"chart_worst_id",chart_worst_id},{"chart_worst_aniso_p95",number(chart_worst_aniso_p95)},
            {"chart_worst_area_p95",number(chart_worst_area_p95)},{"per_chart_aniso_p95",per_chart_aniso_p95},
            {"per_chart_area_p95",per_chart_area_p95}};
}
namespace {
ChartWorst chart_worst_of(const TriangleDistortion& d,MeshIndex n) {
    ChartWorst out;out.per_chart_aniso_p95.assign(n,nan_value);out.per_chart_area_p95.assign(n,nan_value);
    std::vector<std::vector<double>> aniso(n),area(n),weight(n);
    for(std::size_t i=0;i<d.chart.size();++i) {aniso[d.chart[i]].push_back(d.anisotropy[i]);area[d.chart[i]].push_back(d.area[i]);weight[d.chart[i]].push_back(d.area3[i]);}
    for(MeshIndex ch=0;ch<n;++ch) {
        if(aniso[ch].empty())continue;
        out.per_chart_aniso_p95[ch]=WeightedDistribution(aniso[ch],weight[ch]).percentile(95);
        out.per_chart_area_p95[ch]=WeightedDistribution(area[ch],weight[ch]).percentile(95);
    }
    std::int64_t worst=-1;
    for(MeshIndex ch=0;ch<n;++ch) {
        const double v=out.per_chart_aniso_p95[ch];
        if(!std::isnan(v) && (worst<0 || v>out.per_chart_aniso_p95[static_cast<std::size_t>(worst)]))worst=ch;
    }
    if(worst>=0) {
        out.chart_worst_id=worst;out.chart_worst_aniso_p95=out.per_chart_aniso_p95[static_cast<std::size_t>(worst)];
        out.chart_worst_area_p95=out.per_chart_area_p95[static_cast<std::size_t>(worst)];
    }
    return out;
}
Hotspot hotspot_of(const CutSurface& cs,const TriangleDistortion& d,std::int64_t rings,double hot) {
    Hotspot out;if(d.area3.empty())return out;
    const auto& m=cs.mesh();const auto& tri_face=m.triangulation().faces;const auto nf=m.input().face_offsets.size()-1;
    std::vector<double> fa(nf,0.0),fx(nf,0.0),fr(nf,0.0);
    for(std::size_t i=0;i<d.area3.size();++i) {
        fa[tri_face[i]]+=d.area3[i];fx[tri_face[i]]+=d.area3[i]*std::min(d.anisotropy[i],1e6);fr[tri_face[i]]+=d.area3[i]*std::min(d.area[i],1e6);
    }
    std::vector<double> face_aniso(nf),face_area(nf);
    for(std::size_t f=0;f<nf;++f) {face_aniso[f]=fa[f]>0?fx[f]/std::max(fa[f],1e-300):1.0;face_area[f]=fa[f]>0?fr[f]/std::max(fa[f],1e-300):1.0;}
    auto wa=windowed_face_values(cs,face_aniso,fa,rings),wr=windowed_face_values(cs,face_area,fa,rings);
    for(auto& v:wa)if(!std::isfinite(v))v=1.0;
    for(auto& v:wr)if(!std::isfinite(v))v=1.0;
    const auto f=arg_extreme(wa,std::greater<>());const auto ch=cs.tables().chart_of_face[f];
    const double total=detail::numpy_sum(fa);
    const double chart_area=sum_where(fa,[&](std::size_t g) {return cs.tables().chart_of_face[g]==ch;});
    out.hotspot_max=wa[f];out.hotspot_area_max=wr[arg_extreme(wr,std::greater<>())];
    out.hotspot_face=static_cast<std::int64_t>(f);out.hotspot_chart=ch;
    out.hotspot_area_pct=total>0?100.0*sum_where(fa,[&](std::size_t g) {return wa[g]>hot;})/total:nan_value;
    out.hotspot_chart_area_pct=chart_area>0?100.0*sum_where(fa,[&](std::size_t g) {return wa[g]>hot && cs.tables().chart_of_face[g]==ch;})/chart_area:nan_value;
    out.windowed_anisotropy=std::move(wa);
    return out;
}
}
ChartWorst chart_worst(const std::shared_ptr<const CutSurface>& cs,std::span<const Vec2> uv) {
    return chart_worst_of(triangle_distortion(cs,uv),cs->tables().n_charts);
}
std::vector<double> windowed_face_values(const CutSurface& cs,std::span<const double> values,std::span<const double> weights,std::int64_t rings) {
    const auto& m=cs.mesh();const auto nf=m.input().face_offsets.size()-1;
    require(values.size()==nf && weights.size()==nf,"metrics.window","one value and weight per face required");
    const auto& sec=cs.tables().sector_of_corner;const auto& cf=m.tables().corner_face;const auto ns=cs.tables().n_sectors;
    std::vector<double> X(nf),A(weights.begin(),weights.end());
    for(std::size_t f=0;f<nf;++f)X[f]=values[f]*weights[f];
    for(std::int64_t r=0;r<rings;++r) {
        std::vector<double> sx(ns,0.0),sa(ns,0.0);
        for(std::size_t c=0;c<cf.size();++c) {sx[sec[c]]+=X[cf[c]];sa[sec[c]]+=A[cf[c]];}
        std::vector<double> nx(nf,0.0),na(nf,0.0);
        for(std::size_t c=0;c<cf.size();++c) {nx[cf[c]]+=sx[sec[c]];na[cf[c]]+=sa[sec[c]];}
        X=std::move(nx);A=std::move(na);
    }
    std::vector<double> out(nf);
    for(std::size_t f=0;f<nf;++f)out[f]=A[f]>0?X[f]/std::max(A[f],1e-300):nan_value;
    return out;
}
nlohmann::json Hotspot::to_json() const {
    return {{"hotspot_max",number(hotspot_max)},{"hotspot_area_max",number(hotspot_area_max)},{"hotspot_face",hotspot_face},
            {"hotspot_chart",hotspot_chart},{"hotspot_area_pct",number(hotspot_area_pct)},{"hotspot_chart_area_pct",number(hotspot_chart_area_pct)}};
}
Hotspot hotspot(const std::shared_ptr<const CutSurface>& cs,std::span<const Vec2> uv,std::int64_t rings,double hot) {
    return hotspot_of(*cs,triangle_distortion(cs,uv),rings,hot);
}
nlohmann::json chain_statistics(std::span<const std::uint8_t> seams,const std::shared_ptr<const CutSurface>& cs) {
    require(bool(cs),"metrics.cut","cut surface owner required");
    const auto& m=cs->mesh();check_seams(m,seams,"metrics.seams");const auto& t=cs->tables();const auto& g=m.geometry();
    const SeamGraph graph(cs->mesh_owner(),seams,{},{},cs);
    std::uint64_t lane=0,turns=0;std::int64_t edges=0;
    for(const auto& chain:graph.chains()) {
        const std::vector<TransferIndex> verts(chain.verts.begin(),chain.verts.end());
        const auto result=chain_turns(m,verts);lane+=result.lane_changes;turns+=result.turns;edges+=result.edges;
    }
    const double total=std::max(detail::numpy_sum(g.face_area),1e-12);
    std::vector<double> chart_area(t.n_charts,0.0),perimeter(t.n_charts,0.0);
    for(std::size_t f=0;f<t.chart_of_face.size();++f)chart_area[t.chart_of_face[f]]+=g.face_area[f];
    for(auto c:t.boundary_corners)perimeter[t.chart_of_corner[c]]+=g.edge_length[m.tables().corner_edge[c]];
    std::vector<double> areas(t.n_charts),iso(t.n_charts);
    for(std::size_t ch=0;ch<t.n_charts;++ch) {
        areas[ch]=chart_area[ch]/total;iso[ch]=perimeter[ch]*perimeter[ch]/(4*std::numbers::pi*std::max(chart_area[ch],1e-12));
    }
    const double sharp_angle=60*(std::numbers::pi/180);
    std::uint64_t seam_edges=0,sharp=0,both=0,tiny=0;
    for(std::size_t e=0;e<seams.size();++e) {
        const bool s=std::abs(g.edge_dihedral[e])>=sharp_angle;seam_edges+=seams[e]!=0;sharp+=s;both+=s && seams[e];
    }
    for(double a:areas)tiny+=a<0.002;
    const double n_seam=static_cast<double>(std::max<std::uint64_t>(seam_edges,1));
    nlohmann::json out{{"chains",graph.chains().size()},{"seam_edges",seam_edges},{"lane_changes",lane},{"turns_over_30deg",turns},
        {"lane_changes_per_100_edges",detail::python_round(100.0*static_cast<double>(lane)/static_cast<double>(std::max<std::int64_t>(edges,1)),2)},
        {"islands",t.n_charts},{"tiny_islands_under_0.2pct",tiny},
        {"seam_on_sharp_edges_pct",detail::python_round(100.0*static_cast<double>(both)/n_seam,1)},
        {"sharp_edges_seamed_pct",detail::python_round(100.0*static_cast<double>(both)/static_cast<double>(std::max<std::uint64_t>(sharp,1)),1)}};
    if(areas.empty()) {out["island_area_min_pct"]=0;out["island_area_max_pct"]=0;out["isoperimetric_median"]=0;out["isoperimetric_max"]=0;}
    else {
        out["island_area_min_pct"]=detail::python_round(*std::min_element(areas.begin(),areas.end())*100,3);
        out["island_area_max_pct"]=detail::python_round(*std::max_element(areas.begin(),areas.end())*100,1);
        out["isoperimetric_median"]=detail::python_round(detail::numpy_median(iso),2);
        out["isoperimetric_max"]=detail::python_round(*std::max_element(iso.begin(),iso.end()),2);
    }
    return out;
}
nlohmann::json seam_report(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> seams,std::shared_ptr<const CutSurface> cs,
                           SeamReportContext* context,std::optional<std::span<const std::uint8_t>> reference,std::optional<std::span<const Vec2>> uv,
                           const SeamReportOptions& options,const SymmetryService& symmetry) {
    const auto start=std::chrono::steady_clock::now();
    require(bool(mesh),"metrics.mesh","mesh snapshot required");const auto& m=*mesh;check_seams(m,seams,"metrics.seams");
    if(!cs)cs=cut_surface(mesh,seams);
    require(cs->mesh().same_geometry(m),"metrics.cut","cut surface belongs to another mesh");
    SeamReportContext local;auto& meta=context?*context:local;
    auto out=chain_statistics(seams,cs);
    out["seam_length"]=seam_length(m,seams);out["seam_length_norm"]=number(seam_length_norm(m,seams));
    if(reference) {
        out["seam_length_ratio"]=number(seam_length_ratio(m,seams,*reference));out["agreement"]=seam_agreement(m,seams,*reference).to_json();
    } else {out["seam_length_ratio"]=number(nan_value);out["agreement"]=nullptr;}
    if(!meta.exposure)meta.exposure=std::make_shared<const ExposureField>(exposure_field(m,options.K));
    const auto exposure=seam_exposure(m,seams,*meta.exposure,options.hidden_below);
    out["seam_exposure"]=number(exposure.seam_exposure);out["mesh_exposure"]=number(exposure.mesh_exposure);
    out["exposure_ratio"]=number(exposure.exposure_ratio);out["hidden_seam_pct"]=number(exposure.hidden_seam_pct);
    // Mirror plane (symmetry.detect_plane / seam_symmetry): the injected service, else the
    // native T14 module. The source always detects the plane.
    static const SymmetryService native=native_symmetry_service();
    const auto& service=symmetry.available()?symmetry:native;
    std::optional<DetectedPlane> plane;
    if(meta.plane)plane=*meta.plane;
    else {plane=service.detect_plane(m,options.symmetry_tol);meta.plane=plane;}
    out["mesh_symmetry_fidelity"]=plane?plane->score:0.0;
    out["symmetry"]=number(plane?service.seam_symmetry(m,seams,*plane):nan_value);
    out["symmetry_plane"]=plane?plane->description:nlohmann::json(nullptr);
    const auto& t=cs->tables();
    std::vector<std::uint64_t> n_faces(t.n_charts,0);for(auto ch:t.chart_of_face)++n_faces[ch];
    out["islands_under_10_faces"]=static_cast<std::uint64_t>(std::count_if(n_faces.begin(),n_faces.end(),[](auto n) {return n<10;}));
    out["non_disk_charts"]=static_cast<std::uint64_t>(std::count(t.chart_is_disk.begin(),t.chart_is_disk.end(),0));
    out["dangling_seam_ends"]=dangling_seam_ends(m,seams);
    if(uv) {
        check_uv(m,*uv);
        out["packing_efficiency"]=number(packing_efficiency(*cs,*uv));
        const auto distortion=triangle_distortion(cs,*uv);
        const auto worst=chart_worst_of(distortion,t.n_charts);
        out["chart_worst_id"]=worst.chart_worst_id;out["chart_worst_aniso_p95"]=number(worst.chart_worst_aniso_p95);
        out["chart_worst_area_p95"]=number(worst.chart_worst_area_p95);
        const auto hot=hotspot_of(*cs,distortion,2,options.hot).to_json();
        for(auto& [key,value]:hot.items())out[key]=value;
    } else {
        out["packing_efficiency"]=number(nan_value);out["chart_worst_id"]=-1;out["chart_worst_aniso_p95"]=number(nan_value);
        out["chart_worst_area_p95"]=number(nan_value);out["hotspot_max"]=number(nan_value);out["hotspot_area_max"]=number(nan_value);
        out["hotspot_face"]=-1;out["hotspot_chart"]=-1;out["hotspot_area_pct"]=number(nan_value);out["hotspot_chart_area_pct"]=number(nan_value);
    }
    out["report_time_s"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return out;
}
// ------------------------------------------------------------------ Phase B (T32)
nlohmann::json SectorDistortionField::units() {
    return {{"log_area","natural log of the UV/3D area ratio relative to the chart density; signed (negative compressed, positive expanded); "
                        "area units: a uniform linear scale k gives 2 ln k, so ln(1.25) is a 25% area change and an 11.8% linear change"},
            {"log_anisotropy","natural log of the singular value ratio s0/s1 (shape units, >= 0)"},
            {"mass","a third of the 3D area of each measured triangle, per incident cut-surface sector"}};
}
nlohmann::json SectorDistortionField::summary() const {
    double total=0;for(double m:mass)total+=m;
    double area=0;for(double a:chart_area)area+=a;
    return {{"version",version},{"smoothing_steps",smoothing_steps},{"charts",charts},{"sectors",mass.size()},{"triangles",triangles},
            {"mass",total},{"measured_area",area},{"mass_error",mass_error},{"moment_error",moment_error},{"units",units()}};
}
SectorDistortionField sector_distortion_field(const QualitySnapshot& snapshot,std::int64_t smoothing_steps) {
    require(smoothing_steps>=0,"relief.smoothing_steps","smoothing steps must be nonnegative");
    const auto& cs=snapshot.geometry().cut();const auto& ct=cs.tables();const auto ns=ct.n_sectors;const auto& tris=snapshot.triangles();
    SectorDistortionField out;out.smoothing_steps=smoothing_steps;
    out.mass.assign(ns,0.0);out.chart_area.assign(ct.n_charts,0.0);out.chart_mass.assign(ct.n_charts,0.0);
    out.log_area.assign(ns,nan_value);out.log_anisotropy.assign(ns,nan_value);
    // The distortion report's measured triangles and per-chart density.
    auto usable=[](const TriangleQuality& t) {
        return t.finite && !t.sliver && !t.exempt && std::isfinite(t.singular_values[0]) && std::isfinite(t.singular_values[1]) &&
               t.singular_values[1]>0 && t.area3>0;
    };
    std::vector<double> d3(ct.n_charts,0.0),duv(ct.n_charts,0.0);
    for(const auto& t:tris)if(usable(t)) {d3[t.chart]+=t.area3;duv[t.chart]+=std::abs(t.area2);}
    struct Measured {std::array<MeshIndex,3> sectors{};double weight=0,area=0,shape=0;MeshIndex chart=0;};
    std::vector<Measured> measured;measured.reserve(tris.size());
    std::vector<std::uint8_t> chart_seen(ct.n_charts,0);
    for(const auto& t:tris) {
        if(!usable(t))continue;
        const double density=d3[t.chart]>0?duv[t.chart]/d3[t.chart]:1.0;
        const double ratio=t.singular_values[0]*t.singular_values[1]/std::max(density,1e-300);
        const double la=std::log(std::max(ratio,1e-300)),ln=std::log(t.singular_values[0]/t.singular_values[1]);
        if(!std::isfinite(la) || !std::isfinite(ln))continue;
        // TriangleQuality::area3 is the doubled 3D area (the source's convention).
        Measured m;m.chart=t.chart;m.weight=0.5*t.area3;m.area=la;m.shape=ln;
        for(std::size_t k=0;k<3;++k)m.sectors[k]=ct.sector_of_corner[t.corners[k]];
        measured.push_back(m);chart_seen[t.chart]=1;
    }
    out.triangles=measured.size();
    for(MeshIndex ch=0;ch<ct.n_charts;++ch)if(chart_seen[ch])out.charts.push_back(ch);
    std::vector<double> xa(ns,0.0),xn(ns,0.0);
    for(const auto& m:measured) {
        out.chart_area[m.chart]+=m.weight;
        for(auto s:m.sectors) {const double w=m.weight/3.0;out.mass[s]+=w;xa[s]+=w*m.area;xn[s]+=w*m.shape;}
    }
    for(MeshIndex s=0;s<ns;++s)if(out.mass[s]>0) {out.log_area[s]=xa[s]/out.mass[s];out.log_anisotropy[s]=xn[s]/out.mass[s];}
    out.raw_log_area=out.log_area;out.raw_log_anisotropy=out.log_anisotropy;
    // Per chart: sum(mass * value) and its scale sum(mass * |value|).
    auto moments=[&](const std::vector<double>& value) {
        std::vector<double> sum(ct.n_charts,0.0),scale(ct.n_charts,0.0);
        for(MeshIndex s=0;s<ns;++s)if(out.mass[s]>0) {const auto ch=ct.sector_chart[s];sum[ch]+=out.mass[s]*value[s];scale[ch]+=out.mass[s]*std::abs(value[s]);}
        return std::pair(sum,scale);
    };
    const auto area_before=moments(out.log_area),shape_before=moments(out.log_anisotropy);
    // Conservative smoothing on the cut adjacency: sector -> triangle (mean of its three
    // sectors) -> sector (mass-weighted mean of its triangles).
    for(std::int64_t step=0;step<smoothing_steps;++step) {
        std::vector<double> na(ns,0.0),nn(ns,0.0);
        for(const auto& m:measured) {
            double ya=0,yn=0;for(auto s:m.sectors) {ya+=out.log_area[s];yn+=out.log_anisotropy[s];}
            ya/=3.0;yn/=3.0;
            for(auto s:m.sectors) {na[s]+=m.weight/3.0*ya;nn[s]+=m.weight/3.0*yn;}
        }
        for(MeshIndex s=0;s<ns;++s)if(out.mass[s]>0) {out.log_area[s]=na[s]/out.mass[s];out.log_anisotropy[s]=nn[s]/out.mass[s];}
    }
    const auto area_after=moments(out.log_area),shape_after=moments(out.log_anisotropy);
    for(MeshIndex s=0;s<ns;++s)out.chart_mass[ct.sector_chart[s]]+=out.mass[s];
    for(auto ch:out.charts) {
        if(out.chart_area[ch]>0)out.mass_error=std::max(out.mass_error,std::abs(out.chart_mass[ch]-out.chart_area[ch])/out.chart_area[ch]);
        const double sa=std::max(area_before.second[ch],1e-300),sn=std::max(shape_before.second[ch],1e-300);
        out.moment_error=std::max({out.moment_error,std::abs(area_after.first[ch]-area_before.first[ch])/sa,
                                   std::abs(shape_after.first[ch]-shape_before.first[ch])/sn});
    }
    return out;
}
}
namespace cleanseams {
// ------------------------------------------------------------------ Phase B (T39)
nlohmann::json SeamGeometry::summary() const {
    return {{"version",version},{"area",area},{"seam_length",seam_length},{"normalized_length",normalized_length},{"spacing",spacing},
            {"centroid",centroid},{"seam_edges",seam_edges},{"charts",charts},{"disks",disks},{"boundary_loops",boundary_loops},
            {"samples",samples.size()}};
}
SeamGeometry seam_geometry(const CutSurface& cs,std::span<const std::uint8_t> exempt,double spacing) {
    const auto& m=cs.mesh();const auto& t=m.tables();const auto& g=m.geometry();const auto& ct=cs.tables();const auto& p=m.input().positions;
    const auto nf=g.face_area.size(),ne=t.edge_vertices.size();
    require(exempt.empty() || exempt.size()==nf,"metrics.seam_geometry.exempt","one flag per face required");
    require(std::isfinite(spacing) && spacing>0,"metrics.seam_geometry.spacing","spacing must be finite and positive");
    SeamGeometry out;out.spacing=spacing;
    // Normalization of the Section 20 comparisons: the area-weighted centroid and sqrt(area)
    // of the measured faces, so rigid motions and uniform scales leave the samples unchanged.
    Vec3 weighted{0,0,0};
    for(std::size_t f=0;f<nf;++f) {
        if((!exempt.empty() && exempt[f]) || !(g.face_area[f]>0))continue;
        out.area+=g.face_area[f];for(std::size_t k=0;k<3;++k)weighted[k]+=g.face_area[f]*g.face_center[f][k];
    }
    const double root=out.area>0?std::sqrt(out.area):1.0;
    if(out.area>0)for(std::size_t k=0;k<3;++k)out.centroid[k]=weighted[k]/out.area;
    for(std::size_t e=0;e<ne;++e) {
        if(!ct.cut[e] || t.edge_boundary[e])continue;
        ++out.seam_edges;out.seam_length+=g.edge_length[e];
        const auto& a=p[static_cast<std::size_t>(t.edge_vertices[e][0])];const auto& b=p[static_cast<std::size_t>(t.edge_vertices[e][1])];
        // Uniform samples along the edge (endpoints included), at most `spacing` apart in
        // normalized units: different tessellations of one seam give comparable point sets.
        const auto n=static_cast<std::size_t>(std::max(1.0,std::ceil(g.edge_length[e]/(spacing*root))));
        for(std::size_t k=0;k<=n;++k) {
            const double s=static_cast<double>(k)/static_cast<double>(n);Vec3 q{};
            for(std::size_t i=0;i<3;++i)q[i]=(a[i]+(b[i]-a[i])*s-out.centroid[i])/root;
            out.samples.push_back(q);
        }
    }
    out.normalized_length=out.seam_length/root;
    out.charts=ct.n_charts;for(auto d:ct.chart_is_disk)out.disks+=d!=0;
    for(auto n:ct.chart_n_loops)out.boundary_loops+=static_cast<std::uint64_t>(n);
    return out;
}
double seam_hausdorff(const SeamGeometry& a,const SeamGeometry& b) {
    if(a.samples.empty() && b.samples.empty())return 0.0;
    if(a.samples.empty() || b.samples.empty())return std::numeric_limits<double>::infinity();
    auto directed=[](const std::vector<Vec3>& from,const std::vector<Vec3>& to) {
        double worst=0;
        for(const auto& x:from) {
            double best=std::numeric_limits<double>::infinity();
            for(const auto& y:to) {const double d=std::hypot(x[0]-y[0],x[1]-y[1],x[2]-y[2]);if(d<best)best=d;}
            worst=std::max(worst,best);
        }
        return worst;
    };
    return std::max(directed(a.samples,b.samples),directed(b.samples,a.samples));
}
}
