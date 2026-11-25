// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/fields.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/control.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <functional>
#include <optional>
#include <numbers>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;using detail::cross;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
void storage(bool yes,const char* field) {if(!yes)throw ContractError(ErrorCode::memory_limit,field,"exposure workspace exceeds memory budget");}
double median_positive_edge(const MeshSnapshot& m) {
    std::vector<double> positive;for(double l:m.geometry().edge_length)if(l>0)positive.push_back(l);
    if(positive.empty())return 1.0;
    std::sort(positive.begin(),positive.end());const auto n=positive.size();
    return n%2?positive[n/2]:(positive[n/2-1]+positive[n/2])/2.0;
}
// Source _weighted_fraction: float32 weights summed pairwise in float64 per row.
std::vector<double> weighted_fraction(const std::vector<std::uint8_t>& visible,const std::vector<float>& weight,std::size_t rows,std::size_t k) {
    std::vector<double> out(rows),w(k),hit(k);
    for(std::size_t r=0;r<rows;++r) {
        for(std::size_t d=0;d<k;++d) {w[d]=weight[r*k+d];hit[d]=visible[r*k+d]?static_cast<double>(weight[r*k+d]):0.0;}
        const double total=detail::numpy_sum(w),seen=detail::numpy_sum(hit);
        out[r]=total>0?seen/std::max(total,1e-300):0.0;
    }
    return out;
}
}
std::vector<Vec3> fibonacci_directions(std::uint64_t count) {
    std::vector<Vec3> out;out.reserve(static_cast<std::size_t>(count));const double k=static_cast<double>(count);
    const double golden=std::numbers::pi*(1.0+std::sqrt(5.0));
    for(std::uint64_t n=0;n<count;++n) {
        const double i=static_cast<double>(n)+0.5,z=1.0-2.0*i/k,r=std::sqrt(std::max(1.0-z*z,0.0)),theta=golden*i;
        out.push_back({r*std::cos(theta),r*std::sin(theta),z});
    }
    return out;
}
ScreenFrames screen_frames(std::span<const Vec3> directions) {
    ScreenFrames out;
    for(const auto& d:directions) {
        std::size_t axis=0;for(std::size_t a=1;a<3;++a)if(std::abs(d[a])<std::abs(d[axis]))axis=a;
        Vec3 unit{0,0,0};unit[axis]=1;
        auto u=cross(d,unit);u=detail::divide(u,norm(u));out.u.push_back(u);out.w.push_back(cross(d,u));
    }
    return out;
}
namespace {
// The fixed-triangulation lattice of the source's surface_samples, generated on demand.
// Every triangle gets `steps` subdivisions of its longest side at the sample spacing;
// `other` is the face across each side that is a mesh edge (the triangle's face otherwise).
struct Lattice {
    std::vector<std::int64_t> steps;
    std::vector<std::array<TransferIndex,3>> other;
    std::uint64_t total=0; // face centres plus every lattice point
};
Lattice lattice(const MeshSnapshot& m,double spacing,std::uint64_t per_sample,std::uint64_t budget,const char* field) {
    require(std::isfinite(spacing),"exposure.spacing","sample spacing must be finite");
    const auto& t=m.tables();const auto& tri=m.triangulation();
    const auto nf=m.input().face_offsets.size()-1,nt=tri.corners.size();
    Lattice out;out.total=nf;
    if(!nt)return out;
    auto position=[&](MeshIndex corner) {return m.input().positions[static_cast<std::size_t>(m.input().corner_vertices[corner])];};
    out.steps.resize(nt);
    for(std::size_t i=0;i<nt;++i) {
        double longest=-INFINITY;
        for(std::size_t k=0;k<3;++k)longest=std::max(longest,norm(sub(position(tri.corners[i][k]),position(tri.corners[i][(k+1)%3]))));
        const double n=spacing>0?std::ceil(longest/spacing):1.0;
        storage(std::isfinite(n) && n<=1073741824.0,field);
        out.steps[i]=std::max<std::int64_t>(static_cast<std::int64_t>(n),1);
        const auto s=static_cast<std::uint64_t>(out.steps[i]);out.total+=(s+1)*(s+2)/2;
        storage(out.total<=budget/per_sample,field);
    }
    out.other.resize(nt);
    const std::array<std::array<std::size_t,2>,3> sides{{{0,1},{1,2},{2,0}}};
    for(std::size_t i=0;i<nt;++i) {
        const auto face=tri.faces[i];out.other[i]={face,face,face};
        for(std::size_t s=0;s<3;++s) {
            const auto ca=tri.corners[i][sides[s][0]],cb=tri.corners[i][sides[s][1]];
            if(m.corner_next(ca)!=cb && m.corner_next(cb)!=ca)continue;
            const auto c=m.corner_next(ca)==cb?ca:cb;const auto fs=t.edge_faces[t.corner_edge[c]];
            if(fs.size()==2)out.other[i][s]=fs[0]==face?fs[1]:fs[0];
        }
    }
    return out;
}
// The lattice points of triangle `i` in the source's meshgrid(ij) order with i + j <= k,
// as f(point, face_a, face_b, vertex): the side a point lies on is j == 0 -> (0,1),
// i + j == k -> (1,2), i == 0 -> (2,0) (the last assignment wins, as in the source).
template<class F> void lattice_points(const MeshSnapshot& m,const Lattice& l,std::size_t triangle,F&& f) {
    const auto& tri=m.triangulation();const auto k=l.steps[triangle];const double kk=static_cast<double>(k);
    const auto& corners=tri.corners[triangle];const auto face=static_cast<TransferIndex>(tri.faces[triangle]);
    auto position=[&](MeshIndex corner) {return m.input().positions[static_cast<std::size_t>(m.input().corner_vertices[corner])];};
    const auto p0=position(corners[0]),p1=position(corners[1]),p2=position(corners[2]);
    for(std::int64_t i=0;i<=k;++i)for(std::int64_t j=0;i+j<=k;++j) {
        std::int64_t side=-1,corner=-1;
        if(i==0)side=2;
        if(i+j==k)side=1;
        if(j==0)side=0;
        if(i==0 && j==0)corner=0;
        if(i==k && j==0)corner=1;
        if(i==0 && j==k)corner=2;
        const double bi=static_cast<double>(i)/kk,bj=static_cast<double>(j)/kk;Vec3 p;
        for(std::size_t d=0;d<3;++d)p[d]=(p0[d]+bi*(p1[d]-p0[d]))+bj*(p2[d]-p0[d]);
        f(p,face,side>=0?l.other[triangle][static_cast<std::size_t>(side)]:face,
          corner>=0?m.input().corner_vertices[corners[static_cast<std::size_t>(corner)]]:TransferIndex{-1});
    }
}
// Minimum depths per screen cell over all samples and over back-facing samples; an
// empty cell is +inf in both. Dense arrays when the grid is small next to the samples,
// otherwise only the occupied cells (open addressing): the answers are the same, and a
// far-apart part or one huge triangle no longer allocates the whole empty screen.
class Depths {
public:
    Depths(std::uint64_t cells,std::uint64_t samples,std::uint64_t budget):budget_(budget) {
        dense_=cells<=std::max<std::uint64_t>(std::uint64_t{1}<<20,4*samples) && cells<=budget/32;
        if(dense_) {all_.assign(static_cast<std::size_t>(cells),INFINITY);back_.assign(static_cast<std::size_t>(cells),INFINITY);}
        else resize(std::size_t{1}<<12);
    }
    void lower(std::uint64_t key,double z,bool back) {
        if(dense_) {
            auto& a=all_[static_cast<std::size_t>(key)];a=std::min(a,z);
            if(back) {auto& b=back_[static_cast<std::size_t>(key)];b=std::min(b,z);}
            return;
        }
        auto slot=find(key);
        if(keys_[slot]==empty) {
            if(2*(used_+1)>keys_.size()) {resize(2*keys_.size());slot=find(key);}
            keys_[slot]=key;depth_[slot]={INFINITY,INFINITY};++used_;
        }
        auto& d=depth_[slot];d[0]=std::min(d[0],z);if(back)d[1]=std::min(d[1],z);
    }
    std::array<double,2> at(std::uint64_t key) const {
        if(dense_)return {all_[static_cast<std::size_t>(key)],back_[static_cast<std::size_t>(key)]};
        const auto slot=find(key);
        return keys_[slot]==empty?std::array<double,2>{INFINITY,INFINITY}:depth_[slot];
    }
private:
    static constexpr std::uint64_t empty=~std::uint64_t{0};
    bool dense_=true;
    std::vector<double> all_,back_;
    std::vector<std::uint64_t> keys_;
    std::vector<std::array<double,2>> depth_;
    std::uint64_t used_=0,budget_=0;
    std::size_t find(std::uint64_t key) const {
        const auto mask=keys_.size()-1;
        auto slot=static_cast<std::size_t>((key*0x9E3779B97F4A7C15ULL)>>17)&mask;
        while(keys_[slot]!=empty && keys_[slot]!=key)slot=(slot+1)&mask;
        return slot;
    }
    void resize(std::size_t capacity) {
        // Old and new tables are live while rehashing.
        const auto bytes=static_cast<std::uint64_t>(capacity+keys_.size())*(sizeof(std::uint64_t)+2*sizeof(double));
        storage(capacity<=(std::numeric_limits<std::size_t>::max()>>6) && bytes<=budget_/2,"exposure.zbuffer");
        auto keys=std::move(keys_);auto depth=std::move(depth_);
        keys_.assign(capacity,empty);depth_.assign(capacity,{INFINITY,INFINITY});
        for(std::size_t i=0;i<keys.size();++i)if(keys[i]!=empty) {const auto slot=find(keys[i]);keys_[slot]=keys[i];depth_[slot]=depth[i];}
    }
};
}
namespace {
// Source order: face centres, then triangles grouped by ascending step count.
SurfaceSamples materialize(const MeshSnapshot& m,const Lattice& l) {
    const auto& g=m.geometry();const auto nf=m.input().face_offsets.size()-1;
    SurfaceSamples out;
    out.points.reserve(static_cast<std::size_t>(l.total));out.face_a.reserve(static_cast<std::size_t>(l.total));
    out.face_b.reserve(static_cast<std::size_t>(l.total));out.vertex.reserve(static_cast<std::size_t>(l.total));
    for(std::size_t f=0;f<nf;++f) {out.points.push_back(g.face_center[f]);out.face_a.push_back(TransferIndex(f));out.face_b.push_back(TransferIndex(f));out.vertex.push_back(-1);}
    // Source order: triangles grouped by ascending step count, then triangle order.
    std::vector<std::int64_t> sizes(l.steps);std::sort(sizes.begin(),sizes.end());sizes.erase(std::unique(sizes.begin(),sizes.end()),sizes.end());
    for(auto k:sizes)for(std::size_t i=0;i<l.steps.size();++i) {
        if(l.steps[i]!=k)continue;
        lattice_points(m,l,i,[&](const Vec3& p,TransferIndex a,TransferIndex b,TransferIndex v) {
            out.points.push_back(p);out.face_a.push_back(a);out.face_b.push_back(b);out.vertex.push_back(v);
        });
    }
    return out;
}
}
SurfaceSamples surface_samples(const MeshSnapshot& m,double spacing,std::uint64_t budget) {
    // Materialized samples: 48 bytes each plus vector growth.
    return materialize(m,lattice(m,spacing,64,budget,"exposure.samples"));
}
Visibility zbuffer_visibility(const MeshSnapshot& m,std::span<const Vec3> directions,const VisibilityOptions& o,const SurfaceSamples* given) {
    const auto& g=m.geometry();const auto nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1,k=directions.size();
    for(const auto& d:directions)require(std::isfinite(d[0]) && std::isfinite(d[1]) && std::isfinite(d[2]),"exposure.directions","nonfinite direction");
    require(std::isfinite(o.cell_factor) && std::isfinite(o.base) && std::isfinite(o.tan_max) && std::isfinite(o.spacing),"exposure.options","nonfinite visibility option");
    storage(k==0 || (nv+nf)<=o.memory_budget/(8*k),"exposure.visibility");
    Visibility out;out.directions.assign(directions.begin(),directions.end());
    out.vertex.assign(nv*k,0);out.face.assign(nf*k,0);out.vertex_weight.assign(nv*k,0.f);out.face_weight.assign(nf*k,0.f);
    if(nf==0 || k==0)return out;
    const double median=median_positive_edge(m),cell=std::max(o.cell_factor*median,1e-12);
    // Gate A F2: the lattice is bounded by the budget at 32 bytes a sample (the depth, x, y
    // and cell key the source's vectorized z-buffer pass holds per sample and direction).
    // Within an eighth of the budget at 48 bytes a sample it is stored once for every
    // direction; a larger one is regenerated triangle by triangle in each pass, never
    // stored. The points are the same, and their order does not change a cell's minimum.
    std::optional<Lattice> stream;SurfaceSamples own;const SurfaceSamples* samples=given;
    if(!given) {
        stream=lattice(m,o.spacing*median,32,o.memory_budget,"exposure.samples");
        if(stream->total<=o.memory_budget/(8*48)) {own=materialize(m,*stream);samples=&own;stream.reset();}
    }
    const auto ns=samples?static_cast<std::uint64_t>(samples->points.size()):stream->total;
    auto for_each_sample=[&](auto&& f) {
        if(samples) {for(std::size_t i=0;i<samples->points.size();++i)f(samples->points[i],samples->face_a[i],samples->face_b[i],samples->vertex[i]);return;}
        for(std::size_t f_=0;f_<nf;++f_)f(g.face_center[f_],TransferIndex(f_),TransferIndex(f_),TransferIndex{-1});
        for(std::size_t i=0;i<stream->steps.size();++i)lattice_points(m,*stream,i,f);
    };
    const auto frames=screen_frames(directions);
    const auto& cv=m.input().corner_vertices;const auto& corner_face=m.tables().corner_face;
    // GB performance: the directions are independent passes, each writing only its own
    // column d of the tables, so with a worker pool they run in parallel, each with its own
    // scratch and z-buffer (inline in direction order without one). The values are those of
    // the serial loop; a failing pass reports the lowest failing direction, as it did. Passes
    // run in parallel only while the simultaneous scratch (the dense z-buffer bound and the
    // per-vertex/face arrays of every thread) stays within a quarter of the budget.
    const auto* pool=active_workers();
    const std::uint64_t pass_bytes=16*std::max<std::uint64_t>(std::uint64_t{1}<<20,4*ns)+13*8*static_cast<std::uint64_t>(nv+nf);
    const bool parallel=pool && !in_worker_task() && pass_bytes<=o.memory_budget/4/std::max<std::uint64_t>(pool->threads(),1);
    auto each_direction=[&](std::size_t n,const std::function<void(std::size_t)>& pass) {
        if(parallel)parallel_for(n,pass);
        else for(std::size_t d=0;d<n;++d)pass(d);
    };
    each_direction(k,[&](std::size_t d) {
        std::vector<double> vx(nv),vy(nv),vz(nv),cx(nf),cy(nf),cz(nf),cosf(nf),tanf(nf),cosv(nv),tanv(nv);
        std::vector<std::uint8_t> front_f(nf),back_f(nf),front_v(nv);
        const auto& dir=directions[d];const auto& u=frames.u[d];const auto& w=frames.w[d];
        for(std::size_t i=0;i<nv;++i) {const auto& p=m.input().positions[i];vx[i]=dot(p,u);vy[i]=dot(p,w);vz[i]=-dot(p,dir);}
        for(std::size_t f=0;f<nf;++f) {
            const auto& c=g.face_center[f];cx[f]=dot(c,u);cy[f]=dot(c,w);cz[f]=-dot(c,dir);
            const double nd=dot(g.face_normal[f],dir);cosf[f]=std::max(nd,0.0);
            front_f[f]=cosf[f]>1e-9;back_f[f]=nd< -1e-9;
            const double tangent=std::sqrt(std::max(1.0-cosf[f]*cosf[f],0.0))/std::max(cosf[f],1e-12);
            tanf[f]=front_f[f]?std::min(tangent,o.tan_max):0.0;
        }
        std::fill(cosv.begin(),cosv.end(),0.0);std::fill(tanv.begin(),tanv.end(),0.0);
        for(std::size_t c=0;c<cv.size();++c) {
            const auto v=static_cast<std::size_t>(cv[c]);const auto f=corner_face[c];
            cosv[v]=std::max(cosv[v],cosf[f]);tanv[v]=std::max(tanv[v],tanf[f]);
        }
        for(std::size_t v=0;v<nv;++v)front_v[v]=cosv[v]>1e-9;
        // The screen extent of the samples (every sample, as the source's ox.min/ox.max).
        double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
        for_each_sample([&](const Vec3& p,TransferIndex,TransferIndex,TransferIndex) {
            const double x=dot(p,u),y=dot(p,w);
            xmin=std::min(xmin,x);xmax=std::max(xmax,x);ymin=std::min(ymin,y);ymax=std::max(ymax,y);
        });
        const double xspan=(xmax-xmin)/cell,yspan=(ymax-ymin)/cell;
        // Cell keys ix * nrow + iy must fit: beyond 2^62 cells the source's dense buffers
        // (16 bytes a cell) cannot exist on any machine either.
        storage(std::isfinite(xspan) && std::isfinite(yspan) && (xspan+2.0)*(yspan+2.0)<0x1p62,"exposure.zbuffer");
        const auto ncol=static_cast<std::int64_t>(xspan)+2,nrow=static_cast<std::int64_t>(yspan)+2;
        // Source np.clip(((x - xmin) / cell).astype(int64), 0, n - 1). Points outside the
        // sampled extent (unused vertices far away) gave quotients beyond int64, and
        // converting those is undefined behavior (UBSan, T25): clamp before truncating.
        // Same cells as the saturating conversion of the recorded arm64 host; NaN -> 0.
        auto cell_index=[](double q,std::int64_t count) {
            if(!(q>=0))return std::int64_t{0};
            return q>=static_cast<double>(count-1)?count-1:static_cast<std::int64_t>(q);
        };
        auto key=[&](double x,double y) {
            const auto ix=cell_index((x-xmin)/cell,ncol),iy=cell_index((y-ymin)/cell,nrow);
            return static_cast<std::uint64_t>(ix)*static_cast<std::uint64_t>(nrow)+static_cast<std::uint64_t>(iy);
        };
        Depths depths(static_cast<std::uint64_t>(ncol)*static_cast<std::uint64_t>(nrow),ns,o.memory_budget);
        for_each_sample([&](const Vec3& p,TransferIndex face_a,TransferIndex face_b,TransferIndex vertex) {
            const double x=dot(p,u),y=dot(p,w),z=-dot(p,dir);
            const bool vertex_front=vertex>=0 && front_v[static_cast<std::size_t>(vertex)];
            depths.lower(key(x,y),z,back_f[static_cast<std::size_t>(face_a)] && back_f[static_cast<std::size_t>(face_b)] && !vertex_front);
        });
        const double noise=o.base*cell;
        for(std::size_t v=0;v<nv;++v) {
            const auto z=depths.at(key(vx[v],vy[v]));const double tol=o.base*cell+cell*tanv[v];
            out.vertex[v*k+d]=front_v[v] && z[0]>=vz[v]-tol && z[1]>=vz[v]-noise;
            out.vertex_weight[v*k+d]=static_cast<float>(cosv[v]);
        }
        for(std::size_t f=0;f<nf;++f) {
            const auto z=depths.at(key(cx[f],cy[f]));const double tol=o.base*cell+cell*tanf[f];
            out.face[f*k+d]=front_f[f] && z[0]>=cz[f]-tol && z[1]>=cz[f]-noise;
            out.face_weight[f*k+d]=static_cast<float>(cosf[f]);
        }
    });
    return out;
}
ExposureField exposure_field(const MeshSnapshot& m,std::uint64_t count,std::optional<std::span<const Vec3>> directions,const VisibilityOptions& o) {
    ExposureField out;
    if(directions)out.directions.assign(directions->begin(),directions->end());else out.directions=fibonacci_directions(count);
    auto visibility=std::make_shared<Visibility>(zbuffer_visibility(m,out.directions,o));
    const auto k=out.directions.size(),nv=m.input().positions.size(),nf=m.input().face_offsets.size()-1;
    out.vertex=weighted_fraction(visibility->vertex,visibility->vertex_weight,nv,k);
    out.face=weighted_fraction(visibility->face,visibility->face_weight,nf,k);
    const auto& edges=m.tables().edge_vertices;out.edge.resize(edges.size());
    for(std::size_t e=0;e<edges.size();++e)out.edge[e]=0.5*(out.vertex[edges[e][0]]+out.vertex[edges[e][1]]);
    out.visibility=std::move(visibility);out.geometry_key=m.geometry_key();
    return out;
}
std::shared_ptr<const ExposureField> ExposureCache::get(const MeshSnapshot& m,std::uint64_t count,const VisibilityOptions& o) {
    const Key key{m.geometry_key(),count,o};
    {
        std::lock_guard lock(mutex_);
        if(auto found=entries_.find(key);found!=entries_.end()) {++hits_;return found->second;}
    }
    auto field=std::make_shared<const ExposureField>(exposure_field(m,count,std::nullopt,o));
    std::lock_guard lock(mutex_);++misses_;
    if(!entries_.contains(key)) {
        if(max_entries_==0)return field;
        while(entries_.size()>=max_entries_) {entries_.erase(order_.front());order_.erase(order_.begin());}
        entries_.emplace(key,field);order_.push_back(key);
    }
    return entries_.at(key);
}
std::size_t ExposureCache::size() const {std::lock_guard lock(mutex_);return entries_.size();}
double mesh_exposure(const MeshSnapshot& m,const ExposureField& field) {
    const auto& lengths=m.geometry().edge_length;require(field.edge.size()==lengths.size(),"exposure.field","field does not match the mesh");
    std::vector<double> weighted(lengths.size());for(std::size_t e=0;e<lengths.size();++e)weighted[e]=field.edge[e]*lengths[e];
    const double total=detail::numpy_sum(lengths);
    return total>0?detail::numpy_sum(weighted)/total:std::numeric_limits<double>::quiet_NaN();
}
SeamExposure seam_exposure(const MeshSnapshot& m,std::span<const std::uint8_t> seams,const ExposureField& field,double hidden_below) {
    const auto& lengths=m.geometry().edge_length;const auto nan=std::numeric_limits<double>::quiet_NaN();
    require(seams.size()==lengths.size() && field.edge.size()==lengths.size(),"exposure.seams","seam mask or field does not match the mesh");
    std::vector<double> seam_lengths,seam_weighted,hidden,weighted(lengths.size());
    for(std::size_t e=0;e<lengths.size();++e) {
        weighted[e]=field.edge[e]*lengths[e];
        if(!seams[e])continue;
        seam_lengths.push_back(lengths[e]);seam_weighted.push_back(field.edge[e]*lengths[e]);
        if(field.edge[e]<hidden_below)hidden.push_back(lengths[e]);
    }
    const double ls=detail::numpy_sum(seam_lengths),lm=detail::numpy_sum(lengths);
    SeamExposure out;
    out.mesh_exposure=lm>0?detail::numpy_sum(weighted)/lm:nan;
    out.seam_exposure=ls>0?detail::numpy_sum(seam_weighted)/ls:nan;
    const double share=ls>0?detail::numpy_sum(hidden)/ls:nan;
    out.exposure_ratio=ls>0 && out.mesh_exposure>0?out.seam_exposure/out.mesh_exposure:nan;
    out.hidden_seam_pct=100.0*share;
    return out;
}
Vec3 view_direction(const Vec3& vector) {
    const double length=norm(vector);
    if(!std::isfinite(vector[0]) || !std::isfinite(vector[1]) || !std::isfinite(vector[2]) || !(length>0))
        throw ContractError(ErrorCode::invalid_option,"view_direction","view direction must be a signed axis or a finite nonzero vector");
    return detail::divide(vector,length);
}
Vec3 view_direction(std::string_view axis) {
    if(axis.size()!=2 || (axis[0]!='+' && axis[0]!='-') || (axis[1]!='X' && axis[1]!='Y' && axis[1]!='Z'))
        throw ContractError(ErrorCode::invalid_option,"view_direction","view direction must be a signed axis or a finite nonzero vector");
    Vec3 v{0,0,0};v[static_cast<std::size_t>(axis[1]-'X')]=axis[0]=='+'?1.0:-1.0;return view_direction(v);
}
double seam_on_front_pct(const MeshSnapshot& m,std::span<const std::uint8_t> seams,std::optional<Vec3> front) {
    if(!front)return std::numeric_limits<double>::quiet_NaN();
    const auto d=view_direction(*front);const auto& t=m.tables();const auto& g=m.geometry();
    require(seams.size()==t.edge_vertices.size(),"exposure.seams","seam mask does not match the mesh");
    std::vector<double> lengths,facing;
    for(std::size_t e=0;e<seams.size();++e) {
        if(!seams[e] || t.edge_boundary[e])continue;
        const auto& a=g.vert_normal[t.edge_vertices[e][0]];const auto& b=g.vert_normal[t.edge_vertices[e][1]];
        const Vec3 normal{(a[0]+b[0])/2.0,(a[1]+b[1])/2.0,(a[2]+b[2])/2.0};
        lengths.push_back(g.edge_length[e]);if(dot(normal,d)>1e-9)facing.push_back(g.edge_length[e]);
    }
    return 100*detail::numpy_sum(facing)/std::max(detail::numpy_sum(lengths),1e-300);
}
SeamCostField seam_cost(const MeshSnapshot& m,const ExposureField& field,const SeamCostOptions& o) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    require(field.edge.size()==ne,"exposure.field","field does not match the mesh");
    for(double x:{o.feature_angle,o.valley_angle,o.valley_factor,o.hard_factor,o.hide_strength})
        if(!std::isfinite(x))throw ContractError(ErrorCode::invalid_option,"seam_cost","nonfinite seam cost option");
    SeamCostField out;out.exposure=field.edge;out.valley.resize(ne);out.hard.assign(ne,0);out.cost=field.edge;
    const double valley=o.valley_angle*(std::numbers::pi/180.0),hard=o.feature_angle*(std::numbers::pi/180.0);
    for(std::size_t e=0;e<ne;++e) {
        out.valley[e]=g.edge_dihedral[e]<=-valley;
        if(o.feature_angle>0)out.hard[e]=std::abs(g.edge_dihedral[e])>=hard;
        if(out.valley[e])out.cost[e]*=o.valley_factor;
        if(out.hard[e])out.cost[e]*=o.hard_factor;
    }
    if(o.hide_dir) {
        auto d=view_direction(*o.hide_dir);const double n=norm(d);
        if(n>0) {
            d=detail::divide(d,n);out.facing.resize(ne);
            for(std::size_t e=0;e<ne;++e) {
                const auto& a=g.vert_normal[t.edge_vertices[e][0]];const auto& b=g.vert_normal[t.edge_vertices[e][1]];
                const Vec3 normal{0.5*(a[0]+b[0]),0.5*(a[1]+b[1]),0.5*(a[2]+b[2])};
                out.facing[e]=std::clamp(dot(normal,d),-1.0,1.0);
                out.cost[e]=out.cost[e]*(1.0-o.hide_strength*0.5*(1.0-out.facing[e]));
            }
        }
    }
    for(auto& c:out.cost)c=std::clamp(c,0.0,1.0);
    return out;
}
}
