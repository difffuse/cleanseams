// SPDX-License-Identifier: GPL-3.0-or-later
// Source core/symmetry.py.
#include "cleanseams/symmetry.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;
constexpr double inf=std::numeric_limits<double>::infinity();
constexpr double nan=std::numeric_limits<double>::quiet_NaN();
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::size_t face_count(const MeshSnapshot& m) {return m.input().face_offsets.size()-1;}
Vec3 position(const MeshSnapshot& m,std::size_t v) {return m.input().positions[v];}
MeshIndex corner_vertex(const MeshSnapshot& m,std::size_t c) {return static_cast<MeshIndex>(m.input().corner_vertices[c]);}
// numpy.einsum("ij,ij->i") on 3-vectors (detail::numpy_einsum_dot3).
double edot(const Vec3& a,const Vec3& b) {return detail::numpy_einsum_dot3(a,b);}
// The source's (X * w).T @ X and d.T @ d products run in the reference build's BLAS,
// whose accumulation order is not reproduced; a fused chain in row order is used.
std::array<std::array<double,3>,3> gram(std::span<const Vec3> x,std::span<const double> w) {
    std::array<std::array<double,3>,3> out{};
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
        double s=0;for(std::size_t k=0;k<x.size();++k)s=std::fma(w.empty()?x[k][i]:x[k][i]*w[k],x[k][j],s);
        out[i][j]=s;
    }
    return out;
}
// numpy.median of the positive edge lengths, or 1.
double positive_edge_median(const MeshSnapshot& m) {
    std::vector<double> positive;for(double l:m.geometry().edge_length)if(l>0)positive.push_back(l);
    return positive.empty()?1.0:detail::numpy_median(std::move(positive));
}
}
// ------------------------------------------------------------------ grid hash
NearestPoints nearest_points(std::span<const Vec3> points,std::span<const Vec3> queries,double radius) {
    NearestPoints out;out.index.assign(queries.size(),-1);out.distance.assign(queries.size(),inf);
    if(points.empty() || queries.empty() || !(radius>0))return out;
    const double cell=radius;
    Vec3 lo{inf,inf,inf};for(const auto& p:points)for(std::size_t d=0;d<3;++d)lo[d]=std::min(lo[d],p[d]);
    for(std::size_t d=0;d<3;++d)lo[d]-=cell;
    using Key=std::array<std::int64_t,3>;
    // Source np.floor(...).astype(np.int64): the arm64 conversion saturates (NaN -> 0);
    // converting an out-of-range quotient in C++ is undefined behaviour, so clamp first
    // (Gate A F32, the class of the T25 exposure fix). Saturated queries fall outside.
    auto index=[](double q) {
        constexpr double limit=9223372036854775808.0; // 2^63
        if(std::isnan(q))return std::int64_t{0};
        if(q>=limit)return std::numeric_limits<std::int64_t>::max();
        if(q<-limit)return std::numeric_limits<std::int64_t>::min();
        return static_cast<std::int64_t>(q);
    };
    auto grid=[&](const Vec3& p) {Key k;for(std::size_t d=0;d<3;++d)k[d]=index(std::floor((p[d]-lo[d])/cell));return k;};
    // Cells sorted by (x, y, z) then point index: the source's stable argsort of the
    // linear key (which never aliases: two empty guard layers surround the points).
    std::vector<std::pair<Key,std::size_t>> cells(points.size());Key dims{0,0,0};
    for(std::size_t i=0;i<points.size();++i) {cells[i]={grid(points[i]),i};for(std::size_t d=0;d<3;++d)dims[d]=std::max(dims[d],cells[i].first[d]);}
    for(auto& d:dims)d=d>std::numeric_limits<std::int64_t>::max()-3?std::numeric_limits<std::int64_t>::max():d+3;
    std::sort(cells.begin(),cells.end());
    for(std::size_t q=0;q<queries.size();++q) {
        const auto g=grid(queries[q]);bool inside=true;
        for(std::size_t d=0;d<3;++d)inside=inside && g[d]>=0 && g[d]<dims[d]-1;
        if(!inside)continue;
        double best=inf;TransferIndex winner=-1;bool unordered=false,any=false;
        for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz) {
            const Key k{g[0]+dx,g[1]+dy,g[2]+dz};
            auto it=std::lower_bound(cells.begin(),cells.end(),std::pair<Key,std::size_t>{k,0});
            for(;it!=cells.end() && it->first==k;++it) {
                const double d=norm(sub(queries[q],points[it->second]));any=true;
                if(std::isnan(d))unordered=true;
                // The last candidate at the minimum distance wins, as the source's
                // repeated fancy assignment does.
                if(d<=best) {best=d;winner=static_cast<TransferIndex>(it->second);}
            }
        }
        // A NaN distance makes the source's running minimum NaN: no candidate wins.
        if(!any || unordered || !(best<=radius))continue;
        out.index[q]=winner;out.distance[q]=best;
    }
    return out;
}
// --------------------------------------------------------- point to triangle
std::vector<Vec3> closest_point_on_triangles(std::span<const Vec3> p,std::span<const Vec3> a,std::span<const Vec3> b,std::span<const Vec3> c) {
    require(p.size()==a.size() && a.size()==b.size() && b.size()==c.size(),"symmetry.triangles","point and triangle counts differ");
    std::vector<Vec3> out(p.size());
    auto lerp=[](const Vec3& base,double t,const Vec3& dir) {return Vec3{base[0]+t*dir[0],base[1]+t*dir[1],base[2]+t*dir[2]};};
    for(std::size_t i=0;i<p.size();++i) {
        const auto ab=sub(b[i],a[i]),ac=sub(c[i],a[i]),ap=sub(p[i],a[i]);
        const double d1=edot(ab,ap),d2=edot(ac,ap);
        const auto bp=sub(p[i],b[i]);const double d3=edot(ab,bp),d4=edot(ac,bp);
        const auto cp=sub(p[i],c[i]);const double d5=edot(ab,cp),d6=edot(ac,cp);
        const double vc=d1*d4-d3*d2,vb=d5*d2-d1*d6,va=d3*d6-d5*d4;
        if(d1<=0 && d2<=0) {out[i]=a[i];continue;}
        if(d3>=0 && d4<=d3) {out[i]=b[i];continue;}
        if(d6>=0 && d5<=d6) {out[i]=c[i];continue;}
        if(vc<=0 && d1>=0 && d3<=0) {out[i]=lerp(a[i],d1/(d1-d3),ab);continue;}
        if(vb<=0 && d2>=0 && d6<=0) {out[i]=lerp(a[i],d2/(d2-d6),ac);continue;}
        if(va<=0 && d4-d3>=0 && d5-d6>=0) {out[i]=lerp(b[i],(d4-d3)/((d4-d3)+(d5-d6)),sub(c[i],b[i]));continue;}
        const double denom=va+vb+vc,v=vb/denom,w=vc/denom;
        const auto inner=lerp(lerp(a[i],v,ab),w,ac);
        out[i]=std::isfinite(inner[0]) && std::isfinite(inner[1]) && std::isfinite(inner[2])?inner:a[i];
    }
    return out;
}
SurfaceDistance surface_distance(const MeshSnapshot& m,std::span<const Vec3> queries,std::optional<double> radius_in) {
    const auto nv=m.input().positions.size();const auto nq=queries.size();
    const double radius=radius_in?*radius_in:2.0*positive_edge_median(m);
    std::vector<std::uint8_t> used(nv,0);for(auto v:m.input().corner_vertices)used[static_cast<std::size_t>(v)]=1;
    std::vector<MeshIndex> used_ids;std::vector<Vec3> used_points;
    for(std::size_t v=0;v<nv;++v)if(used[v]) {used_ids.push_back(static_cast<MeshIndex>(v));used_points.push_back(position(m,v));}
    const auto near=nearest_points(used_points,queries,radius);
    SurfaceDistance out;out.distance.assign(nq,inf);out.vertex.assign(nq,-1);out.foot.assign(queries.begin(),queries.end());
    for(std::size_t q=0;q<nq;++q)if(near.index[q]>=0)out.vertex[q]=used_ids[static_cast<std::size_t>(near.index[q])];
    const auto& tri=m.triangulation().corners;
    if(tri.empty()) {
        for(std::size_t q=0;q<nq;++q)if(out.vertex[q]>=0) {out.distance[q]=near.distance[q];out.foot[q]=position(m,static_cast<std::size_t>(out.vertex[q]));}
        return out;
    }
    // vertex -> triangles, ascending triangle id, a triangle repeated per incidence
    std::vector<std::vector<std::size_t>> incident(nv);
    for(std::size_t t=0;t<tri.size();++t)for(auto c:tri[t])incident[corner_vertex(m,c)].push_back(t);
    for(std::size_t q=0;q<nq;++q) {
        if(out.vertex[q]<0)continue;
        const auto& around=incident[static_cast<std::size_t>(out.vertex[q])];
        std::vector<Vec3> P(around.size(),queries[q]),A,B,C;
        for(auto t:around) {A.push_back(position(m,corner_vertex(m,tri[t][0])));B.push_back(position(m,corner_vertex(m,tri[t][1])));C.push_back(position(m,corner_vertex(m,tri[t][2])));}
        const auto cp=closest_point_on_triangles(P,A,B,C);
        double best=inf;bool unordered=false;std::optional<std::size_t> winner;
        for(std::size_t k=0;k<cp.size();++k) {
            const double d=norm(sub(queries[q],cp[k]));if(std::isnan(d))unordered=true;
            if(d<=best) {best=d;winner=k;}
        }
        if(unordered) {best=nan;winner.reset();}
        if(winner)out.foot[q]=cp[*winner];
        out.distance[q]=best;
        if(!std::isfinite(best)) {out.distance[q]=near.distance[q];out.foot[q]=position(m,static_cast<std::size_t>(out.vertex[q]));}
    }
    return out;
}
// ------------------------------------------------------------------- planes
double MirrorPlane::signed_distance(const Vec3& p) const {return dot(p,normal)-offset;}
Vec3 MirrorPlane::reflect(const Vec3& p) const {
    const double twice=2.0*signed_distance(p);return {p[0]-twice*normal[0],p[1]-twice*normal[1],p[2]-twice*normal[2]};
}
nlohmann::json MirrorPlane::to_json() const {
    return {{"normal",{normal[0],normal[1],normal[2]}},{"offset",offset},{"score",score},{"far",far},{"source",source}};
}
double bbox_diagonal(const MeshSnapshot& m) {
    const auto& v=m.input().positions;if(v.empty())return 0.0;
    Vec3 lo=v[0],hi=v[0];for(const auto& p:v)for(std::size_t d=0;d<3;++d) {lo[d]=std::min(lo[d],p[d]);hi[d]=std::max(hi[d],p[d]);}
    return norm(sub(hi,lo));
}
namespace {
Vec3 unit(Vec3 n) {double l=norm(n);if(l==0)l=1.0;return {n[0]/l,n[1]/l,n[2]/l};}
}
detail::FaceCovariance detail::face_centre_covariance(const MeshSnapshot& m) {
    FaceCovariance out;const auto nf=face_count(m);
    if(nf==0)return out;
    const auto& g=m.geometry();
    std::vector<double> w(g.face_area.begin(),g.face_area.end());double total=detail::numpy_sum(w);
    if(total<=0) {w.assign(nf,1.0);total=static_cast<double>(nf);}
    auto& centroid=out.centroid;for(std::size_t f=0;f<nf;++f)for(std::size_t d=0;d<3;++d)centroid[d]+=g.face_center[f][d]*w[f];
    for(auto& x:centroid)x/=total;
    std::vector<Vec3> x(nf);for(std::size_t f=0;f<nf;++f)x[f]=sub(g.face_center[f],centroid);
    out.matrix=gram(x,w);for(auto& row:out.matrix)for(auto& value:row)value/=total;
    return out;
}
std::vector<MirrorPlane> candidate_planes(const MeshSnapshot& m) {
    std::vector<MirrorPlane> out;const auto nf=face_count(m);
    if(nf==0)return out;
    const auto [cov,centroid]=detail::face_centre_covariance(m);
    const auto eig=detail::lapack_syevd3(cov);
    for(std::size_t i=0;i<3;++i) {
        const auto n=unit({eig.vectors[0][i],eig.vectors[1][i],eig.vectors[2][i]});
        out.push_back({n,dot(n,centroid),0.0,"pca"+std::to_string(i),0.0});
    }
    const auto& v=m.input().positions;Vec3 lo=v[0],hi=v[0];
    for(const auto& p:v)for(std::size_t d=0;d<3;++d) {lo[d]=std::min(lo[d],p[d]);hi[d]=std::max(hi[d],p[d]);}
    for(std::size_t i=0;i<3;++i) {
        Vec3 n{0,0,0};n[i]=1.0;
        out.push_back({n,0.5*(hi[i]+lo[i]),0.0,std::string("bbox_")+"xyz"[i],0.0});
    }
    return out;
}
PlaneScore score_plane(const MeshSnapshot& m,const MirrorPlane& plane,double tol,std::optional<double> far_tol) {
    PlaneScore out;const auto& v=m.input().positions;
    if(v.empty())return out;
    const double diag=bbox_diagonal(m),eps=tol*diag,far_eps=(far_tol?*far_tol:5.0*tol)*diag;
    std::vector<Vec3> q(v.size());for(std::size_t i=0;i<v.size();++i)q[i]=plane.reflect(v[i]);
    const double reach=2.0*m.mean_edge_length();
    auto sd=surface_distance(m,q,far_eps>reach?far_eps:reach);
    out.matched.resize(v.size());std::size_t matched=0,far=0;
    for(std::size_t i=0;i<v.size();++i) {out.matched[i]=sd.distance[i]<=eps;matched+=out.matched[i];far+=sd.distance[i]>far_eps;}
    out.score=static_cast<double>(matched)/static_cast<double>(v.size());out.far=static_cast<double>(far)/static_cast<double>(v.size());
    out.foot=std::move(sd.foot);
    return out;
}
MirrorPlane refine_plane(const MeshSnapshot& m,const MirrorPlane& plane,std::span<const std::uint8_t> matched,std::span<const Vec3> foot) {
    const auto& v=m.input().positions;
    require(matched.size()==v.size() && foot.size()==v.size(),"symmetry.refine","matched/foot do not match the mesh vertices");
    std::vector<Vec3> P,S;for(std::size_t i=0;i<v.size();++i)if(matched[i]) {P.push_back(v[i]);S.push_back(foot[i]);}
    if(P.size()<3)return plane;
    std::vector<Vec3> d(P.size());for(std::size_t i=0;i<P.size();++i)d[i]=sub(P[i],S[i]);
    const auto eig=detail::lapack_syevd3(gram(d,{}));
    Vec3 n{eig.vectors[0][2],eig.vectors[1][2],eig.vectors[2][2]};
    if(dot(n,plane.normal)<0)n={-n[0],-n[1],-n[2]};
    n=unit(n);
    std::vector<double> along(P.size());
    for(std::size_t i=0;i<P.size();++i) {
        const Vec3 mid{0.5*(P[i][0]+S[i][0]),0.5*(P[i][1]+S[i][1]),0.5*(P[i][2]+S[i][2])};along[i]=dot(mid,n);
    }
    return {n,detail::numpy_sum(along)/static_cast<double>(along.size()),0.0,plane.source,0.0};
}
std::vector<MirrorPlane> detect_planes(const MeshSnapshot& m,double tol,double accept,double max_far,bool refine) {
    std::vector<MirrorPlane> found;const double diag=bbox_diagonal(m);
    for(const auto& candidate:candidate_planes(m)) {
        const auto scored=score_plane(m,candidate,tol);
        MirrorPlane plane{candidate.normal,candidate.offset,scored.score,candidate.source,scored.far};
        const bool any=std::any_of(scored.matched.begin(),scored.matched.end(),[](auto x) {return x!=0;});
        if(refine && any && scored.score>=accept-0.25) {
            const auto refined=refine_plane(m,candidate,scored.matched,scored.foot);
            const auto again=score_plane(m,refined,tol);
            if(again.score>=scored.score)plane={refined.normal,refined.offset,again.score,candidate.source,again.far};
        }
        if(plane.score<accept || plane.far>max_far)continue;
        bool duplicate=false;
        for(auto& other:found) {
            if(std::abs(dot(other.normal,plane.normal))>0.995 && std::abs(std::abs(other.offset)-std::abs(plane.offset))<=tol*diag+1e-12) {
                duplicate=true;
                if(plane.score>other.score)other=plane;
                break;
            }
        }
        if(!duplicate)found.push_back(plane);
    }
    std::stable_sort(found.begin(),found.end(),[](const MirrorPlane& a,const MirrorPlane& b) {return -a.score<-b.score;});
    return found;
}
std::optional<MirrorPlane> detect_plane(const MeshSnapshot& m,double tol,double accept,double max_far,bool refine) {
    auto planes=detect_planes(m,tol,accept,max_far,refine);
    if(planes.empty())return std::nullopt;
    return planes.front();
}
// --------------------------------------------------------------------- maps
MirrorMaps mirror_maps(const MeshSnapshot& m,const MirrorPlane& plane,double tol) {
    const auto& v=m.input().positions;const auto& t=m.tables();const auto nf=face_count(m);
    const double eps=tol*bbox_diagonal(m);
    std::vector<Vec3> q(v.size());for(std::size_t i=0;i<v.size();++i)q[i]=plane.reflect(v[i]);
    MirrorMaps out;out.vertex=nearest_points(v,q,eps<1e-12?1e-12:eps).index;
    out.edge.assign(t.edge_vertices.size(),-1);out.face.assign(nf,-1);
    for(std::size_t e=0;e<t.edge_vertices.size();++e) {
        const auto a=out.vertex[t.edge_vertices[e][0]],b=out.vertex[t.edge_vertices[e][1]];
        if(a>=0 && b>=0)out.edge[e]=m.edge_index(a,b);
    }
    // Sorted vertex tuples; a repeated tuple keeps the last face, as the source's dict.
    std::map<std::vector<TransferIndex>,TransferIndex> key;
    auto face_vertices=[&](std::size_t f) {
        std::vector<TransferIndex> ids;for(auto c:m.face_corners(static_cast<MeshIndex>(f)))ids.push_back(m.input().corner_vertices[c]);return ids;
    };
    for(std::size_t f=0;f<nf;++f) {auto ids=face_vertices(f);std::sort(ids.begin(),ids.end());key[ids]=static_cast<TransferIndex>(f);}
    for(std::size_t f=0;f<nf;++f) {
        auto image=face_vertices(f);bool missing=false;
        for(auto& id:image) {id=out.vertex[static_cast<std::size_t>(id)];missing=missing || id<0;}
        if(missing)continue;
        std::sort(image.begin(),image.end());
        const auto found=key.find(image);if(found!=key.end())out.face[f]=found->second;
    }
    return out;
}
bool is_combinatorial(const MirrorMaps& maps) {
    const auto& sv=maps.vertex;
    auto negative=[](const std::vector<TransferIndex>& x) {return std::any_of(x.begin(),x.end(),[](auto i) {return i<0;});};
    if(sv.empty() || negative(sv) || negative(maps.edge) || negative(maps.face))return false;
    for(std::size_t i=0;i<sv.size();++i)if(sv[static_cast<std::size_t>(sv[i])]!=static_cast<TransferIndex>(i))return false;
    return true;
}
// ------------------------------------------------------------- seam symmetry
SeamSymmetry seam_symmetry(const MeshSnapshot& m,std::span<const std::uint8_t> seams,const std::optional<MirrorPlane>& plane,double tol_edges,std::int64_t samples) {
    const auto& t=m.tables();const auto& g=m.geometry();const auto ne=t.edge_vertices.size();
    require(seams.size()==ne,"symmetry.seams","edge mask does not match the mesh edges");
    SeamSymmetry out;out.mask.assign(ne,0);
    std::vector<MeshIndex> edges;for(std::size_t e=0;e<ne;++e)if(seams[e])edges.push_back(static_cast<MeshIndex>(e));
    if(edges.empty() || !plane)return out;
    // numpy.linspace(0, 1, n): k * (1 / (n - 1)) with the last sample exactly 1.
    const auto n=static_cast<std::size_t>(std::max<std::int64_t>(samples,2));const double step=1.0/static_cast<double>(n-1);
    std::vector<double> ts(n);for(std::size_t k=0;k<n;++k)ts[k]=static_cast<double>(k)*step;ts[n-1]=1.0;
    std::vector<Vec3> sample_points,mids;
    for(auto e:edges) {
        const auto a=position(m,t.edge_vertices[e][0]),b=position(m,t.edge_vertices[e][1]);const auto ba=sub(b,a);
        for(double s:ts)sample_points.push_back({a[0]+s*ba[0],a[1]+s*ba[1],a[2]+s*ba[2]});
        mids.push_back(plane->reflect({0.5*(a[0]+b[0]),0.5*(a[1]+b[1]),0.5*(a[2]+b[2])}));
    }
    const double tol=tol_edges*positive_edge_median(m);
    const auto near=nearest_points(sample_points,mids,tol);
    std::vector<double> lengths,kept;
    for(std::size_t i=0;i<edges.size();++i) {
        lengths.push_back(g.edge_length[edges[i]]);
        if(near.distance[i]<=tol) {kept.push_back(g.edge_length[edges[i]]);out.mask[edges[i]]=1;}
    }
    const double total=detail::numpy_sum(lengths);
    out.score=total>0?detail::numpy_sum(kept)/total:nan;
    return out;
}
// ------------------------------------------------------ LAPACK dsyevd, n = 3
namespace detail {
namespace {
template<class A> auto& ix(A& array,int i) {return array[static_cast<std::size_t>(i)];}
constexpr double lapack_eps=0x1p-53;             // dlamch('E')
constexpr double lapack_safmin=0x1p-1022;         // dlamch('S')
double lapack_sign(double a,double b) {return b>=0?std::abs(a):-std::abs(a);}
double dlapy2(double x,double y) {
    const double xa=std::abs(x),ya=std::abs(y),w=std::max(xa,ya),z=std::min(xa,ya);
    if(z==0)return w;
    const double q=z/w;return w*std::sqrt(q*q+1.0);
}
// dlartg (LAPACK 3.9 and earlier; the scaling branches only apply near under/overflow).
void dlartg(double f,double g,double& cs,double& sn,double& r) {
    if(g==0) {cs=1;sn=0;r=f;return;}
    if(f==0) {cs=0;sn=1;r=g;return;}
    r=std::sqrt(f*f+g*g);cs=f/r;sn=g/r;
    if(std::abs(f)>std::abs(g) && cs<0) {cs=-cs;sn=-sn;r=-r;}
}
void dlaev2(double a,double b,double c,double& rt1,double& rt2,double& cs1,double& sn1) {
    const double sm=a+c,df=a-c,adf=std::abs(df),tb=b+b,ab=std::abs(tb);
    const double acmx=std::abs(a)>std::abs(c)?a:c,acmn=std::abs(a)>std::abs(c)?c:a;
    double rt;
    if(adf>ab) {const double q=ab/adf;rt=adf*std::sqrt(q*q+1.0);}
    else if(adf<ab) {const double q=adf/ab;rt=ab*std::sqrt(q*q+1.0);}
    else rt=ab*std::sqrt(2.0);
    int sgn1;
    if(sm<0) {rt1=0.5*(sm-rt);sgn1=-1;rt2=(acmx/rt1)*acmn-(b/rt1)*b;}
    else if(sm>0) {rt1=0.5*(sm+rt);sgn1=1;rt2=(acmx/rt1)*acmn-(b/rt1)*b;}
    else {rt1=0.5*rt;rt2=-0.5*rt;sgn1=1;}
    int sgn2;double cs;
    if(df>=0) {cs=df+rt;sgn2=1;} else {cs=df-rt;sgn2=-1;}
    if(std::abs(cs)>ab) {const double ct=-tb/cs;sn1=1.0/std::sqrt(1.0+ct*ct);cs1=ct*sn1;}
    else if(ab==0) {cs1=1;sn1=0;}
    else {const double tn=-cs/tb;cs1=1.0/std::sqrt(1.0+tn*tn);sn1=tn*cs1;}
    if(sgn1==sgn2) {const double tn=cs1;cs1=-sn1;sn1=tn;}
}
using Z3=std::array<std::array<double,3>,3>;
// dlasr('R', 'V', direct) on columns first..first+count-1 of Z.
void dlasr(Z3& z,bool forward,const double* c,const double* s,int first,int count) {
    auto apply=[&](int j) {
        const double ct=c[j],st=s[j];
        if(ct==1 && st==0)return;
        for(int i=0;i<3;++i) {
            const double temp=ix(ix(z,i),first+j+1);
            ix(ix(z,i),first+j+1)=ct*temp-st*ix(ix(z,i),first+j);
            ix(ix(z,i),first+j)=st*temp+ct*ix(ix(z,i),first+j);
        }
    };
    if(forward)for(int j=0;j<count-1;++j)apply(j);
    else for(int j=count-2;j>=0;--j)apply(j);
}
// dsteqr('I') on the 3x3 tridiagonal (d, e): implicit QL or QR per unreduced block.
Z3 dsteqr(std::array<double,3>& d,std::array<double,2>& e) {
    constexpr int n=3;Z3 z{};for(int i=0;i<3;++i)ix(ix(z,i),i)=1.0;
    const double eps2=lapack_eps*lapack_eps;const int nmaxit=n*30;int jtot=0,l1=0;
    std::array<double,3> cw{},sw{};
    while(l1<=n-1) {
        if(l1>0)ix(e,l1-1)=0.0;
        int m=n-1;
        for(int mm=l1;mm<n-1;++mm) {
            const double tst=std::abs(ix(e,mm));
            if(tst==0) {m=mm;break;}
            if(tst<=(std::sqrt(std::abs(ix(d,mm)))*std::sqrt(std::abs(ix(d,mm+1))))*lapack_eps) {ix(e,mm)=0.0;m=mm;break;}
        }
        int l=l1;const int lsv=l;int lend=m;const int lendsv=lend;l1=m+1;
        if(lend==l)continue;
        double anorm=0;for(int i=l;i<=lend;++i)anorm=std::max(anorm,std::abs(ix(d,i)));for(int i=l;i<lend;++i)anorm=std::max(anorm,std::abs(ix(e,i)));
        if(anorm==0)continue;
        if(std::abs(ix(d,lend))<std::abs(ix(d,l))) {lend=lsv;l=lendsv;}
        if(lend>l) {
            // QL iteration
            while(true) {
                int mm=lend;
                if(l!=lend)for(int k=l;k<lend;++k) {const double tst=ix(e,k)*ix(e,k);if(tst<=(eps2*std::abs(ix(d,k)))*std::abs(ix(d,k+1))+lapack_safmin) {mm=k;break;}}
                m=mm;
                if(m<lend)ix(e,m)=0.0;
                double p=ix(d,l);
                if(m==l) {ix(d,l)=p;++l;if(l<=lend)continue;break;}
                if(m==l+1) {
                    double rt1,rt2,c,s;dlaev2(ix(d,l),ix(e,l),ix(d,l+1),rt1,rt2,c,s);
                    dlasr(z,false,&c,&s,l,2);
                    ix(d,l)=rt1;ix(d,l+1)=rt2;ix(e,l)=0.0;l+=2;
                    if(l<=lend)continue;
                    break;
                }
                if(jtot==nmaxit)break;
                ++jtot;
                double g=(ix(d,l+1)-p)/(ix(e,l)*2.0),r=dlapy2(g,1.0);
                g=ix(d,m)-p+ix(e,l)/(g+lapack_sign(r,g));
                double s=1.0,c=1.0;p=0.0;
                for(int i=m-1;i>=l;--i) {
                    const double f=s*ix(e,i),b=c*ix(e,i);
                    dlartg(g,f,c,s,r);
                    if(i!=m-1)ix(e,i+1)=r;
                    g=ix(d,i+1)-p;r=(ix(d,i)-g)*s+c*2.0*b;p=s*r;ix(d,i+1)=g+p;g=c*r-b;
                    ix(cw,i-l)=c;ix(sw,i-l)=-s;
                }
                dlasr(z,false,cw.data(),sw.data(),l,m-l+1);
                ix(d,l)-=p;ix(e,l)=g;
            }
        } else {
            // QR iteration
            while(true) {
                int mm=lend;
                if(l!=lend)for(int k=l;k>lend;--k) {const double tst=ix(e,k-1)*ix(e,k-1);if(tst<=(eps2*std::abs(ix(d,k)))*std::abs(ix(d,k-1))+lapack_safmin) {mm=k;break;}}
                m=mm;
                if(m>lend)ix(e,m-1)=0.0;
                double p=ix(d,l);
                if(m==l) {ix(d,l)=p;--l;if(l>=lend)continue;break;}
                if(m==l-1) {
                    double rt1,rt2,c,s;dlaev2(ix(d,l-1),ix(e,l-1),ix(d,l),rt1,rt2,c,s);
                    dlasr(z,true,&c,&s,l-1,2);
                    ix(d,l-1)=rt1;ix(d,l)=rt2;ix(e,l-1)=0.0;l-=2;
                    if(l>=lend)continue;
                    break;
                }
                if(jtot==nmaxit)break;
                ++jtot;
                double g=(ix(d,l-1)-p)/(ix(e,l-1)*2.0),r=dlapy2(g,1.0);
                g=ix(d,m)-p+ix(e,l-1)/(g+lapack_sign(r,g));
                double s=1.0,c=1.0;p=0.0;
                for(int i=m;i<l;++i) {
                    const double f=s*ix(e,i),b=c*ix(e,i);
                    dlartg(g,f,c,s,r);
                    if(i!=m)ix(e,i-1)=r;
                    g=ix(d,i)-p;r=(ix(d,i+1)-g)*s+c*2.0*b;p=s*r;ix(d,i)=g+p;g=c*r-b;
                    ix(cw,i-m)=c;ix(sw,i-m)=s;
                }
                dlasr(z,true,cw.data(),sw.data(),m,l-m+1);
                ix(d,l)-=p;ix(e,l-1)=g;
            }
        }
        if(jtot>=nmaxit)throw ContractError(ErrorCode::numerical_failure,"symmetry.eigen","symmetric eigen decomposition did not converge");
    }
    // Selection sort, ascending, swapping eigenvector columns.
    for(int ii=1;ii<n;++ii) {
        const int i=ii-1;int k=i;double p=ix(d,i);
        for(int j=ii;j<n;++j)if(ix(d,j)<p) {k=j;p=ix(d,j);}
        if(k!=i) {ix(d,k)=ix(d,i);ix(d,i)=p;for(int r=0;r<3;++r)std::swap(ix(ix(z,r),i),ix(ix(z,r),k));}
    }
    return z;
}
}
SymmetricEigen3 lapack_syevd3(const std::array<std::array<double,3>,3>& lower) {
    for(const auto& row:lower)for(double x:row)
        if(!std::isfinite(x))throw ContractError(ErrorCode::numerical_failure,"symmetry.eigen","nonfinite matrix");
    auto a=lower;
    // dsytd2 (lower), i = 1: reflector annihilating a(3,1) against a(2,1).
    double alpha=ix(ix(a,1),0);const double x=ix(ix(a,2),0),xnorm=std::abs(x);double tau=0;
    if(xnorm!=0) {
        const double beta=-lapack_sign(dlapy2(alpha,xnorm),alpha);
        tau=(beta-alpha)/beta;ix(ix(a,2),0)=x*(1.0/(alpha-beta));alpha=beta;
    }
    std::array<double,3> d{};std::array<double,2> e{};
    ix(e,0)=alpha;
    if(tau!=0) {
        const std::array<double,2> v{1.0,ix(ix(a,2),0)};
        // dsymv (lower, 2x2 trailing block, alpha = tau, beta = 0)
        std::array<double,2> y{0.0,0.0};
        for(int j=0;j<2;++j) {
            const double temp1=tau*ix(v,j);double temp2=0;
            ix(y,j)=ix(y,j)+temp1*ix(ix(a,1+j),1+j);
            for(int i=j+1;i<2;++i) {ix(y,i)=ix(y,i)+temp1*ix(ix(a,1+i),1+j);temp2=temp2+ix(ix(a,1+i),1+j)*ix(v,i);}
            ix(y,j)=ix(y,j)+tau*temp2;
        }
        const double shift=tau*-0.5*(ix(y,0)*ix(v,0)+ix(y,1)*ix(v,1));
        for(int i=0;i<2;++i)ix(y,i)=ix(y,i)+shift*ix(v,i);
        // dsyr2 (lower, alpha = -1)
        for(int j=0;j<2;++j) {
            if(ix(v,j)==0 && ix(y,j)==0)continue;
            const double temp1=-1.0*ix(y,j),temp2=-1.0*ix(v,j);
            for(int i=j;i<2;++i)ix(ix(a,1+i),1+j)=ix(ix(a,1+i),1+j)+ix(v,i)*temp1+ix(y,i)*temp2;
        }
    }
    ix(d,0)=ix(ix(a,0),0);ix(e,1)=ix(ix(a,2),1);ix(d,1)=ix(ix(a,1),1);ix(d,2)=ix(ix(a,2),2);
    auto z=dsteqr(d,e);
    // dorm2r: Q = H(1) H(2) with tau(2) = 0, applied to rows 2..3.
    if(tau!=0) {
        const std::array<double,2> v{1.0,ix(ix(a,2),0)};const int lastv=ix(v,1)!=0?2:1;
        std::array<double,3> work{};
        for(int j=0;j<3;++j) {double temp=0;for(int i=0;i<lastv;++i)temp=temp+ix(ix(z,1+i),j)*ix(v,i);ix(work,j)=ix(work,j)+temp;}
        for(int j=0;j<3;++j) {
            if(ix(work,j)==0)continue;
            const double temp=-tau*ix(work,j);
            for(int i=0;i<lastv;++i)ix(ix(z,1+i),j)=ix(ix(z,1+i),j)+ix(v,i)*temp;
        }
    }
    SymmetricEigen3 out;out.values=d;out.vectors=z;
    return out;
}
}
}
