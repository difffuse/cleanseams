// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/axis.hpp"
#include "cleanseams/cache.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
namespace cleanseams {
namespace {
using detail::sub;using detail::dot;using detail::norm;using detail::cross;
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
std::size_t edge_count(const MeshSnapshot& m) {return m.tables().edge_vertices.size();}
std::size_t face_count(const MeshSnapshot& m) {return m.input().face_offsets.size()-1;}
void check_mask(const MeshSnapshot& m,std::span<const std::uint8_t> mask,const char* field) {
    require(mask.size()==edge_count(m),field,"edge mask does not match the mesh edges");
}
double ptp(const std::vector<double>& values,std::span<const MeshIndex> ids) {
    double low=INFINITY,high=-INFINITY;for(auto i:ids) {low=std::min(low,values[i]);high=std::max(high,values[i]);}
    return high-low;
}
std::vector<MeshIndex> sorted_unique(std::vector<MeshIndex> ids) {
    std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());return ids;
}
std::vector<MeshIndex> chart_vertices(const MeshSnapshot& m,std::span<const MeshIndex> faces) {
    std::vector<MeshIndex> out;
    for(auto f:faces)for(auto c:m.face_corners(f))out.push_back(static_cast<MeshIndex>(m.input().corner_vertices[c]));
    return sorted_unique(std::move(out));
}
Eigen::Matrix3d cross_matrix(const Eigen::Vector3d& p) {
    Eigen::Matrix3d m;m<<0,-p[2],p[1],p[2],0,-p[0],-p[1],p[0],0;return m;
}
// numpy.linalg.pinv: singular values above rcond times the largest are inverted.
Eigen::Matrix3d pinv(const Eigen::Matrix3d& a,double rcond) {
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(a,Eigen::ComputeFullU|Eigen::ComputeFullV);
    const auto& s=svd.singularValues();const double cutoff=rcond*s.maxCoeff();
    Eigen::Vector3d inv;for(int i=0;i<3;++i)inv[i]=s[i]>cutoff?1.0/s[i]:0.0;
    return svd.matrixV()*inv.asDiagonal()*svd.matrixU().transpose();
}
// Symmetric eigen decomposition, ascending eigenvalues (numpy.linalg.eigh).
template<class M> std::pair<Eigen::VectorXd,Eigen::MatrixXd> eigh(const M& a) {
    const Eigen::MatrixXd symmetric=(a+a.transpose())/2.0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(symmetric);
    if(solver.info()!=Eigen::Success)throw ContractError(ErrorCode::numerical_failure,"axis.fit","eigen decomposition failed");
    return {solver.eigenvalues(),solver.eigenvectors()};
}
struct FitMemo {
    std::map<std::pair<std::vector<MeshIndex>,double>,std::optional<RevolutionFrame>> frames;
    FitStats stats;
    const std::optional<RevolutionFrame>& get(const MeshSnapshot& m,std::span<const MeshIndex> faces,double limit) {
        std::pair key{std::vector<MeshIndex>(faces.begin(),faces.end()),limit};
        if(auto found=frames.find(key);found!=frames.end()) {++stats.reused;return found->second;}
        ++stats.fits;return frames.emplace(std::move(key),revolution_frame(m,faces,limit)).first->second;
    }
};
}
std::optional<LinearComplexFit> fit_linear_complex(std::span<const Vec3> p_in,std::span<const Vec3> n_in,
                                                   std::optional<std::span<const double>> weights_in,double residual_limit,double pitch_limit) {
    require(p_in.size()==n_in.size(),"axis.fit","positions and normals must have matching (n, 3) shapes");
    std::vector<double> weights(p_in.size(),1.0);
    if(weights_in) {
        require(weights_in->size()==p_in.size(),"axis.fit","weights must be finite and non-negative");
        for(std::size_t i=0;i<weights.size();++i) {
            const double w=(*weights_in)[i];require(std::isfinite(w) && w>=0,"axis.fit","weights must be finite and non-negative");weights[i]=w;
        }
    }
    std::vector<Vec3> P,N;std::vector<double> w;
    for(std::size_t i=0;i<p_in.size();++i) {
        const auto& p=p_in[i];const auto& n=n_in[i];
        const bool finite=std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) && std::isfinite(n[0]) && std::isfinite(n[1]) && std::isfinite(n[2]);
        if(finite && weights[i]>0 && norm(n)>1e-12) {P.push_back(p);N.push_back(n);w.push_back(weights[i]);}
    }
    if(P.size()<6)return std::nullopt;
    const auto n=P.size();const double total=detail::numpy_sum(w);for(auto& x:w)x/=total;
    Vec3 centre{0,0,0};for(std::size_t i=0;i<n;++i)for(std::size_t d=0;d<3;++d)centre[d]+=w[i]*P[i][d];
    double spread=0;for(std::size_t i=0;i<n;++i) {const auto q=sub(P[i],centre);spread+=w[i]*dot(q,q);}
    const double scale=std::sqrt(spread);
    if(scale<=1e-14)return std::nullopt;
    const auto rows=static_cast<Eigen::Index>(n);Eigen::MatrixXd X(rows,3),Nn(rows,3),q(rows,6);
    std::vector<Vec3> xs(n),ns(n);
    for(std::size_t i=0;i<n;++i) {
        const auto x=detail::divide(sub(P[i],centre),scale);const auto nn=detail::divide(N[i],norm(N[i]));const auto c=cross(x,nn);
        xs[i]=x;ns[i]=nn;const auto r=static_cast<Eigen::Index>(i);
        for(std::size_t d=0;d<3;++d) {const auto k=static_cast<Eigen::Index>(d);X(r,k)=x[d];Nn(r,k)=nn[d];q(r,k)=c[d];q(r,3+k)=nn[d];}
    }
    const Eigen::Map<const Eigen::VectorXd> wv(w.data(),static_cast<Eigen::Index>(n));
    const Eigen::MatrixXd M=(q.array().colwise()*wv.array()).matrix().transpose()*q;
    Eigen::Matrix<double,6,6> normaliser=Eigen::Matrix<double,6,6>::Zero();
    const double radial=(X.array().square().rowwise().sum()*wv.array()).sum();
    normaliser.topLeftCorner<3,3>()=(1+radial)*Eigen::Matrix3d::Identity()-(Nn.array().colwise()*wv.array()).matrix().transpose()*Nn-
        (X.array().colwise()*wv.array()).matrix().transpose()*X;
    normaliser.topRightCorner<3,3>()=cross_matrix(X.transpose()*wv);
    normaliser.bottomLeftCorner<3,3>()=normaliser.topRightCorner<3,3>().transpose();
    normaliser.bottomRightCorner<3,3>()=Eigen::Matrix3d::Identity();
    // A pure cylinder also has a translational field along its axis: remove the free
    // minimum-norm translation instead of mistaking it for screw pitch.
    const Eigen::Matrix3d translation=-pinv(M.bottomRightCorner<3,3>(),1e-3)*M.bottomLeftCorner<3,3>();
    Eigen::Matrix<double,6,3> B;B.topRows<3>()=Eigen::Matrix3d::Identity();B.bottomRows<3>()=translation;
    const Eigen::Matrix3d reduced=B.transpose()*M*B,normal=B.transpose()*normaliser*B;
    const auto [values,vectors]=eigh(normal);
    if(values.minCoeff()<=1e-12)return std::nullopt;
    Eigen::MatrixXd scaled=vectors;for(int k=0;k<3;++k)scaled.col(k)/=std::sqrt(values[k]);
    const Eigen::Matrix3d whitening=scaled*vectors.transpose();
    const auto [eigenvalues,modes]=eigh(whitening*reduced*whitening);
    // A sphere and other ambiguous samples admit several rotation axes.
    if(eigenvalues[1]<=std::max(1e-9,5*std::max(0.0,eigenvalues[0])))return std::nullopt;
    const Eigen::VectorXd field=B*whitening*modes.col(0);
    Vec3 omega{field[0],field[1],field[2]},velocity{field[3],field[4],field[5]};
    const double length=norm(omega);
    if(length<=1e-12)return std::nullopt;
    omega=detail::divide(omega,length);velocity=detail::divide(velocity,length);
    std::size_t axis=0;for(std::size_t d=1;d<3;++d)if(std::abs(omega[d])>std::abs(omega[axis]))axis=d;
    if(omega[axis]<0) {omega=detail::mul(omega,-1.0);velocity=detail::mul(velocity,-1.0);}
    const double pitch=dot(omega,velocity);
    double along=0,magnitude=0;
    for(std::size_t i=0;i<n;++i) {
        const auto& x=xs[i];const auto& nn=ns[i];const auto tangent=detail::add(cross(omega,x),velocity);
        const double a=dot(nn,tangent);along+=w[i]*(a*a);magnitude+=w[i]*dot(tangent,tangent);
    }
    const double residual=std::sqrt(along/std::max(magnitude,1e-15));
    if(residual>residual_limit || std::abs(pitch)>pitch_limit)return std::nullopt;
    LinearComplexFit out;out.point=detail::add(centre,detail::mul(cross(omega,velocity),scale));out.direction=omega;
    out.pitch=pitch*scale;out.residual=residual;out.scale=scale;
    for(int k=0;k<3;++k)out.eigenvalues[static_cast<std::size_t>(k)]=eigenvalues[k];
    return out;
}
std::optional<RevolutionFrame> revolution_frame(const MeshSnapshot& m,std::optional<std::span<const MeshIndex>> faces_in,double residual_limit) {
    const auto& g=m.geometry();const auto& t=m.tables();const auto nf=face_count(m);
    std::vector<MeshIndex> faces;
    if(faces_in) {for(auto f:*faces_in) {require(f<nf,"axis.faces","face index out of range");faces.push_back(f);}}
    else for(std::size_t f=0;f<nf;++f)faces.push_back(static_cast<MeshIndex>(f));
    std::vector<double> weights;for(auto f:faces)weights.push_back(g.face_area[f]);
    if(faces.empty())return std::nullopt;
    const double total=detail::numpy_sum(weights);
    if(total<=0)return std::nullopt;
    Vec3 average{0,0,0};for(std::size_t i=0;i<faces.size();++i)for(std::size_t d=0;d<3;++d)average[d]+=g.face_normal[faces[i]][d]*weights[i];
    if(norm(detail::divide(average,total))>.9)return std::nullopt;
    std::vector<Vec3> positions(g.face_center),normals(g.face_normal);
    std::vector<std::uint8_t> selected(nf,0);for(auto f:faces)selected[f]=1;
    std::vector<std::array<MeshIndex,2>> pairs;
    for(std::size_t e=0;e<edge_count(m);++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()==2 && selected[pair[0]] && selected[pair[1]])pairs.push_back({pair[0],pair[1]});
    }
    if(!pairs.empty()) {
        std::vector<MeshIndex> a,b;
        for(const auto& p:pairs) {a.push_back(p[0]);b.push_back(p[1]);}
        for(const auto& p:pairs) {a.push_back(p[1]);b.push_back(p[0]);}
        for(int round=0;round<2;++round) {
            std::vector<double> mass(g.face_area);std::vector<Vec3> P(nf),N(nf);
            for(std::size_t f=0;f<nf;++f) {P[f]=detail::mul(positions[f],mass[f]);N[f]=detail::mul(normals[f],mass[f]);}
            for(std::size_t k=0;k<a.size();++k) {
                const auto area=g.face_area[b[k]];
                for(std::size_t d=0;d<3;++d) {P[a[k]][d]+=positions[b[k]][d]*area;}
            }
            for(std::size_t k=0;k<a.size();++k) {
                const auto area=g.face_area[b[k]];
                for(std::size_t d=0;d<3;++d) {N[a[k]][d]+=normals[b[k]][d]*area;}
            }
            for(std::size_t k=0;k<a.size();++k)mass[a[k]]+=g.face_area[b[k]];
            for(std::size_t f=0;f<nf;++f) {
                positions[f]=detail::divide(P[f],std::max(mass[f],1e-300));normals[f]=detail::divide(N[f],std::max(norm(N[f]),1e-300));
            }
        }
    }
    std::vector<Vec3> P,N;for(auto f:faces) {P.push_back(positions[f]);N.push_back(normals[f]);}
    auto fit=fit_linear_complex(P,N,std::span<const double>(weights),residual_limit);
    if(!fit)return std::nullopt;
    RevolutionFrame out;static_cast<LinearComplexFit&>(out)=*fit;
    const auto& d=fit->direction;std::size_t axis=0;for(std::size_t k=1;k<3;++k)if(std::abs(d[k])<std::abs(d[axis]))axis=k;
    Vec3 helper{0,0,0};helper[axis]=1;
    auto u=cross(d,helper);u=detail::divide(u,norm(u));const auto v=cross(d,u);
    for(const auto& p:m.input().positions) {
        const auto delta=sub(p,fit->point);const double h=dot(delta,d);
        out.theta.push_back(std::atan2(dot(delta,v),dot(delta,u)));out.height.push_back(h);
        out.radius.push_back(norm(sub(delta,detail::mul(d,h))));
    }
    return out;
}
TracePath trace_path(const MeshSnapshot& m,std::span<const std::uint8_t> cut,MeshIndex vertex,MeshIndex edge) {
    check_mask(m,cut,"axis.cut");const auto& t=m.tables();require(edge<edge_count(m),"axis.edge","edge index out of range");
    TracePath out;std::set<MeshIndex> seen{vertex};
    while(!cut[edge]) {
        out.edges.push_back(edge);
        const auto arrival=m.edge_other(edge,vertex);
        const bool blocked=std::any_of(t.vert_edges[arrival].begin(),t.vert_edges[arrival].end(),[&](auto e) {return cut[e]!=0;});
        if(seen.contains(arrival) || blocked) {out.complete=true;return out;}
        seen.insert(arrival);
        if(!t.vert_regular[arrival] || !t.fans[arrival])return out;
        const auto following=t.fans[arrival]->straight(edge);
        if(following<0)return out;
        vertex=arrival;edge=static_cast<MeshIndex>(following);
    }
    out.complete=true;return out;
}
std::vector<std::uint8_t> trace_directed(const MeshSnapshot& m,std::span<const std::uint8_t> cut,std::span<const DirectedSeed> seeds) {
    check_mask(m,cut,"axis.cut");const auto& t=m.tables();
    std::vector<std::uint8_t> result(cut.begin(),cut.end());
    for(const auto& seed:seeds) {
        if(seed.edge<0 || static_cast<std::uint64_t>(seed.edge)>=edge_count(m))continue;
        const auto edge=static_cast<MeshIndex>(seed.edge);const auto& ends=t.edge_vertices[edge];
        if((seed.vertex!=ends[0] && seed.vertex!=ends[1]) || result[edge])continue;
        const auto path=trace_path(m,result,static_cast<MeshIndex>(seed.vertex),edge);
        if(path.complete)for(auto e:path.edges)result[e]=1;
    }
    return result;
}
std::vector<std::uint8_t> shoulder_seed(const CutSurface& cs,MeshIndex chart,const RevolutionFrame& frame) {
    const auto& m=cs.mesh();const auto& ct=cs.tables();const auto& t=m.tables();const auto& g=m.geometry();
    require(chart<ct.n_charts,"axis.chart","chart index out of range");
    const auto ne=edge_count(m),nf=face_count(m);std::vector<std::uint8_t> cuts(ne,0);
    // A closed handle has no caps: its horizontal tangent bands are not shoulders.
    if(ct.chart_loops[chart].empty() && ct.chart_euler[chart]!=2)return cuts;
    const auto faces=ct.chart_faces[chart];
    std::vector<double> alignment(nf);for(std::size_t f=0;f<nf;++f)alignment[f]=std::abs(dot(g.face_normal[f],frame.direction));
    double high=-INFINITY,low=INFINITY;for(auto f:faces) {high=std::max(high,alignment[f]);low=std::min(low,alignment[f]);}
    if(!(high>std::cos(std::numbers::pi/6) && low<std::cos(std::numbers::pi/3)))return cuts;
    std::vector<std::uint8_t> selected(nf,0);for(auto f:faces)selected[f]=1;
    const double half=std::sqrt(0.5);
    for(std::size_t e=0;e<ne;++e) {
        const auto pair=t.edge_faces[e];
        if(pair.size()==2 && selected[pair[0]] && selected[pair[1]])cuts[e]=(alignment[pair[0]]>half)!=(alignment[pair[1]]>half);
    }
    // Retain closed rims only: an incomplete shoulder cannot introduce a loose stub.
    for(bool changed=true;changed;) {
        std::vector<std::uint64_t> degree(m.input().positions.size(),0);
        for(std::size_t e=0;e<ne;++e)if(cuts[e])for(auto v:t.edge_vertices[e])++degree[v];
        changed=false;
        for(std::size_t e=0;e<ne;++e)if(cuts[e] && (degree[t.edge_vertices[e][0]]==1 || degree[t.edge_vertices[e][1]]==1)) {cuts[e]=0;changed=true;}
    }
    return cuts;
}
std::optional<DirectedSeed> generatrix_seed(const CutSurface& cs,MeshIndex chart,const RevolutionFrame& frame,
                                            std::optional<std::span<const double>> cost,std::optional<Vec3> hide_dir) {
    const auto& m=cs.mesh();const auto& ct=cs.tables();const auto& t=m.tables();const auto& g=m.geometry();
    require(chart<ct.n_charts,"axis.chart","chart index out of range");
    if(cost)require(cost->size()==edge_count(m),"axis.cost","cost does not match the mesh edges");
    if(ct.chart_loops[chart].size()!=2)return std::nullopt;
    std::vector<MeshIndex> candidates;
    for(auto c:ct.loops[ct.chart_loops[chart][0]])candidates.push_back(static_cast<MeshIndex>(m.input().corner_vertices[c]));
    candidates=sorted_unique(std::move(candidates));
    if(candidates.size()>64) {
        // numpy.linspace(0, n - 1, 64).astype(int): i * step truncated, exact endpoint.
        std::vector<MeshIndex> sampled;const double step=static_cast<double>(candidates.size()-1)/63.0;
        for(std::size_t i=0;i<64;++i)sampled.push_back(candidates[i==63?candidates.size()-1:static_cast<std::size_t>(static_cast<double>(i)*step)]);
        candidates=std::move(sampled);
    }
    const auto vertices_in_chart=chart_vertices(m,ct.chart_faces[chart]);
    const double required_span=0.8*ptp(frame.height,vertices_in_chart);
    std::optional<std::tuple<double,double,MeshIndex,MeshIndex>> best;
    for(auto vertex:candidates)for(auto edge:t.vert_edges[vertex]) {
        if(ct.cut[edge])continue;
        const auto owners=t.edge_faces[edge];
        if(std::none_of(owners.begin(),owners.end(),[&](auto f) {return ct.chart_of_face[f]==chart;}))continue;
        const auto path=trace_path(m,ct.cut,vertex,edge);
        if(!path.complete || path.edges.empty())continue;
        std::vector<MeshIndex> verts;for(auto e:path.edges)for(auto v:t.edge_vertices[e])verts.push_back(v);
        verts=sorted_unique(std::move(verts));
        if(ptp(frame.height,verts)<required_span)continue;
        std::vector<double> lengths,weighted;for(auto e:path.edges)lengths.push_back(g.edge_length[e]);
        double score=0;
        if(cost) {
            for(std::size_t k=0;k<path.edges.size();++k)weighted.push_back((*cost)[path.edges[k]]*lengths[k]);
            score=detail::numpy_sum(weighted)/detail::numpy_sum(lengths);
        } else if(hide_dir) {
            std::vector<double> facing;
            for(auto v:verts) {
                auto radial=sub(m.input().positions[v],frame.point);radial=sub(radial,detail::mul(frame.direction,dot(radial,frame.direction)));
                facing.push_back(dot(radial,*hide_dir));
            }
            score=detail::numpy_sum(facing)/static_cast<double>(facing.size());
        } else {
            for(std::size_t k=0;k<path.edges.size();++k)weighted.push_back(std::abs(g.edge_dihedral[path.edges[k]])*lengths[k]);
            score=detail::numpy_sum(weighted)/detail::numpy_sum(lengths);
        }
        const std::tuple<double,double,MeshIndex,MeshIndex> key{score,detail::numpy_sum(lengths),vertex,edge};
        if(!best || key<*best)best=key;
    }
    if(!best)return std::nullopt;
    return DirectedSeed{std::get<2>(*best),std::get<3>(*best)};
}
RevolutionCuts revolution_cuts(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> cut,double residual_limit,
                               std::optional<std::span<const double>> cost,std::optional<Vec3> hide_dir,const Deadline& deadline) {
    require(bool(mesh),"axis.mesh","missing mesh snapshot");const auto& m=*mesh;check_mask(m,cut,"axis.cut");
    RevolutionCuts out;FitMemo memo;std::vector<std::uint8_t> current(cut.begin(),cut.end());
    out.shoulders.assign(edge_count(m),0);
    {
        const auto cs_owner=cut_surface(mesh,current);const CutSurface& cs=*cs_owner;const auto& ct=cs.tables();
        for(MeshIndex chart=0;chart<ct.n_charts;++chart) {
            if(deadline.expired())break;
            // A singly open bottle can have several narrow lip/shoulder bands: keep its layout.
            const auto faces=ct.chart_faces[chart];
            if(faces.size()<12 || ct.chart_loops[chart].size()==1)continue;
            if(const auto& frame=memo.get(m,faces,residual_limit)) {
                const auto rim=shoulder_seed(cs,chart,*frame);
                for(std::size_t e=0;e<rim.size();++e)out.shoulders[e]=out.shoulders[e] || rim[e];
            }
        }
    }
    for(std::size_t e=0;e<current.size();++e)current[e]=current[e] || out.shoulders[e];
    const auto cs_owner=cut_surface(mesh,current);const CutSurface& cs=*cs_owner;const auto& ct=cs.tables();
    for(MeshIndex chart=0;chart<ct.n_charts;++chart) {
        if(deadline.expired())break;
        if(ct.chart_loops[chart].size()!=2)continue;
        const auto faces=ct.chart_faces[chart];
        const auto& frame=memo.get(m,faces,residual_limit);
        if(!frame)continue;
        const auto seed=generatrix_seed(cs,chart,*frame,cost,hide_dir);
        if(!seed)continue;
        auto proposed=trace_directed(m,current,std::span<const DirectedSeed>(&*seed,1));
        const auto trial_owner=cut_surface(mesh,proposed);const CutSurface& trial=*trial_owner;std::set<MeshIndex> affected;
        for(auto f:faces)affected.insert(trial.tables().chart_of_face[f]);
        if(affected.size()!=1 || !trial.tables().chart_is_disk[*affected.begin()])continue;
        RevolutionRecord record;record.seed=*seed;record.residual=frame->residual;record.pitch=frame->pitch;
        record.direction=frame->direction;record.point=frame->point;
        for(std::size_t e=0;e<proposed.size();++e)record.edges+=proposed[e] && !current[e];
        out.records.push_back(record);current=std::move(proposed);
    }
    out.cut=std::move(current);out.stats=memo.stats;
    return out;
}
LatheProfileCuts lathe_profile_cuts(const std::shared_ptr<const MeshSnapshot>& mesh,std::span<const std::uint8_t> initial_cut,
                                    std::optional<std::span<const double>> cost,std::optional<Vec3> hide_dir,const Deadline& deadline) {
    require(bool(mesh),"axis.mesh","missing mesh snapshot");const auto& m=*mesh;check_mask(m,initial_cut,"axis.cut");
    const auto& t=m.tables();const auto ne=edge_count(m),nv=m.input().positions.size();
    std::vector<std::uint8_t> original(ne);for(std::size_t e=0;e<ne;++e)original[e]=initial_cut[e] || t.edge_boundary[e];
    const auto cs_owner=cut_surface(mesh,original);const CutSurface& cs=*cs_owner;const auto& ct=cs.tables();
    LatheProfileCuts out;out.cut=original;out.protected_rims.assign(ne,0);FitMemo memo;
    Vec3 low{INFINITY,INFINITY,INFINITY},high{-INFINITY,-INFINITY,-INFINITY};
    for(const auto& p:m.input().positions)for(std::size_t d=0;d<3;++d) {low[d]=std::min(low[d],p[d]);high[d]=std::max(high[d],p[d]);}
    const double scale=std::max(nv?norm(sub(high,low)):0.0,1e-300);
    for(MeshIndex chart=0;chart<ct.n_charts;++chart) {
        if(deadline.expired())break;
        const auto faces=ct.chart_faces[chart];
        if(ct.chart_loops[chart].size()!=1 || faces.size()<24)continue;
        std::vector<MeshIndex> vertices,edges;
        for(auto f:faces)for(auto c:m.face_corners(f)) {vertices.push_back(static_cast<MeshIndex>(m.input().corner_vertices[c]));edges.push_back(t.corner_edge[c]);}
        vertices=sorted_unique(std::move(vertices));edges=sorted_unique(std::move(edges));
        if(std::any_of(edges.begin(),edges.end(),[&](auto e) {return original[e] && !t.edge_boundary[e];}))continue;
        const auto& frame=memo.get(m,faces,1e-4);
        if(!frame)continue;
        const auto& height=frame->height;const auto& radius=frame->radius;
        std::vector<MeshIndex> order;for(auto v:vertices)if(radius[v]>scale*1e-4)order.push_back(v);
        std::stable_sort(order.begin(),order.end(),[&](auto a,auto b) {return height[a]<height[b];});
        std::vector<std::vector<MeshIndex>> rows;
        for(std::size_t i=0;i<order.size();++i) {
            if(i==0 || height[order[i]]-height[order[i-1]]>scale*1e-5)rows.emplace_back();
            rows.back().push_back(order[i]);
        }
        if(rows.size()<4 || rows.size()>64)continue;
        std::vector<double> levels,radii;
        for(const auto& row:rows) {
            std::vector<double> h,r;for(auto v:row) {h.push_back(height[v]);r.push_back(radius[v]);}
            levels.push_back(detail::numpy_sum(h)/static_cast<double>(h.size()));radii.push_back(detail::numpy_sum(r)/static_cast<double>(r.size()));
        }
        if(std::any_of(rows.begin(),rows.end(),[&](const auto& row) {return ptp(radius,row)>scale*1e-4;}))continue;
        // A closed base contributes one pole or a flat first row; every other row is a complete ring.
        std::vector<std::int64_t> row_id(nv,-1);
        for(std::size_t i=0;i<rows.size();++i)for(auto v:rows[i])row_id[v]=static_cast<std::int64_t>(i);
        std::vector<std::vector<MeshIndex>> ring_edges;bool valid=true;
        for(std::size_t i=0;i<rows.size() && valid;++i) {
            std::vector<MeshIndex> ring;std::vector<std::uint64_t> degree(nv,0);
            for(auto e:edges) {
                const auto& ends=t.edge_vertices[e];
                if(row_id[ends[0]]==static_cast<std::int64_t>(i) && row_id[ends[1]]==static_cast<std::int64_t>(i)) {ring.push_back(e);++degree[ends[0]];++degree[ends[1]];}
            }
            if(rows[i].size()>1 && (rows[i].size()<8 || std::any_of(rows[i].begin(),rows[i].end(),[&](auto v) {return degree[v]!=2;})))valid=false;
            ring_edges.push_back(std::move(ring));
        }
        if(!valid)continue;
        // Orient the profile from its closed end towards its open mouth.
        std::set<std::int64_t> mouth_rows;
        for(auto c:ct.loops[ct.chart_loops[chart][0]])mouth_rows.insert(row_id[static_cast<std::size_t>(m.input().corner_vertices[c])]);
        const auto last=static_cast<std::int64_t>(rows.size())-1;
        if(mouth_rows.size()!=1 || (*mouth_rows.begin()!=0 && *mouth_rows.begin()!=last))continue;
        if(*mouth_rows.begin()==0) {
            std::reverse(levels.begin(),levels.end());std::reverse(radii.begin(),radii.end());std::reverse(ring_edges.begin(),ring_edges.end());
        }
        const auto segments=levels.size()-1;std::vector<double> arc(segments);std::vector<std::uint8_t> straight(segments);
        for(std::size_t i=0;i<segments;++i) {
            const double dz=std::abs(levels[i+1]-levels[i]),dr=radii[i+1]-radii[i];arc[i]=std::hypot(dz,dr);
            straight[i]=std::abs(dr)/std::max(arc[i],1e-300)<std::sin(20*(std::numbers::pi/180.0));
        }
        std::vector<std::pair<std::size_t,std::size_t>> groups;std::optional<std::size_t> start;
        for(std::size_t i=0;i<=segments;++i) {
            const bool flag=i<segments && straight[i];
            if(flag && !start)start=i;
            if(!flag && start) {groups.push_back({*start,i});start.reset();}
        }
        if(groups.empty())continue;
        auto range_sum=[&](const std::vector<double>& v,std::size_t a,std::size_t b) {return detail::numpy_sum(std::span<const double>(v).subspan(a,b-a));};
        std::size_t body=0;double best=-INFINITY;
        for(std::size_t k=0;k<groups.size();++k) {
            const auto [a,b]=groups[k];std::vector<double> terms;
            for(std::size_t i=a;i<b;++i)terms.push_back(arc[i]*(radii[i]+radii[i+1]));
            const double value=detail::numpy_sum(terms);if(value>best) {best=value;body=k;}
        }
        const auto [body_start,body_end]=groups[body];
        const double arc_total=detail::numpy_sum(arc);
        if(body_start==0 || range_sum(arc,body_start,body_end)<.25*arc_total)continue;
        std::set<std::int64_t> boundaries{static_cast<std::int64_t>(body_start)};
        if(static_cast<std::int64_t>(body_end)<last)boundaries.insert(static_cast<std::int64_t>(body_end));
        const double body_radius=range_sum(radii,body_start,body_end+1)/static_cast<double>(body_end+1-body_start);
        std::optional<std::pair<std::size_t,std::size_t>> neck;double neck_arc=-INFINITY;
        for(const auto& [a,b]:groups) {
            if(!(a>body_end))continue;const double length=range_sum(arc,a,b);
            if(!(length>=.08*arc_total && range_sum(radii,a,b+1)/static_cast<double>(b+1-a)<.75*body_radius))continue;
            if(length>neck_arc) {neck_arc=length;neck=std::pair{a,b};}
        }
        if(neck)boundaries.insert(static_cast<std::int64_t>(neck->first));
        std::vector<std::uint8_t> rims(ne,0);
        for(auto b:boundaries)for(auto e:ring_edges[static_cast<std::size_t>(b)])rims[e]=1;
        auto candidate=out.cut;for(std::size_t e=0;e<ne;++e)candidate[e]=candidate[e] || rims[e];
        const auto partition_owner=cut_surface(mesh,candidate);const CutSurface& partition=*partition_owner;std::set<MeshIndex> parts;
        for(auto f:faces)parts.insert(partition.tables().chart_of_face[f]);
        bool failed=false;
        for(auto part:parts) {
            if(partition.tables().chart_is_disk[part])continue;
            const auto seed=generatrix_seed(partition,part,*frame,cost,hide_dir);
            if(!seed) {failed=true;break;}
            candidate=trace_directed(m,candidate,std::span<const DirectedSeed>(&*seed,1));
        }
        const auto final_surface_owner=cut_surface(mesh,candidate);const CutSurface& final_surface=*final_surface_owner;std::set<MeshIndex> affected;
        for(auto f:faces)affected.insert(final_surface.tables().chart_of_face[f]);
        if(failed || std::any_of(affected.begin(),affected.end(),[&](auto c) {return !final_surface.tables().chart_is_disk[c];}))continue;
        out.cut=std::move(candidate);
        for(std::size_t e=0;e<ne;++e)out.protected_rims[e]=out.protected_rims[e] || rims[e];
        out.handled_faces.insert(out.handled_faces.end(),faces.begin(),faces.end());
        LatheRecord record;record.charts=affected.size();record.rows=rows.size();record.rim_rows.assign(boundaries.begin(),boundaries.end());
        record.residual=frame->residual;record.direction=frame->direction;out.records.push_back(std::move(record));
    }
    out.stats=memo.stats;
    return out;
}
}
