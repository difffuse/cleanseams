// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/cache.hpp"
#include "detail.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <set>
namespace cleanseams {
using namespace quality_detail;
TriangleFrame triangle_frame(const std::array<Vec3,3>& p) {
    const auto e1=sub(p[1],p[0]),e2=sub(p[2],p[0]);
    const double length=std::sqrt(dot(e1,e1));
    const auto x=div(e1,std::max(length,1e-300)),n=cross(e1,e2);
    TriangleFrame frame;frame.double_area=std::sqrt(dot(n,n));
    const auto y=cross(div(n,std::max(frame.double_area,1e-300)),x);
    // Source triangle_frames uses einsum("ij,ij->i") for these projections.
    frame.coordinates[1][0]=length;frame.coordinates[2]={detail::numpy_einsum_dot3(e2,x),detail::numpy_einsum_dot3(e2,y)};
    return frame;
}
QualityGeometry::QualityGeometry(std::shared_ptr<const CutSurface> cut,std::uint64_t budget):cut_(std::move(cut)),memory_budget_(budget) {
    require(bool(cut_),"quality.cut","cut surface owner required");detail::note_geometry_build(*cut_);
    const auto& mesh=cut_->mesh();const auto& input=mesh.input();const auto& tri=mesh.triangulation();
    const auto nf=input.face_offsets.size()-1,nc=input.corner_vertices.size(),nt=tri.corners.size();
    const auto allowance=validate_mesh_counts(input.positions.size(),nf,nc,input.edge_vertices.size(),budget);
    // The mesh tables plus the measured charge of the frames, chart rows and one full
    // snapshot with its report sort/aggregation scratch (Gate A F3). A run checks the
    // simultaneous mesh, cut, quality and seam tables up front (check_run_memory).
    const auto own=quality_table_bytes({input.positions.size(),nf,nc,input.edge_vertices.size()});
    if(own>budget-allowance)throw ContractError(ErrorCode::memory_limit,"quality.memory_budget","quality storage exceeds the memory budget");
    frames_.reserve(nt);
    for(const auto& c:tri.corners) {
        std::array<Vec3,3> p;for(std::size_t k=0;k<3;++k) p[k]=input.positions[static_cast<std::size_t>(input.corner_vertices[c[k]])];
        frames_.push_back(triangle_frame(p));
    }
    chart_triangles_=cut_detail::rows(cut_->tables().n_charts,[&](auto add){for(std::size_t t=0;t<nt;++t) add(cut_->tables().chart_of_face[tri.faces[t]],static_cast<MeshIndex>(t));});
    if(nf) {
        const auto area=detail::numpy_sum(mesh.geometry().face_area); // source np.mean(face_area): pairwise sum
        const double mean=area/static_cast<double>(nf)*2/std::max(1.0,static_cast<double>(nc)/static_cast<double>(nf)-2);
        sliver_threshold_=1e-6*std::max(mean,1e-300);
    }
}
std::vector<std::int8_t> chart_reference_orientation(const QualityGeometry& geometry,CornerUvView uv) {
    const auto& cs=geometry.cut();const auto& mesh=cs.mesh();validate_view(mesh,uv);
    std::vector<std::int8_t> result(cs.tables().n_charts,1);
    for(MeshIndex ch=0;ch<result.size();++ch) {
        bool any=false,negative=true;
        for(auto f:cs.tables().chart_faces[ch]) if(!mesh.geometry().face_degenerate_kind[f]) {
            any=true;negative=negative && signed_area(mesh,uv,f)<0;
        }
        if(any && negative) result[ch]=-1;
    }
    return result;
}
QualitySnapshot::QualitySnapshot(std::shared_ptr<const QualityGeometry> geometry,CornerUvView uv,const QualityPolicy& policy,
    std::optional<std::span<const MeshIndex>> charts):geometry_(std::move(geometry)) {
    require(bool(geometry_),"quality.geometry","geometry owner required");
    const auto& cs=geometry_->cut();const auto& mesh=cs.mesh();const auto& ct=cs.tables();const auto& mt=mesh.triangulation();
    validate_view(mesh,uv);const auto nf=mesh.input().face_offsets.size()-1;
    require(policy.exempt_input_degenerate || policy.exempt_faces.empty(),"quality.exempt_faces","raw measurement cannot add exemptions");
    require(policy.exempt_faces.size()<=geometry_->memory_budget()/sizeof(TransferIndex),"quality.exempt_faces","constraint input exceeds memory allowance");
    std::set<TransferIndex> extras;
    for(auto f:policy.exempt_faces) {require(f>=0 && static_cast<std::uint64_t>(f)<nf,"quality.exempt_faces","face out of range");extras.insert(f);}
    policy_.exempt_input_degenerate=policy.exempt_input_degenerate;policy_.uv_revision=policy.uv_revision;
    policy_.exempt_faces.assign(extras.begin(),extras.end());
    const auto& signs=policy.reference_orientation;
    require(signs.empty() || signs.size()==ct.n_charts,"quality.reference_orientation","one sign per chart required");
    for(auto sign:signs) require(sign==1 || sign==-1,"quality.reference_orientation","sign must be +1 or -1");
    policy_.reference_orientation=signs;
    auto get=[&](MeshIndex c){auto p=point(uv,c);if(!signs.empty()) p[0]*=signs[ct.chart_of_corner[c]];return p;};
    auto exempt=[&](MeshIndex f){return (policy_.exempt_input_degenerate && mesh.geometry().face_degenerate_kind[f]) || std::binary_search(policy_.exempt_faces.begin(),policy_.exempt_faces.end(),static_cast<TransferIndex>(f));};
    if(charts) {
        require(charts->size()<=geometry_->memory_budget()/sizeof(MeshIndex),"quality.charts","selection input exceeds memory allowance");
        std::set<MeshIndex> selected;
        for(auto ch:*charts) require(ch<ct.n_charts,"quality.charts","chart out of range");
        selected.insert(charts->begin(),charts->end());
        selected_charts_.assign(selected.begin(),selected.end());
    } else {selected_charts_.resize(ct.n_charts);std::iota(selected_charts_.begin(),selected_charts_.end(),0);}
    std::vector<MeshIndex> selected_faces,selected_triangles;
    for(auto ch:selected_charts_) {
        const auto faces=ct.chart_faces[ch],tris=geometry_->chart_triangles()[ch];
        selected_faces.insert(selected_faces.end(),faces.begin(),faces.end());selected_triangles.insert(selected_triangles.end(),tris.begin(),tris.end());
    }
    std::sort(selected_faces.begin(),selected_faces.end());std::sort(selected_triangles.begin(),selected_triangles.end());
    triangles_.reserve(selected_triangles.size());polygons_.reserve(selected_faces.size());
    for(auto f:selected_faces) {
        PolygonQuality p;p.face=f;p.chart=ct.chart_of_face[f];p.finite=true;p.exempt=exempt(f);p.degenerate_input=mesh.geometry().face_degenerate_kind[f]!=0;
        const auto range=mesh.face_corners(f);
        for(auto c:range) {const auto a=get(c),b=get(mesh.corner_next(c));p.signed_area+=a[0]*b[1]-b[0]*a[1];p.finite=p.finite && finite(a);}
        p.signed_area=p.finite?.5*p.signed_area:quality_nan;
        const auto sign=signs.empty()?1:signs[p.chart];p.raw_signed_area=p.signed_area*sign;
        const bool finite_area=std::isfinite(p.signed_area);
        p.flipped=finite_area && p.signed_area<=1e-12*mesh.geometry().face_area[f];
        p.raw_flipped=finite_area && p.raw_signed_area<=1e-12*mesh.geometry().face_area[f];
        if(range.size()==4) {
            const auto c=*range.begin();const auto a=get(c),b=get(c+1),d=get(c+2),e=get(c+3);
            const bool mixed0=(orient(a,b,d)<0)!=(orient(a,d,e)<0),mixed1=(orient(b,d,e)<0)!=(orient(b,e,a)<0);
            p.flipped=p.flipped || (finite_area && mixed0 && mixed1);
            auto raw_a=a,raw_b=b,raw_d=d,raw_e=e;
            for(auto pnt:{&raw_a,&raw_b,&raw_d,&raw_e}) (*pnt)[0]*=sign;
            const bool raw_mixed0=(orient(raw_a,raw_b,raw_d)<0)!=(orient(raw_a,raw_d,raw_e)<0);
            const bool raw_mixed1=(orient(raw_b,raw_d,raw_e)<0)!=(orient(raw_b,raw_e,raw_a)<0);
            p.raw_flipped=p.raw_flipped || (finite_area && raw_mixed0 && raw_mixed1);
        }
        polygons_.push_back(p);
    }
    std::vector<std::uint8_t> signs_of_tri(polygons_.size(),0);
    for(auto id:selected_triangles) {
        TriangleQuality t;t.corners=mt.corners[id];t.face=mt.faces[id];t.chart=ct.chart_of_face[t.face];t.exempt=exempt(t.face);
        const auto& frame=geometry_->frames()[id];t.area3=frame.double_area;t.sliver=t.area3<=geometry_->sliver_threshold();
        const auto p0=get(t.corners[0]),p1=get(t.corners[1]),p2=get(t.corners[2]);t.finite=finite(p0)&&finite(p1)&&finite(p2);
        const double u00=p1[0]-p0[0],u01=p2[0]-p0[0],u10=p1[1]-p0[1],u11=p2[1]-p0[1];
        t.area2=u00*u11-u01*u10;
        t.raw_area2=t.area2*(signs.empty()?1:signs[t.chart]);
        const auto& q=frame.coordinates;const double q00=q[1][0]-q[0][0],q01=q[2][0]-q[0][0],q10=q[1][1]-q[0][1],q11=q[2][1]-q[0][1];
        const double det=q00*q11-q01*q10;
        if(std::abs(det)>1e-16) {
            const double i00=q11/det,i01=-q01/det,i10=-q10/det,i11=q00/det;
            const double a=u00*i00+u01*i10,b=u00*i01+u01*i11,c=u10*i00+u11*i10,d=u10*i01+u11*i11;
            const double plus=std::sqrt(std::max((a+d)*(a+d)+(b-c)*(b-c),0.0));
            const double minus=std::sqrt(std::max((a-d)*(a-d)+(b+c)*(b+c),0.0));
            t.singular_values={.5*(plus+minus),.5*std::abs(plus-minus)};
        }
        const auto local=static_cast<std::size_t>(std::lower_bound(selected_faces.begin(),selected_faces.end(),t.face)-selected_faces.begin());
        if(t.area2<0) signs_of_tri[local]|=1;if(t.area2>0) signs_of_tri[local]|=2;
        triangles_.push_back(t);
    }
    for(std::size_t i=0;i<polygons_.size();++i) {
        auto& p=polygons_[i];const auto range=mesh.face_corners(p.face);
        if(range.size()>4 && signs_of_tri[i]==3 && std::isfinite(p.signed_area) && (!p.flipped || !p.raw_flipped)) {
            std::vector<Vec2> outline;outline.reserve(range.size());for(auto c:range) outline.push_back(get(c));
            const bool self_crosses=polygon_self_intersects(outline);p.flipped=p.flipped || self_crosses;p.raw_flipped=p.raw_flipped || self_crosses;
        }
    }
    chart_triangles_=cut_detail::rows(ct.n_charts,[&](auto add){for(std::size_t i=0;i<triangles_.size();++i) add(triangles_[i].chart,static_cast<MeshIndex>(i));});
    chart_polygons_=cut_detail::rows(ct.n_charts,[&](auto add){for(std::size_t i=0;i<polygons_.size();++i) add(polygons_[i].chart,static_cast<MeshIndex>(i));});
}
}
