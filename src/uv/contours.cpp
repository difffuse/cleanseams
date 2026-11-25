// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/contours.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/control.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include "cleanseams/detail/numpy_compat.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <exception>
#include <new>
#include <optional>
#include <numbers>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);}
double degrees(double radians) {return radians*(180.0/std::numbers::pi);}
double norm2(const Vec2& a,const Vec2& b) {const double x=b[0]-a[0],y=b[1]-a[1];return std::sqrt(x*x+y*y);}
// Python float modulo for a positive divisor.
double python_mod(double x,double y) {
    if(y==0)throw ContractError(ErrorCode::numerical_failure,"contours.arc_anchors","float modulo by zero");
    double r=std::fmod(x,y);
    if(r!=0 && ((r<0)!=(y<0)))r+=y;
    return r;
}
std::size_t cyclic(std::int64_t x,std::size_t n) {const auto m=static_cast<std::int64_t>(n);return static_cast<std::size_t>(((x%m)+m)%m);}
std::vector<MeshIndex> chart_corners(const CutSurface& cut,MeshIndex chart) {
    std::vector<MeshIndex> ids;
    for(auto f:cut.tables().chart_faces[chart])for(auto c:cut.mesh().face_corners(f))ids.push_back(c);
    return ids;
}
std::vector<double> flat(std::span<const Vec2> v) {std::vector<double> out;for(const auto& p:v)out.insert(out.end(),p.begin(),p.end());return out;}
}
ContourMeasurement measure_chart(const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const Vec2> uv) {
    require(bool(geometry),"contours.geometry","missing geometry snapshot");
    const auto& cut=geometry->cut();require(chart<cut.tables().n_charts,"contours.chart","chart index out of range");
    require(uv.size()==cut.tables().sector_of_corner.size(),"contours.uv","UV array does not match the mesh corners");
    const auto ids=chart_corners(cut,chart);std::vector<Vec2> values;values.reserve(ids.size());
    for(auto c:ids)values.push_back(uv[c]);
    const std::array<MeshIndex,1> selected{chart};const CornerUvView view{values,ids};
    ContourMeasurement out;out.report=QualitySnapshot(geometry,view,{},std::span<const MeshIndex>(selected)).report();
    out.crossings=boundary_crossings(cut,view,{},std::span<const MeshIndex>(selected)).counts[chart];
    return out;
}
std::vector<std::int64_t> cut_degree(const CutSurface& cut) {
    const auto& m=cut.mesh();std::vector<std::int64_t> degree(m.input().positions.size(),0);
    const auto& ct=cut.tables();
    for(std::size_t e=0;e<ct.cut.size();++e)if(ct.cut[e])for(auto v:m.tables().edge_vertices[e])++degree[v];
    return degree;
}
std::vector<std::string> loop_anchor_flags(const CutSurface& cut,std::span<const MeshIndex> loop,std::span<const std::int64_t> degree,
                                           const std::set<MeshIndex>& pinned,const std::set<MeshIndex>& manual,const BoundaryOptions& o) {
    const auto& m=cut.mesh();const auto& ct=cut.tables();const auto& corner_vertices=m.input().corner_vertices;
    const auto& positions=m.input().positions;const auto n=loop.size();
    require(degree.size()==positions.size(),"contours.degree","cut degree does not match the mesh vertices");
    std::vector<std::string> flags(n);
    auto vertex=[&](MeshIndex c) {return static_cast<MeshIndex>(corner_vertices[c]);};
    for(std::size_t k=0;k<n;++k) {
        const auto c=loop[k];const auto v=vertex(c);
        if(manual.contains(v)) {flags[k]="manual";continue;}
        bool pin=pinned.contains(c);
        for(auto cc:m.tables().vert_corners[v])if(!pin && ct.sector_of_corner[cc]==ct.sector_of_corner[c] && pinned.contains(cc))pin=true;
        if(pin) {flags[k]="pin";continue;}
        if(degree[v]!=2) {flags[k]="node";continue;}
        // Geometric corner in 3D.
        const auto& p=positions[v];
        const auto a=detail::sub(positions[vertex(loop[(k+n-1)%n])],p),b=detail::sub(positions[vertex(loop[(k+1)%n])],p);
        const double la=detail::norm(a),lb=detail::norm(b);
        if(la>0 && lb>0) {
            const double c2=std::max(-1.0,std::min(1.0,detail::dot(a,b)/(la*lb)));
            if(180.0-degrees(std::acos(c2))>o.corner_angle)flags[k]="corner";
        }
    }
    std::vector<std::size_t> anchors;for(std::size_t k=0;k<n;++k)if(!flags[k].empty())anchors.push_back(k);
    if(anchors.size()>3) {
        std::vector<std::size_t> demoted;const auto count=anchors.size();
        for(std::size_t i=0;i<count;++i) {
            const auto a=static_cast<std::int64_t>(anchors[i]);
            const auto before=static_cast<std::int64_t>(anchors[(i+count-1)%count]),after=static_cast<std::int64_t>(anchors[(i+1)%count]);
            if(flags[anchors[i]]=="node" && cyclic(a-before,n)<2 && cyclic(after-a,n)<2)demoted.push_back(anchors[i]);
        }
        for(auto a:demoted)flags[a].clear();
    }
    return flags;
}
bool arc_is_loop_straight(const MeshSnapshot& m,std::span<const MeshIndex> verts,const BoundaryOptions& o) {
    const auto& t=m.tables();
    for(std::size_t k=1;k+1<verts.size();++k) {
        const auto a=verts[k-1],v=verts[k],b=verts[k+1];
        const auto e_in=m.edge_index(a,v),e_out=m.edge_index(v,b);
        require(e_in>=0 && e_out>=0,"contours.arc","arc step is not a mesh edge");
        if(t.vert_regular[v]) {
            if(t.fans[v]->straight(static_cast<MeshIndex>(e_in))!=e_out)return false;
        } else {
            const auto d_in=detail::mul(m.tangent_dir(v,static_cast<MeshIndex>(e_in)),-1.0),d_out=m.tangent_dir(v,static_cast<MeshIndex>(e_out));
            const double c=std::max(-1.0,std::min(1.0,detail::dot(d_in,d_out)));
            if(degrees(std::acos(c))>180.0-o.straight_angle)return false;
        }
    }
    return true;
}
bool arc_near_chord(std::span<const Vec2> p,double tol_edges) {
    std::vector<double> positive;
    for(std::size_t k=0;k+1<p.size();++k) {const double s=norm2(p[k],p[k+1]);if(s>0)positive.push_back(s);}
    const double h=positive.empty()?0.0:detail::numpy_median(positive);
    if(h<=0)return false;
    const auto& a=p.front();const auto& b=p.back();const Vec2 d{b[0]-a[0],b[1]-a[1]};
    const double length=std::sqrt(d[0]*d[0]+d[1]*d[1]);
    if(length<=0)return false;
    const Vec2 normal{-d[1]/length,d[0]/length};
    double worst=-INFINITY;
    for(const auto& q:p)worst=std::max(worst,std::abs((q[0]-a[0])*normal[0]+(q[1]-a[1])*normal[1]));
    return worst<=tol_edges*h;
}
bool chord_is_clear(std::span<const Vec2> p,const Vec2& a,const Vec2& b,std::span<const std::size_t> arc_positions) {
    const auto n=p.size();std::vector<std::uint8_t> in_arc(n,0);for(auto k:arc_positions)in_arc.at(k)=1;
    for(std::size_t k=0;k<n;++k)if(!in_arc[k] && segments_cross(a,b,p[k],p[(k+1)%n]))return false;
    return true;
}
std::optional<std::vector<Vec2>> straight_target(std::span<const Vec3> points,const Vec2& a,const Vec2& b) {
    std::vector<double> arc{0.0};
    for(std::size_t k=0;k+1<points.size();++k)arc.push_back(arc.back()+detail::norm(detail::sub(points[k+1],points[k])));
    if(arc.back()<=1e-12)return std::nullopt;
    std::vector<Vec2> out;
    for(double s:arc) {const double t=s/arc.back();out.push_back({(1-t)*a[0]+t*b[0],(1-t)*a[1]+t*b[1]});}
    return out;
}
std::optional<std::vector<Vec2>> smooth_target(std::span<const Vec2> p,double radius,double band) {
    const auto n=p.size();
    if(n<3)return std::vector<Vec2>(p.begin(),p.end());
    std::vector<double> seg,positive;
    for(std::size_t k=0;k+1<n;++k) {seg.push_back(norm2(p[k],p[k+1]));if(seg.back()>0)positive.push_back(seg.back());}
    const double h=positive.empty()?0.0:detail::numpy_median(positive);
    if(h<=0)return std::nullopt;
    const auto resampled=resample_polyline(p,h*0.5);const auto& res=resampled.points;
    auto smooth=gaussian_smooth_open(res,radius/0.5);
    const double limit=band*h;
    for(std::size_t i=0;i<res.size();++i) {
        Vec2 disp{smooth[i][0]-res[i][0],smooth[i][1]-res[i][1]};const double dn=std::sqrt(disp[0]*disp[0]+disp[1]*disp[1]);
        if(dn>limit) {const double scale=limit/dn;disp[0]*=scale;disp[1]*=scale;}
        smooth[i]={res[i][0]+disp[0],res[i][1]+disp[1]};
    }
    // Re-parametrise the smoothed curve by its own arc length; originals keep their
    // relative arc-length parameters.
    std::vector<double> arc2{0.0},arc_orig{0.0},xs,ys;
    for(std::size_t i=0;i+1<smooth.size();++i)arc2.push_back(arc2.back()+norm2(smooth[i],smooth[i+1]));
    for(double s:seg)arc_orig.push_back(arc_orig.back()+s);
    for(const auto& q:smooth) {xs.push_back(q[0]);ys.push_back(q[1]);}
    std::vector<Vec2> out;
    for(double s:arc_orig) {
        const double tt=s/arc_orig.back()*arc2.back();
        out.push_back({detail::numpy_interp(tt,arc2,xs),detail::numpy_interp(tt,arc2,ys)});
    }
    out.front()=p.front();out.back()=p.back();
    return out;
}
InteriorSolver default_interior_solver(const BoundaryOptions& o,const InteriorSolverSettings& settings) {
    // A negative count runs no ARAP iteration, as the source's range() (flatten.arap).
    const auto iterations=static_cast<std::uint64_t>(std::max<std::int64_t>(0,o.arap_iterations));
    return [iterations,settings](const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const SectorPin> pins,
                                 std::span<const Vec2>) {
        ChartSolveRequest r;r.geometry=geometry;r.chart=chart;r.sector_pins.assign(pins.begin(),pins.end());r.iterations=iterations;
        r.policy=settings.policy;r.validate=false;r.memory_budget=settings.memory_budget;r.pin_completion=settings.pin_completion;
        auto solved=solve_arap(r);
        throw_if_cancelled(settings.cancel);
        return std::pair{std::move(solved.corner_ids),std::move(solved.corner_uv)};
    };
}
ContourValidator default_validator(const BoundaryOptions& o) {
    return [o](const std::shared_ptr<const QualityGeometry>& geometry,MeshIndex chart,std::span<const Vec2> uv,const ContourMeasurement& baseline) {
        auto m=measure_chart(geometry,chart,uv);const auto& r=m.report;const auto& b=baseline.report;
        const bool ok=!(r.nonfinite || r.flipped || m.crossings || r.collapsed_uv) &&
            !(r.anisotropy_p95>std::min(b.anisotropy_p95*o.p95_growth,o.abs_anisotropy_p95)) &&
            !(r.area_p95>std::min(b.area_p95*o.p95_growth,o.abs_area_p95)) &&
            !(r.anisotropy_max>b.anisotropy_max*o.max_growth);
        return std::pair{ok,std::move(m)};
    };
}
namespace {
ContourOutcome clean_impl(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,MeshIndex chart,const ContourRequest& q) {
    require(bool(geometry),"contours.geometry","missing geometry snapshot");
    const auto& cut=geometry->cut();const auto& ct=cut.tables();const auto& m=cut.mesh();const auto& o=q.options;
    require(chart<ct.n_charts,"contours.chart","chart index out of range");
    require(uv.size()==ct.sector_of_corner.size(),"contours.uv","UV array does not match the mesh corners");
    for(auto c:q.pinned_corners)require(c<ct.sector_of_corner.size(),"contours.pins","pinned corner out of range");
    // "off" is resolved by the options (clean_contours). Any mode other than auto and panels
    // smooths every arc, as the source's clean_chart does.
    ContourOutcome out;out.uv.assign(uv.begin(),uv.end());auto& res=out.result;res.chart=chart;
    auto finish=[&](const char* status,const char* reason) {res.status=status;res.reason=reason;return std::move(out);};
    if(q.protect)return finish("protected","protected_chart");
    if(!ct.chart_is_disk[chart])return finish("skipped","not_a_disk");
    if(ct.chart_loops[chart].size()!=1)return finish("skipped","loop_count");
    const auto loop=ct.loops[ct.chart_loops[chart][0]];const auto n=loop.size();
    if(n<4)return finish("skipped","too_small");
    const auto baseline=measure_chart(geometry,chart,uv);res.before=baseline;
    if(baseline.report.nonfinite || baseline.report.flipped || baseline.crossings)return finish("skipped","invalid_input_map");
    const auto default_check=default_validator(o);const auto& validate=q.validate?*q.validate:default_check;
    const auto& corner_vertices=m.input().corner_vertices;
    auto vertex=[&](MeshIndex c) {return static_cast<MeshIndex>(corner_vertices[c]);};
    const bool has_pins=std::any_of(q.pinned_corners.begin(),q.pinned_corners.end(),[&](auto c) {return ct.chart_of_corner[c]==chart;});
    const auto corners=chart_corners(cut,chart);
    const bool has_manual=std::any_of(corners.begin(),corners.end(),[&](auto c) {return q.manual_anchors.contains(vertex(c));});
    // A topological lattice admits an exact rectangle, recognized once and still validated.
    if((o.mode=="auto" || o.mode=="panels") && !has_pins && !has_manual) {
        GridStats stats;const auto lattice=grid_topology(cut,chart,&stats);res.grid_recognitions=stats.recognitions;
        if(lattice) {
            bool straight=true;
            for(const auto& arc:lattice->arcs) {
                std::vector<MeshIndex> verts;for(auto c:arc)verts.push_back(vertex(c));
                if(!arc_is_loop_straight(m,verts,o)) {straight=false;break;}
            }
            if(straight)if(auto map=grid_uv(*lattice,cut)) {
                std::vector<Vec2> proposed(uv.begin(),uv.end());
                for(std::size_t i=0;i<map->corner_ids.size();++i)proposed[map->corner_ids[i]]=map->uv[i];
                auto [ok,measured]=validate(geometry,chart,proposed,baseline);
                if(ok) {
                    res.intensity=1.0;res.after=std::move(measured);
                    for(std::size_t i=0;i<4;++i) {
                        std::vector<Vec2> points;for(auto c:lattice->arcs[i])points.push_back(uv[c]);
                        ContourArc arc;arc.index=static_cast<std::int64_t>(i);arc.mode="straight";arc.n=lattice->arcs[i].size();
                        arc.jag_before=arc_jaggedness(points,std::nullopt,true).rms;arc.jag_after=0.0;res.arcs.push_back(std::move(arc));
                    }
                    res.grid=std::move(*map);out.uv=std::move(proposed);
                    return finish("optimized","grid_rectangle");
                }
            }
        }
    }
    const auto degree=cut_degree(cut);
    auto flags=loop_anchor_flags(cut,loop,degree,q.pinned_corners,q.manual_anchors,o);
    std::vector<Vec2> p;std::vector<MeshIndex> verts;std::vector<Vec3> p3;
    for(auto c:loop) {p.push_back(uv[c]);verts.push_back(vertex(c));p3.push_back(m.input().positions[vertex(c)]);}
    std::vector<std::size_t> anchors;for(std::size_t k=0;k<n;++k)if(!flags[k].empty())anchors.push_back(k);
    if(anchors.size()<3) {
        // Closed loop with too few anchors: add some by arc length, keeping the existing ones.
        std::vector<double> arc{0.0};for(std::size_t k=0;k<n;++k)arc.push_back(arc.back()+norm2(p[k],p[(k+1)%n]));
        const auto want=std::max<std::int64_t>(o.min_closed_anchors,3);const auto start=anchors.empty()?0:anchors[0];
        for(std::int64_t j=0;j<want;++j) {
            const double value=python_mod(arc[start]+arc[n]*static_cast<double>(j)/static_cast<double>(want),arc[n]);
            std::size_t best=0;for(std::size_t k=1;k<n;++k)if(std::abs(arc[k]-value)<std::abs(arc[best]-value))best=k;
            if(flags[best].empty())flags[best]="arc";
        }
        anchors.clear();for(std::size_t k=0;k<n;++k)if(!flags[k].empty())anchors.push_back(k);
    }
    res.flags=flags;
    auto target=p;const auto count=anchors.size();
    for(std::size_t i=0;i<count;++i) {
        const auto a=anchors[i],b=anchors[(i+1)%count];
        std::vector<std::size_t> positions{a};for(auto k=a;k!=b;) {k=(k+1)%n;positions.push_back(k);}
        if(positions.size()<3)continue;
        std::vector<Vec2> puv;std::vector<Vec3> pxyz;std::vector<MeshIndex> vlist;
        for(auto k:positions) {puv.push_back(p[k]);pxyz.push_back(p3[k]);vlist.push_back(verts[k]);}
        std::string mode=o.mode;
        if(mode=="auto") {
            const bool straight=arc_is_loop_straight(m,vlist,o) && arc_near_chord(puv,o.straight_tolerance) && chord_is_clear(p,p[a],p[b],positions);
            mode=straight?"straight":"smooth";
        } else if(mode=="panels")mode=chord_is_clear(p,p[a],p[b],positions)?"straight":"smooth";
        else mode="smooth";
        ContourArc arc;arc.index=static_cast<std::int64_t>(i);arc.n=positions.size();arc.positions=positions;
        arc.jag_before=arc_jaggedness(puv,std::nullopt,mode=="straight").rms;
        auto t=mode=="straight"?straight_target(pxyz,p[a],p[b]):smooth_target(puv,o.smooth_radius,o.max_displacement);
        if(!t) {arc.mode="degenerate";res.arcs.push_back(std::move(arc));continue;}
        arc.target=*t;
        bool pinned=false;
        for(std::size_t j=1;j+1<positions.size();++j)if(flags[positions[j]]=="pin")pinned=true;
        if(pinned) {arc.mode="pinned";res.arcs.push_back(std::move(arc));continue;}
        for(std::size_t j=0;j<positions.size();++j)target[positions[j]]=(*t)[j];
        arc.mode=mode;res.arcs.push_back(std::move(arc));
    }
    if(detail::numpy_allclose(flat(target),flat(p)))return finish("kept","already_clean");
    const auto default_solver=q.solve?InteriorSolver():default_interior_solver(o);const auto& solve=q.solve?*q.solve:default_solver;
    // PIN-02: user pins on sectors off the contour keep their incoming UV (lowest corner first).
    std::map<MeshIndex,std::pair<MeshIndex,Vec2>> held;
    if(q.pin_policy==PinCompletion::preserve) {
        std::set<MeshIndex> loop_sectors;for(auto c:loop)loop_sectors.insert(ct.sector_of_corner[c]);
        for(auto c:q.pinned_corners) {
            const auto sector=ct.sector_of_corner[c];
            if(ct.chart_of_corner[c]==chart && !loop_sectors.contains(sector) && !held.contains(sector))held[sector]={c,uv[c]};
        }
    }
    for(double s:o.intensities) {
        std::vector<Vec2> candidate(n);std::vector<SectorPin> pins;
        for(std::size_t k=0;k<n;++k) {
            candidate[k]={p[k][0]+s*(target[k][0]-p[k][0]),p[k][1]+s*(target[k][1]-p[k][1])};
            pins.push_back({ct.sector_of_corner[loop[k]],candidate[k]});
        }
        for(const auto& [sector,pin]:held)pins.push_back({sector,pin.second});
        const auto [ids,values]=solve(geometry,chart,pins,uv);++res.interior_solves;
        require(ids.size()==values.size(),"contours.solver","interior solver returned mismatched arrays");
        std::vector<Vec2> updated(uv.begin(),uv.end());
        for(std::size_t i=0;i<ids.size();++i)updated.at(ids[i])=values[i];
        // The pinned sectors land exactly on their pins whatever the solver did.
        for(std::size_t k=0;k<n;++k) {
            const auto sector=ct.sector_of_corner[loop[k]];
            for(auto c:m.tables().vert_corners[verts[k]])if(ct.sector_of_corner[c]==sector)updated[c]=candidate[k];
        }
        for(const auto& [sector,pin]:held) {
            for(auto c:m.tables().vert_corners[vertex(pin.first)])if(ct.sector_of_corner[c]==sector)updated[c]=pin.second;
        }
        auto [ok,measured]=validate(geometry,chart,updated,baseline);
        if(ok) {
            res.intensity=s;res.after=std::move(measured);
            std::vector<Vec2> pn;for(auto c:loop)pn.push_back(updated[c]);
            for(auto& arc:res.arcs) {
                if(arc.positions.empty())continue;
                std::vector<Vec2> points;for(auto k:arc.positions)points.push_back(pn[k]);
                arc.jag_after=arc_jaggedness(points,std::nullopt,arc.mode=="straight").rms;
            }
            out.uv=std::move(updated);res.status="optimized";
            return out;
        }
        res.reason="intensity_"+detail::python_float_repr(s)+"_failed";
    }
    return finish("kept","deformation_budget");
}
}
ContourOutcome clean_chart(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,MeshIndex chart,const ContourRequest& q) {
    try {return clean_impl(geometry,uv,chart,q);}
    catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,"contours","contour cleanup allocation failed");}
}
CleanAllOutcome clean_all(const std::shared_ptr<const QualityGeometry>& geometry,std::span<const Vec2> uv,const ContourRequest& q,
                          std::optional<std::span<const MeshIndex>> charts,const Deadline& deadline,
                          const std::function<void(const nlohmann::json&)>& progress) {
    require(bool(geometry),"contours.geometry","missing geometry snapshot");
    const auto& ct=geometry->cut().tables();
    CleanAllOutcome out;out.uv.assign(uv.begin(),uv.end());
    std::vector<MeshIndex> selected;
    if(charts)selected.assign(charts->begin(),charts->end());else for(MeshIndex ch=0;ch<ct.n_charts;++ch)selected.push_back(ch);
    std::size_t first_serial=0;
    if(q.parallel && !deadline.at() && selected.size()>1 && active_workers() && !in_worker_task()) {
        // GB performance: a chart's cleanup reads and writes only its own corners, so every
        // chart (up to the first out-of-range index, where the serial loop stops) is cleaned
        // in parallel from the incoming map and committed in chart order, with the serial
        // loop's progress events and first failure. A result that changed another chart's
        // corners (never expected) is discarded and the serial loop resumes at that chart.
        std::size_t valid=0;while(valid<selected.size() && selected[valid]<ct.n_charts)++valid;
        std::vector<std::size_t> work;for(std::size_t i=0;i<valid;++i)if(!ct.chart_faces[selected[i]].empty())work.push_back(i);
        std::vector<std::optional<ContourOutcome>> slots(work.size());std::vector<std::exception_ptr> errors(work.size());
        parallel_for(work.size(),[&](std::size_t k) {
            try {slots[k]=clean_chart(geometry,uv,selected[work[k]],q);} catch(...) {errors[k]=std::current_exception();}
        });
        first_serial=valid;
        for(std::size_t k=0;k<work.size();++k) {
            if(errors[k])std::rethrow_exception(errors[k]);
            const auto ch=selected[work[k]];auto& one=*slots[k];
            bool local=one.uv.size()==uv.size();
            for(std::size_t c=0;c<uv.size() && local;++c)local=ct.chart_of_corner[c]==ch || one.uv[c]==uv[c];
            if(!local) {first_serial=work[k];break;}
            for(std::size_t c=0;c<uv.size();++c)if(ct.chart_of_corner[c]==ch)out.uv[c]=one.uv[c];
            out.charts.push_back(std::move(one.result));
            if(progress)progress({{"stage","contours"},{"status","running"},{"done",out.charts.size()},{"total",ct.n_charts}});
        }
    }
    for(std::size_t i=first_serial;i<selected.size();++i) {
        const auto ch=selected[i];
        if(deadline.expired())break;
        require(ch<ct.n_charts,"contours.charts","chart index out of range");
        if(ct.chart_faces[ch].empty())continue;
        auto one=clean_chart(geometry,out.uv,ch,q);
        out.uv=std::move(one.uv);out.charts.push_back(std::move(one.result));
        if(progress)progress({{"stage","contours"},{"status","running"},{"done",out.charts.size()},{"total",ct.n_charts}});
    }
    return out;
}
}
