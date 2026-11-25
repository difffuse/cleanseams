// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/grid.hpp"
#include "cleanseams/errors.hpp"
#include "cleanseams/detail/mesh_math.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <new>
namespace cleanseams {
namespace {
void require(bool yes,const char* field,const char* message) {
    if(!yes)throw ContractError(ErrorCode::invalid_input,field,message);
}
std::optional<GridDescriptor> recognize(const CutSurface& cut,MeshIndex chart) {
    const auto& ct=cut.tables();const auto& mesh=cut.mesh();const auto& offsets=mesh.input().face_offsets;
    require(chart<ct.n_charts,"grid.chart","chart index out of range");
    if(!ct.chart_is_disk[chart] || ct.chart_loops[chart].size()!=1)return std::nullopt;
    const auto faces=ct.chart_faces[chart];
    if(faces.empty() || std::any_of(faces.begin(),faces.end(),[&](auto f) {return offsets[f+1]-offsets[f]!=4;}))return std::nullopt;
    GridDescriptor g;g.chart=chart;g.faces.assign(faces.begin(),faces.end());
    std::vector<std::uint64_t> count(ct.n_sectors,0);
    for(auto f:faces) {
        std::array<MeshIndex,4> sec{};
        for(std::size_t k=0;k<4;++k) {
            const auto c=static_cast<MeshIndex>(offsets[f]+k);g.corners.push_back(c);sec[k]=ct.sector_of_corner[c];++count[sec[k]];
        }
        g.face_sectors.push_back(sec);
    }
    std::vector<MeshIndex> loop(ct.loops[ct.chart_loops[chart][0]].begin(),ct.loops[ct.chart_loops[chart][0]].end());
    const auto n=loop.size();
    std::vector<MeshIndex> boundary_sectors;for(auto c:loop)boundary_sectors.push_back(ct.sector_of_corner[c]);
    auto unique=boundary_sectors;std::sort(unique.begin(),unique.end());
    if(std::unique(unique.begin(),unique.end())!=unique.end())return std::nullopt;
    std::vector<std::uint8_t> boundary(ct.n_sectors,0);for(auto s:boundary_sectors)boundary[s]=1;
    std::uint64_t present=0;
    for(MeshIndex s=0;s<ct.n_sectors;++s) {
        if(!count[s])continue;++present;
        if(!boundary[s] && count[s]!=4)return std::nullopt;
    }
    for(auto s:boundary_sectors)if(count[s]!=1 && count[s]!=2)return std::nullopt;
    auto anchors=[&] {std::vector<std::size_t> a;for(std::size_t k=0;k<n;++k)if(count[ct.sector_of_corner[loop[k]]]==1)a.push_back(k);return a;};
    if(anchors().size()!=4)return std::nullopt;
    std::rotate(loop.begin(),loop.begin()+static_cast<std::ptrdiff_t>(anchors()[0]),loop.end());
    for(std::size_t k=0;k<n;++k)boundary_sectors[k]=ct.sector_of_corner[loop[k]];
    const auto at=anchors();
    std::array<std::uint64_t,4> sizes{};for(std::size_t k=0;k<4;++k)sizes[k]=(k+1<4?at[k+1]:n)-at[k];
    if(sizes[0]!=sizes[2] || sizes[1]!=sizes[3])return std::nullopt;
    const auto nu=sizes[0],nv=sizes[1];
    if((nu+1)*(nv+1)!=present)return std::nullopt;
    std::vector<GridCoordinate> ij(ct.n_sectors,GridCoordinate{-1,-1});
    std::size_t k=0;const auto su=static_cast<std::int64_t>(nu),sv=static_cast<std::int64_t>(nv);
    for(std::int64_t i=0;i<su;++i)ij[boundary_sectors[k++]]={i,0};
    for(std::int64_t j=0;j<sv;++j)ij[boundary_sectors[k++]]={su,j};
    for(std::int64_t i=0;i<su;++i)ij[boundary_sectors[k++]]={su-i,sv};
    for(std::int64_t j=0;j<sv;++j)ij[boundary_sectors[k++]]={0,sv-j};
    // Source breadth-first completion: a face with one unknown corner closes its parallelogram.
    std::vector<std::vector<std::size_t>> incident(ct.n_sectors);
    for(std::size_t f=0;f<g.face_sectors.size();++f)for(auto s:g.face_sectors[f])incident[s].push_back(f);
    std::deque<std::size_t> queue;for(std::size_t f=0;f<g.face_sectors.size();++f)queue.push_back(f);
    while(!queue.empty()) {
        const auto f=queue.front();queue.pop_front();const auto& sec=g.face_sectors[f];
        std::size_t missing=0,which=0;for(std::size_t m=0;m<4;++m)if(ij[sec[m]][0]<0) {++missing;which=m;}
        if(missing!=1)continue;
        const auto& a=ij[sec[(which+3)%4]];const auto& b=ij[sec[(which+1)%4]];const auto& c=ij[sec[(which+2)%4]];
        ij[sec[which]]={a[0]+b[0]-c[0],a[1]+b[1]-c[1]};
        for(auto face:incident[sec[which]])queue.push_back(face);
    }
    std::vector<GridCoordinate> grid;
    for(MeshIndex s=0;s<ct.n_sectors;++s)if(count[s]) {
        const auto& p=ij[s];
        if(p[0]<0 || p[1]<0 || p[0]>su || p[1]>sv)return std::nullopt;
        g.sectors.push_back(s);g.ij.push_back(p);grid.push_back(p);
    }
    std::sort(grid.begin(),grid.end());if(std::unique(grid.begin(),grid.end())!=grid.end())return std::nullopt;
    for(const auto& sec:g.face_sectors) {
        std::array<GridCoordinate,4> steps{};
        for(std::size_t m=0;m<4;++m) {
            const auto& p=ij[sec[m]];const auto& q=ij[sec[(m+1)%4]];steps[m]={q[0]-p[0],q[1]-p[1]};
            if(std::abs(steps[m][0])+std::abs(steps[m][1])!=1)return std::nullopt;
        }
        if(steps[0][0]*steps[1][1]-steps[0][1]*steps[1][0]!=1)return std::nullopt;
    }
    g.nu=nu;g.nv=nv;
    for(std::size_t side=0;side<4;++side) {
        const auto a=at[side],b=side+1<4?at[side+1]:n;
        g.arcs[side].assign(loop.begin()+static_cast<std::ptrdiff_t>(a),loop.begin()+static_cast<std::ptrdiff_t>(b));
        g.arcs[side].push_back(loop[b%n]);
    }
    g.chart_key=cut.chart_key(chart).digest();
    return g;
}
}
GridCoordinate GridDescriptor::coordinate(MeshIndex sector) const {
    const auto found=std::lower_bound(sectors.begin(),sectors.end(),sector);
    require(found!=sectors.end() && *found==sector,"grid.sector","sector is not part of the lattice");
    return ij[static_cast<std::size_t>(found-sectors.begin())];
}
std::optional<GridDescriptor> grid_topology(const CutSurface& cut,MeshIndex chart,GridStats* stats) {
    try {
        if(stats)++stats->recognitions;
        return recognize(cut,chart);
    } catch(const std::bad_alloc&) {throw ContractError(ErrorCode::memory_limit,"grid.topology","grid recognition allocation failed");}
}
bool is_grid_chart(const CutSurface& cut,MeshIndex chart) {return grid_topology(cut,chart).has_value();}
std::optional<GridMap> grid_uv(const GridDescriptor& g,const CutSurface& cut,std::span<const CornerPin> pins) {
    const auto& ct=cut.tables();const auto& mesh=cut.mesh();
    require(g.chart<ct.n_charts && cut.chart_key(g.chart).digest()==g.chart_key,"grid.descriptor","descriptor belongs to another chart or cut");
    require(g.nu>0 && g.nv>0,"grid.descriptor","empty lattice");
    std::array<std::vector<double>,4> arc_t;std::array<double,4> lengths{};
    for(std::size_t side=0;side<4;++side) {
        const auto& arc=g.arcs[side];std::vector<double> cumulative{0.0};
        for(std::size_t k=1;k<arc.size();++k) {
            const auto& p=mesh.input().positions[static_cast<std::size_t>(mesh.input().corner_vertices[arc[k-1]])];
            const auto& q=mesh.input().positions[static_cast<std::size_t>(mesh.input().corner_vertices[arc[k]])];
            cumulative.push_back(cumulative.back()+detail::norm(detail::sub(q,p)));
        }
        lengths[side]=cumulative.back();
        if(!(lengths[side]>0))return std::nullopt;
        for(auto& c:cumulative)c/=lengths[side];
        arc_t[side]=std::move(cumulative);
    }
    GridMap out;out.lu=(lengths[0]+lengths[2])/2;out.lv=(lengths[1]+lengths[3])/2;out.grid={g.nu,g.nv};
    out.corner_ids=g.corners;out.uv.reserve(g.corners.size());
    const auto nu=static_cast<std::int64_t>(g.nu),nv=static_cast<std::int64_t>(g.nv);
    for(auto c:g.corners) {
        const auto [i,j]=g.coordinate(ct.sector_of_corner[c]);
        const double s=double(i)/double(nu),t=double(j)/double(nv);
        const auto at=[](const std::vector<double>& v,std::int64_t k) {return v[static_cast<std::size_t>(k)];};
        const double x=out.lu*((1-t)*at(arc_t[0],i)+t*(1-at(arc_t[2],nu-i)));
        const double y=out.lv*((1-s)*(1-at(arc_t[3],nv-j))+s*at(arc_t[1],j));
        out.uv.push_back({x,y});
    }
    for(const auto& pin:pins) {
        require(pin.corner>=0 && static_cast<std::uint64_t>(pin.corner)<ct.chart_of_corner.size(),"grid.pins","corner index out of range");
        if(ct.chart_of_corner[static_cast<std::size_t>(pin.corner)]!=g.chart)continue;
        const auto found=std::lower_bound(out.corner_ids.begin(),out.corner_ids.end(),static_cast<MeshIndex>(pin.corner));
        if(found==out.corner_ids.end() || *found!=static_cast<MeshIndex>(pin.corner) ||
           out.uv[static_cast<std::size_t>(found-out.corner_ids.begin())]!=pin.uv)
            throw ContractError(ErrorCode::constraint_conflict,"grid.pins","the grid map would move a pinned corner");
    }
    return out;
}
}
