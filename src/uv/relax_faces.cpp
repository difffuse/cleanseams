// SPDX-License-Identifier: GPL-3.0-or-later
#include "layout_postprocess_internal.hpp"
#include <algorithm>
namespace cleanseams::repair_detail {
std::optional<std::vector<Vec2>> relax_selected_faces(const PipelineInput& in,const CutSurface& cut,
    const std::function<bool()>& stopped,nlohmann::json& diagnostic) {
    const auto& ct=cut.tables();const auto& mesh=*in.mesh;
    diagnostic["reason"]="no_selected_faces";
    if(in.selected_faces.empty() || std::none_of(in.selected_faces.begin(),in.selected_faces.end(),[](auto x){return x!=0;}))return std::nullopt;
    if(stopped())return std::nullopt;
    std::vector<std::uint8_t> fixed(ct.n_sectors,0),moving(ct.n_sectors,0);
    std::vector<Vec2> position(ct.n_sectors);
    std::vector<std::pair<MeshIndex,MeshIndex>> edges;
    edges.reserve(in.uv->size());
    for(MeshIndex f=0;f<in.selected_faces.size();++f) {
        if(f%1024==0 && stopped())return std::nullopt;
        const auto corners=mesh.face_corners(f);
        for(std::size_t i=0;i<corners.size();++i) {
            const auto c=corners[static_cast<std::ptrdiff_t>(i)],a=ct.sector_of_corner[c],
                b=ct.sector_of_corner[corners[static_cast<std::ptrdiff_t>((i+1)%corners.size())]];
            position[a]=(*in.uv)[c];
            if(in.selected_faces[f])moving[a]=1;else fixed[a]=1;
            if(a!=b)edges.emplace_back(a,b); // repeated polygon edges retain their weight
        }
    }
    for(const auto& [c,p]:in.pins)fixed[ct.sector_of_corner[c]]=1;
    for(std::size_t s=0;s<moving.size();++s)moving[s]=moving[s] && !fixed[s];
    diagnostic["movable_sectors"]=std::count(moving.begin(),moving.end(),1);
    diagnostic["reason"]="no_movable_sectors";
    if(std::none_of(moving.begin(),moving.end(),[](auto x){return x!=0;}))return std::nullopt;
    std::vector<std::uint8_t> settled(moving.size());
    for(std::size_t s=0;s<moving.size();++s)settled[s]=!moving[s];
    std::vector<Vec2> summed(position.size());std::vector<double> total(position.size());
    for(int round=0;round<12;++round) {
        if(stopped())return std::nullopt;
        std::fill(summed.begin(),summed.end(),Vec2{0,0});std::fill(total.begin(),total.end(),0);
        // Accumulate both directions from the previous round, then update together.
        for(int reverse=0;reverse<2;++reverse)for(std::size_t e=0;e<edges.size();++e) {
            if(e%1024==0 && stopped())return std::nullopt;
            const auto a=reverse?edges[e].second:edges[e].first,b=reverse?edges[e].first:edges[e].second;
            const double weight=settled[b]?1:.25;total[a]+=weight;
            for(std::size_t d=0;d<2;++d)summed[a][d]+=weight*position[b][d];
        }
        for(std::size_t s=0;s<position.size();++s) {
            if(moving[s] && total[s]>0)for(std::size_t d=0;d<2;++d)position[s][d]=summed[s][d]/total[s];
            settled[s]=!moving[s] || total[s]>0;
        }
        diagnostic["rounds"]=round+1;
    }
    if(stopped())return std::nullopt;
    auto proposal=*in.uv;
    for(std::size_t c=0;c<proposal.size();++c)if(moving[ct.sector_of_corner[c]])proposal[c]=position[ct.sector_of_corner[c]];
    return proposal;
}
}
