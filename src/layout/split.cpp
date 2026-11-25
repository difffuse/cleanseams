// SPDX-License-Identifier: GPL-3.0-or-later
// Phase B (T32): records of the baseline's broad over-budget split, so that enhanced_v1
// can evaluate relief slits before it (Sections 14.3 and 21). Nothing here changes a split
// or any legacy_010/native_baseline decision; the records only describe what the split did.
#include "cleanseams/layout.hpp"
#include "cleanseams/errors.hpp"
#include <algorithm>
#include <cmath>
#include <set>
namespace cleanseams {
nlohmann::json SplitRecord::summary() const {
    return {{"stage",stage},{"chart",chart},{"faces",faces.size()},{"added_edges",added_edges},
            {"anisotropy_p95",std::isfinite(anisotropy_p95)?nlohmann::json(anisotropy_p95):nlohmann::json(nullptr)},
            {"chart_key",surface?surface->chart_key(chart).digest():std::string()}};
}
std::vector<SplitRecord> record_split(const std::string& stage,const std::shared_ptr<const CutSurface>& before,std::span<const Vec2> uv,
                                      const CutSurface& after,std::span<const MeshIndex> over_charts,std::span<const double> anisotropy_p95) {
    if(!before)throw ContractError(ErrorCode::invalid_input,"split.before","missing surface before the split");
    const auto& m=before->mesh();const auto& t=m.tables();const auto& bt=before->tables();const auto& at=after.tables();
    if(!after.mesh().same_geometry(m))throw ContractError(ErrorCode::invalid_input,"split.after","surfaces of different meshes");
    if(uv.size()!=bt.sector_of_corner.size())throw ContractError(ErrorCode::invalid_input,"split.uv","one UV per corner required");
    if(!anisotropy_p95.empty() && anisotropy_p95.size()!=over_charts.size())
        throw ContractError(ErrorCode::invalid_input,"split.anisotropy_p95","one measurement per over-budget chart");
    std::vector<SplitRecord> out;
    for(std::size_t i=0;i<over_charts.size();++i) {
        const auto ch=over_charts[i];
        if(ch>=bt.n_charts)throw ContractError(ErrorCode::invalid_input,"split.over_charts","chart out of range");
        SplitRecord r;r.stage=stage;r.surface=before;r.uv.assign(uv.begin(),uv.end());r.chart=ch;
        r.faces.assign(bt.chart_faces[ch].begin(),bt.chart_faces[ch].end());std::sort(r.faces.begin(),r.faces.end());
        if(!anisotropy_p95.empty())r.anisotropy_p95=anisotropy_p95[i];
        std::set<MeshIndex> added;
        for(auto f:r.faces)
            for(auto c:m.face_corners(f)) {
                const auto e=t.corner_edge[c];
                if(at.cut[e] && !bt.cut[e])added.insert(e);
            }
        if(added.empty())continue;
        r.added_edges.assign(added.begin(),added.end());
        out.push_back(std::move(r));
    }
    return out;
}
std::vector<MeshIndex> split_edges_present(const SplitRecord& record,std::span<const std::uint8_t> cuts) {
    std::vector<MeshIndex> out;
    for(auto e:record.added_edges)if(e<cuts.size() && cuts[e])out.push_back(e);
    return out;
}
}
