// SPDX-License-Identifier: GPL-3.0-or-later
#include "cleanseams/build_info.hpp"
#include "cleanseams/options.hpp"
#include "cleanseams/control.hpp"
#include "cs_build_info.hpp"
namespace cleanseams {
nlohmann::json build_info() {
    using Json=nlohmann::json;
    return {{"engine","cleanseams"},{"version",CS_ENGINE_VERSION},
        {"source_revision",CS_SOURCE_REVISION},{"source_tree_sha256",CS_SOURCE_TREE_SHA256},
        {"schemas",{{"case",1},{"run_input",1},{"result",1},{"options",1},{"array_descriptor",1}}},
        {"compiler",CS_COMPILER},{"dependencies",Json::parse(CS_DEPENDENCIES_JSON)},{"precision","float64"},
        {"backends",{"native_arap"}},{"behaviors",{"legacy_010","native_baseline","enhanced_v1","repair_v1"}},
        {"default_behavior",PipelineOptions{}.policy},{"default_threads",default_native_threads()},
        {"operations",{"pipeline","postprocess_layout"}},{"host_capabilities",{"solve","pack","region_proposal"}},
        {"repair_operations",{"relax_faces","repair_regions","join_small_islands","repair_singletons","recover_fragments"}},{"join_modes",{"small","fewer"}},{"pipeline_recovery_operations",{"mixed_quad_fallback"}},
        {"features",{{"mesh_analysis",true},{"pipeline",true},{"case_arrays",{"float64-le","int64-le","uint64-le","uint8"}},
            {"enhanced_v1",true},{"case_constraints",true},{"python_engine",false},{"supplied_layout_validation",true},{"host_method_selection",true}}}};
}
}
