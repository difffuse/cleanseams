// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "cleanseams/layout_postprocess.hpp"
namespace cleanseams {
// Internal Session entry: identical cache/worker/memory/timing scope to the public
// runner, with optional supplied-layout services. Never bypass the run envelope.
PipelineResult run_pipeline_with_services(const PipelineInput&,UnwrapBackend*,const ProgressSink&,const Control&,
    const LayoutPostprocessServices&);
}
