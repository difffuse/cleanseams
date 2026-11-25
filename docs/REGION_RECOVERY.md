# Native mixed-quad recovery (development candidate)

With `operation="pipeline"`, `options.policy="repair_v1"` and a positive
`repair.region_recovery_seconds`, a failed native/host chart solve can use the new
mixed-quad fallback before emergency per-face projection. `legacy_010`,
`native_baseline` and `enhanced_v1` retain their existing behavior. A zero allowance
skips this optional fallback. The existing triangle fallback remains separate.

The implementation in `src/uv/region_recovery.cpp` performs local fault isolation
without importing a Python engine. Eligibility requires 2–256 faces, at most 1,024 corners, one cut-sector
boundary loop, Euler characteristic one, no pinched cut sector and no pinned corners.
An unsupported original vertex fan is still reported. It does not itself invalidate
an otherwise usable disk of cut sectors. This does not redefine the existing global
`chart_is_disk` flag or claim that the original mesh's unsupported fan was repaired.

The native LSCM proposal must contain every original corner exactly once in the expected
order, with finite coordinates. Fixed source triangles, boundary validity and strict
positive-area overlap identify damaged faces. At most `max(3, face_count / 10)` local
faults may be isolated using a convex, density-preserving face projection. Original
source geometry, pins, existing cuts and unrelated/cached UVs remain untouched. A
candidate that needs to change a locked cut is rejected; an emergency that would cut a
locked continuous edge fails explicitly rather than silently violating the lock.

The entire local recovery phase receives at most one second and eight attempts, further
bounded by the caller's allowance and operation deadline. Earlier successful native
solves do not consume this local allowance. Failed attempts share the remainder; it is
not renewed per chart. Native clocks and local solver services can be supplied through
`RepairContext` for deterministic deadline/error testing. They cannot bypass proposal
validation. Cancellation and memory-limit failures propagate.

Before adoption, the repaired region is converted to the requested storage precision
and checked again. Its P95 anisotropy ratio must be at most 2.0, equivalent to
stretch at most 1.0 (ratio minus one). This uses the same step-CDF percentile and original **face** area
repeated for each fixed triangle, not the engine's interpolated triangle-area statistic.
The report gives both units and identifies the weighting and `mesh_snapshot_fixed`
triangulation. Native ear clipping can differ from convex n-gon fans; callers
using another triangulation must independently check their returned layout.
Identical metric values are not guaranteed. A successful solve alone does not
establish the requested fragmentation or quality target.

`deadline_chart_fallback.local_records` records each attempted chart, topology,
acceptance/rejection reason, damaged original-face IDs and final distortion when
available. Solver rejection is reported as `unwrap.status="fallback"` under
`repair_v1`, separately from actual operation timeout. The legacy policies' report
semantics remain unchanged. `build_info().pipeline_recovery_operations` advertises
`mixed_quad_fallback`; the supplied-layout `repair_operations` list now includes
`relax_faces`, `repair_regions` and `join_small_islands` (see [the separate supplied-layout contract](LAYOUT_POSTPROCESS.md)).

UV recovery does not repair unsupported input topology. Consumers must distinguish
`uv_valid` from `topology_valid` and honor whole-result acceptance.
