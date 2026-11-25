# Supplied-layout processing (development candidate)

`operation="postprocess_layout"` and `options.policy="repair_v1"` submit an existing
per-corner UV layout to the native Session or direct pipeline API. Omitting operation
retains the existing pipeline behavior. This operation skips seam generation, path
refinement and global unwrap, even with a host backend. Missing UVs, invalid masks and
unknown operations are rejected before Session worker dispatch. Nonfinite coordinates
are allowed only on explicitly selected faces with a positive region recovery allowance.

This candidate implements supplied-layout validation, selected-face relaxation and
connected-region recovery, protected stored-cut recovery, small/fewer chart joining and bounded singleton reconstruction.
Native fragment selection and acceptance are available through `LayoutPostprocessServices`
and the opt-in Session `region_proposal` capability. Clients without that capability
report `recover_fragments.reason="host_capability_unavailable"` when the allowance is
positive. Other native recovery phases still run. Nonzero host-proposal allowances are supported
only for `postprocess_layout` with `repair_v1`. The separate
[pipeline mixed-quad fallback](REGION_RECOVERY.md) supports a region recovery allowance.

`RepairOptions` contains `relax_faces_seconds`, `region_recovery_seconds`,
`singleton_seconds`, `join_seconds`, `host_proposal_seconds` (finite nonnegative total
phase allowances, default zero), and `storage_precision` (`float64` by default, or
`float32`). Sub-budgets are bounded by the remaining operation deadline. No repair
allowance is implicitly renewed per chart. `selected_faces` is an optional mask in
original face order; an empty mask means no explicit selection. The Python binding
requires a contiguous bool array for that mask; native/wire masks contain zeros and ones.
`minimum_relative_area` is a finite nonnegative fraction (default zero) of the input's
absolute UV area, used as a lower area bound for each isolated selected face.

With zero budgets the original UV values are returned unchanged. Validation uses a copy
converted to the requested storage precision, fixed source triangulation and strict
within-chart overlap checks. Independent charts may overlap before packing. An invalid
layout returns `accepted=false`, `reason="invalid_layout"`; it is not silently repaired.
The report distinguishes skipped `raw_generation`, unchanged `postprocessing` and
`stored_validation`, including precision, pins, protected cuts, topology and overlap.

Charts come from explicit supplied cuts plus exact UV discontinuities along actual mesh
edges. Coincident positions or UVs do not weld disconnected geometry. Explicit cuts
remain cuts even when their UV endpoints coincide unless an explicit joining allowance
accepts a merge across an unlocked boundary between charts. Unrelated internal slits remain cuts. Locked edges retain their supplied
cut state; incompatible UV discontinuities make validation fail. Input-degenerate faces
are reported as exemptions, not silently lost.

A geometry exemption does not permit a self-crossing stored UV polygon. Such
crossings make the layout invalid and appear in `coverage_defect_faces`, even
when the supplied corners have not changed. The existing selected-region recovery
can then repair them before packing, subject to pins and protected cuts. Simple
UV coverage on tiny or otherwise geometry-exempt faces remains unchanged unless
selected for an operation.

Validation separates original input support from the supplied cut surface.
`stored_validation.input_topology_supported` is false when a used original vertex
has unsupported fans; `unsupported_input_vertices` gives its original vertex IDs.
`cut_sectors_valid` reports whether the supplied cut surface has no pinched sectors,
with current sector IDs in `pinched_sectors`. `topology_valid` retains its conservative
meaning: both original support and unpinched sectors. A valid cut surface does not
claim to repair the original mesh.

`uv_coverage_valid` reports finite coordinates, fixed-triangle coverage and the
existing exemption-aware within-chart tests independently of original fan support.
Final `valid`/`accepted` still require coverage, unpinched sectors, pins and protected
cuts. The result's `uv_valid` retains that final validated-layout meaning; it does not
ignore violated pin/cut constraints. Thus `accepted=true`, `uv_valid=true` and
`topology_valid=false` can describe usable UVs on an unsupported original fan.
Legacy pipeline policies and their output semantics are unchanged by this distinction.

Results own `chart_of_face`, `chart_of_corner`, `edge_vertices`, `seams`, `uv`,
`changed_faces` (empty when unchanged), exemptions and geometry fingerprint. Chart IDs
are local to the returned partition. All face and corner mappings use original input
order. Existing Session cancellation, stale-geometry checks, `result_for_apply` and
`close` still apply. Cancellation does not become a rejected optional improvement.

With `relax_faces_seconds > 0`, `selected_faces` selects original faces for twelve
synchronous neighbor-weighted relaxation rounds. The graph uses actual cut sectors,
not rounded UV coordinates. Sectors touching any unselected face or any pinned corner
remain fixed. Repeated original polygon edges keep their neighbor weights. Every round
reads the previous round; initially settled neighbors have weight one and unsettled
neighbors weight one quarter.

The proposed coordinates must pass stored-precision validity, fixed-triangle orientation
and strict within-chart overlap checks. An affected, previously valid chart cannot increase its P95
or maximum anisotropy beyond numerical tolerance (1e-6 relative); a repaired invalid
chart must have P95 anisotropy ratio at most two. This uses the native triangle-area
metric, rather than a face-weighted percentile. Unchanged charts (including exempt
degenerate charts without a distortion metric) do not participate in that comparison.
The original input is retained if any
check fails or the phase expires before adoption. Original cuts remain explicit and only
changed sectors are rounded in the returned buffer; unselected and pinned corner values
remain exactly as supplied. A rounded copy of the entire result is still validated.

`report.relax_faces` reports acceptance, reason, status, movable sectors and completed
rounds. `changed_faces` identifies only original faces with changed returned UVs. Empty
selection and zero allowance do no relaxation. The relaxation allowance includes graph
construction and proposal validation; cancellation propagates. Validation of the retained
incumbent still runs after a deadline. Relaxation skips nonfinite inputs; a separately
budgeted region phase may recover them. Local allowance expiry does not consume another
phase's allowance; both remain bounded by the global deadline.

## Artist-marked cuts in stored UVs

`repair.require_open_protected_cuts=true` additionally requires locked seam edges to
remain open after the requested storage conversion. The default is false, preserving
the existing logical-cut contract. This flag is supported only for `postprocess_layout`
with `repair_v1`; pipeline dispatch rejects it rather than silently ignoring it.
A marked edge means `seams[e] && locked[e]`. Locked continuous edges retain their
existing constraint and are not requested UV cuts.

Every pair of incident edge corners is compared in original vertex endpoint order,
including all incident pairs of a nonmanifold edge. A pair is closed when all four
coordinate differences are <= 1e-7. Boundary edges have
no pair. Validation measures the requested stored precision; a double gap that
collapses in float32 fails. Diagnostics contain `logical_protected_cuts_valid`,
`protected_uv_cuts_required`, `protected_uv_cuts_valid`, `closed_protected_uv_edges`
(original edge IDs) and `protected_cut_recovery_faces` (original face IDs).
`protected_cuts_valid` combines logical protection and the requested stored check.
When the flag is false, the stored check is not requested and its result is vacuously
true. These checks do not certify unsupported original topology as supported.

With a positive region allowance, closed pairs select their second incident face in
stable original corner order. Native per-face isolation runs before explicit selected
region work, sharing the same cumulative `region_recovery_seconds` deadline. It does
not spend a fresh allowance per face or reset the phase timer. The selected faces use
the existing projection, density and minimum-area rules; connected reconstruction is
skipped for this cut-opening proposal. Pins and locked continuous edges may refuse it.
A stored candidate must open all requested cuts, preserve pins and unpinched sectors,
leave every isolated face valid, and introduce no new strict coverage-defective face
IDs. Unrelated existing defects stay unchanged and keep result-level validity false.

`report.recover_protected_cuts` reports this proposal separately from explicit region
work. Zero/expired budgets leave UVs unchanged; cancellation propagates. All
relaxation/region/join/singleton candidates are subject to the same stored-cut check.
The addon must opt in when mapping artist marks; that integration is not yet switched.

## Connected regions and isolated faces

With `region_recovery_seconds > 0`, selected faces first receive individual convex
projections with the source density and optional minimum-area rule. Native adjacency
then identifies connected selected groups without crossing existing explicit or exact
UV cuts. Unknown nonfinite edge endpoints do not establish continuity. At most eight
groups of 2–256 faces, in decreasing face-count order, are considered; each must have
at most 1,024 corners and a valid disk of cut sectors.

Native LSCM must return all original corners in order, finite nonfolded coverage, no
strict within-group overlap, and face-area-weighted step-CDF P95 anisotropy ratio at
most two. Connected maps retain the total isolated area and are placed beyond the
input bounds. They are checked again after storage conversion. Ineligible/rejected
groups retain valid per-face projections. The phase shares one caller-supplied time
allowance across isolation, solves and final validation; an expired phase does not
commit a partial proposal. Zero allowance does no isolation, unlike the old low-level
Python helper's unconditional zero-time projection.

Synchronous C++ callers may pass `LayoutPostprocessServices` to the four-argument
`postprocess_layout` overload to supply a native region solver. The existing three-argument
entry point, Session, CLI and Python binding continue to use native LSCM. This service
also allows deterministic solver-originated cancellation tests. Supplied proposals cannot
bypass mapping, distortion or storage checks; cancellation, constraint and memory-limit
errors propagate rather than becoming failed optional improvements.

The reference orientation is the input's **global** signed fixed-triangle majority,
with positive tie-breaking. Nonfinite coordinates are replaced
by zero only in temporary orientation/density statistics. The original data remains
unchanged until a valid proposal is adopted. This preserves an entirely mirrored layout
without treating a minority flipped chart as a new orientation reference.

Unselected UVs and source geometry stay exact. A selected pinned corner refuses the
region proposal (`selected_pin`), because automatic placement/rescaling must not move
it. Locked interior continuity may survive a connected solve; a final partition that
needs a forbidden cut is rejected. The entire retained result still passes pin, cut,
topology and stored-coverage validation before acceptance.

`report.repair_regions` records acceptance, reason, solver attempts, connected groups,
remaining isolated faces and per-group outcomes. `stored_validation.coverage_checked_faces`
lists changed faces whose fixed UV triangles are checked even when their original 3D
geometry is reported as exempt. A self-crossing input polygon remains an input-geometry
exemption; a simple recovered UV layout does not claim to repair that geometry.

Native/binding results retain nonfinite coordinates when recovery is rejected or expires.
CLI result JSON represents each such coordinate as `null` with an exact `/uv/corner/axis`
entry in `unavailable` (`nan`, `positive_infinity`, or `negative_infinity`). The native JSON
reader restores those values and rejects missing, duplicate or inconsistent metadata.
Accepted or UV-valid results must contain only finite coordinates. JSON input has no
nonfinite numbers; use the existing binary UV array descriptor or the native/binding API.

## Chart joining modes

`join_seconds > 0` enables orientation-preserving similarity joining along actual
manifold mesh edges. `join_mode="small"` (default) allows moving charts of at most
four faces and caps the phase at one second. Explicit `join_mode="fewer"` removes the
face-count restriction and caps the phase at two seconds. Both are further limited
by the supplied `join_seconds` allowance and operation deadline. Charts are visited
in increasing face count and move toward a larger chart (earlier chart ID breaks
equal size ties). Up to eight shared-edge proposals, longest target edge first, are considered
per visited chart. One deadline spans the entire pass; it is not reset per proposal.

A proposal preserves target UV values exactly, copies matching endpoint sectors exactly,
and rounds moved coordinates to the requested precision. An independent rounded copy
also checks unchanged target coordinates. It cannot move a pin, add a discontinuity to
an existing connection, remove an unrelated slit or close an open protected cut via
another edge. The actual stored cut graph must have fewer charts. No moved face may
remain defective, and a merged chart may not add a defective original face relative
to the last accepted incumbent. Existing defects elsewhere, including distant defects
in the target chart, can remain unchanged and reported.
Joining also rejects newly introduced partial boundary contacts: an endpoint inside
another boundary segment, a partial collinear overlap, or a proper crossing. Different
endpoints along a shared line can round into an interior overlap during later packing.
The check uses stored coordinates and exact orientation predicates, not a proximity
tolerance. Exact matching endpoints and whole edges remain legal. A partial contact
already within one incumbent chart is retained as an existing slit; coincident
boundaries of separate unpacked charts do not authorize a new joined contact.
The `boundary_contact` rejection happens before adoption and does not alter UVs.
Enumeration obeys the phase stop callback; incomplete checks never authorize a join.
Whole-layout face-weighted P95 anisotropy ratio cannot increase by more than `1e-5`
from the phase's original layout. The reference is fixed across all proposals, so
individual accepted joins cannot repeatedly consume a new tolerance.
The tolerance is additive in ratio units; ratio-minus-one has the same increment.

`report.join_small_islands` records before/after chart counts, accepted joins, attempted
edge proposals, status and reason, plus `anisotropy_p95_before/after` in ratio units.
Rejection counts identify the failed guard; up to eight validation failures retain
original overlapping/changed-exempt face IDs and constraint results. Full aggregate
quality reports are constructed for retained results rather than every rejected
proposal. Candidate validation measures the entire newly merged chart and all of its
participating incumbent charts, including distant faces and unchanged input-degenerate
faces that still occupy UV space. A geometry exemption cannot hide newly defective
coverage. Pins, cut constraints and cut-sector topology still get global checks.
Every retained result receives both final whole-layout validation and a whole-layout
no-new-defect check before adoption; failure of the latter restores the phase incumbent.
`validation_scope` distinguishes partial diagnostic measurements from final reports.
`coverage_defects_before/after` records the strict phase comparison. Result-level
`accepted` and `uv_valid` remain false when unresolved coverage defects leave the
whole layout invalid, even if the joining phase reports an accepted improvement.
`remaining_invalid_charts` exposes that distinction.

`work_seconds` measures cumulative wall time in cut construction, quality measurement,
candidate validation and boundary-contact checks. `validated_face_visits` records candidate measurement work;
`baseline_face_visits` records the corresponding incumbent measurement work.
Up to 256 `proposal_records` retain original source/target corner IDs, moving-chart
size/first original face, and outcome. These bounded diagnostics do not renew the phase
budget. `build_info().join_modes` advertises the supported mode values.
Already validated joins can survive expiry of a later
optional attempt; no partially evaluated proposal is returned. Cancellation propagates.
A seam-only join is reported as repaired even when `changed_faces` is empty because that
array records coordinate changes. Independent disconnected geometry remains a separate
chart even when its UVs coincide; overlap between separate unpacked charts is permitted.

A valid partial improvement does not guarantee a requested fragmentation target.
Consumers must check both validity and the resulting chart counts under the declared
allowance. Budget expiration and precision changes need separate regression coverage.

The binding's `build_info()` and CLI `--version --json` share source/build identity,
schema versions, supported policies, operations and host capabilities. Check
`repair_operations` as well as the operation name: this candidate lists `relax_faces`
`repair_regions`, `join_small_islands`, `repair_singletons` and `recover_fragments`.
Source identity hashes native headers, sources,
binding/CLI sources, CMakeLists and the dependency lock at configure time, including local
edits.

## Singleton reconstruction

`singleton_seconds > 0` enables best-effort native reconstruction with a 45-second
phase cap, further limited by the supplied allowance and remaining operation deadline.
The same deadline covers all faces, target charts, solves and repeated passes. Real
mesh-edge cut sectors determine singletons; coincident UVs on disconnected geometry
never establish adjacency.

For each singleton, the engine first tries uniform harmonic attachment to adjacent
charts, expanding at most four rings and 128 faces while fixing sectors outside the
ring. A boundary-leaf attempt then matches one shared edge by an orientation-preserving
similarity and pins its free sectors, retaining the surrounding chart's constraints.
If attachment fails, it tries connected patches of 2–64 faces: complete small target
charts, neighbour pairs and up to three expanding rings. Native LSCM followed by a
convex-boundary fallback reconstructs disk patches. The patch preserves the previous
coverage area, follows the original global orientation and is placed beyond the
incumbent UV bounds. Source face/corner IDs and immutable geometry are retained. Numerical constraints,
edge alignment, source-area scaling and distortion measurements consistently use
the requested stored precision. Original return values remain separate: unchanged
representable coverage and explicit pins retain their exact supplied coordinates,
including doubles that are not exactly representable in float32.

Faces excluded from the 3D distortion metric do not become free solver variables
when their existing UV coverage is valid. Before singleton attachment, the engine
identifies source-degenerate faces whose stored fixed triangles have the expected
orientation and whose UV polygon is simple. Their sectors become exact harmonic
constraints, like the outside ring. Free patch reconstruction skips groups containing
these covered faces (`covered_exempt_patch` in rejection diagnostics). This prevents
an optional neighbour attachment from flattening healthy UVs that a percentile
metric cannot assess. The mask updates after accepted work, so newly recovered
coverage receives the same protection. Invalid or collapsed UVs remain eligible for
reconstruction; the rule does not claim to repair the source geometry.

Protected cuts exclude adjacency before solving. Locked continuity cannot be broken
by a patch boundary. Harmonic solves include user pins and exact outside-sector
constraints; conflicting constraints refuse the proposal. Automatic rescaling/placement
skips patches containing user pins. A candidate must reduce singleton count, connect
the requested singleton, preserve pins/cuts, and introduce no newly defective original face in any changed or repartitioned
chart. The singleton being attached must itself become valid. Defect measurement
includes fixed triangle folds/collapse, nonfinite values, native quality/boundary
failures and strict overlap participants, including geometry-exempt UV coverage.
The comparison uses the last accepted incumbent, so a previously repaired face
cannot regress merely because it was defective at phase start. Pins, cuts and
cut-sector topology remain mandatory global constraints. Existing unrelated defects
can remain in the same target chart and are still reported.

The whole-layout face-weighted step-CDF P95 ratio must be finite and at most
`max(previous_accepted_ratio, 1.25) + 1e-5`, equivalent to the source stretch envelope
`max(previous_stretch, .25) + 1e-5`. Unlike similarity joining, singleton reconstruction
updates that comparison after each accepted improvement, as the source does. A patch
also satisfies `max(previous_patch_ratio, 2.0) + 1e-5` before rescaling. Native fixed
triangulation can differ from the source's triangulation; metrics are not claimed
bitwise equivalent.

`report.repair_singletons` contains before/after counts, `attached`, `rebuilt_patches`,
solver attempts, rejection counts, stored ratio measurements and phase allowance.
`coverage_defects_before/after` lists defective original face IDs under the strict
candidate check. Whole-result `stored_validation.coverage_defect_faces` follows
the normal final exemption policy; phase measurements additionally include
unchanged geometry-exempt coverage. The masks are not claims of repaired 3D geometry.
`remaining_faces` lists all remaining singleton original face IDs; `unavoidable_faces`
only lists faces with no eligible mesh-edge neighbour (disconnected or protected).
A timeout does not prove that its unattempted faces are unavoidable. Cancellation
propagates; ordinary numerical or constraint rejection retains the incumbent, and
an accepted improvement survives later optional timeout. Zero allowance does no work.

The phase's `accepted` field describes a fragmentation improvement. The result-level
`accepted` and `uv_valid` still describe **whole-layout** validity after final validation.
An unrelated invalid component can therefore leave both result-level flags false even
when the phase improved another chart; `remaining_invalid_charts` makes that explicit.
Such a result must not be represented as a validated layout ready to apply.


## Native fragment proposals

A positive `host_proposal_seconds` enables an optional phase after selected-region
recovery and before singleton reconstruction/joining. `region_proposal(mesh,
original_corner_ids, stopped)` returns proposed UVs, and `fragment_unwrap(HostSolveRequest, stopped)`
returns a fixed-cut solve. Both return `optional<vector<Vec2>>`; `nullopt` means the
operation failed. The projection service must not be passed as a fixed-cut solver.
Callbacks are synchronous, use only their supplied scratch arrays, and must not retain
borrowed maps. No service is installed by default. The default unwrap method is
`ANGLE_BASED`; `fragment_unwrap_method` may instead
select `CONFORMAL` or `MINIMUM_STRETCH`. `SMART_PROJECT` is rejected as an unwrap method.

Native selection uses singleton **UV** charts, manifold face adjacency, inclusive
4–512-face connected groups, largest first with stable original-face ordering on ties,
and at most eight groups. This recovery does not use the separate warning detector's
45-degree threshold. Faces touching a protected cut are excluded; groups with pins are
rejected before callbacks. Explicit seam flags alone do not make continuous UVs eligible.
Native code creates compact region arrays in original face/corner order, derives cuts
from the projection and constructs an unwrap mesh with one vertex per cut sector.
Both callbacks retain the original-corner map. They never choose whether to adopt UVs.

After solving, native checks isolate exceptional damaged faces, require fewer charts,
strict local coverage and a face-area-weighted P95 anisotropy ratio at most 1.25
(stretch 0.25, where stretch is the ratio minus one). Independent regions may overlap before packing; inside a
recovered region, the inherited recovery rule additionally rejects cross-chart overlap.
Isolation uses a greedy cover of unique overlapping face pairs, selecting the highest
degree and breaking ties by original region-face order. It does not automatically isolate
both sides. Strict native overlap predicates remain in force. A truncated native pair
report (more than 4,096 triangle witnesses in either overlap pass) rejects the optional
proposal as `overlap_pair_budget`; an incomplete graph never proves coverage.
The proposal retains the selected original absolute triangle-based UV area and is parked
beyond the incumbent X extent. It is rounded to the requested precision, then coverage,
chart reduction, pins, cuts and the same stretch threshold are checked again. Fixed
native source triangles are used consistently. A uniformly mirrored host map is
normalized for local checks and returned with the incumbent's reference orientation.
Unselected coordinates are copied exactly, without rounding them in place.

`report.recover_fragments` records before/after chart counts, recovered regions/faces,
rejections, proposal calls, unwrap calls, and per-group reasons/quality. Malformed,
nonfinite or failed optional host answers reject that group. Host-failure exceptions
are treated likewise; cancellation and other native contract/resource errors propagate.
There is one cumulative phase budget. Expiry between projection and unwrap sends no
unwrap; earlier validated groups survive optional expiry. Direct callbacks receive the
cumulative `stopped` predicate and must honor it during waits. Native Session polls it
while awaiting a reply and retires the pending request on expiry. A running Blender
operator cannot be preempted; its late reply is stale and must not be applied. Rejected/unfinished groups preserve
their incumbent coordinates exactly. Final whole-layout validity is reported separately
from recovery success, so an unrelated pre-existing defect can still make the result
invalid.

See [SESSIONS.md](SESSIONS.md#optional-fragment-host-operations) for the distinct
projection/unwrap protocol and Python example. Fragment recovery success does not
hide unrelated pre-existing strict coverage defects in the final layout.


Singleton harmonic attachments require a complete constrained chart boundary.
Existing boundary edges on the receiving chart retain both endpoint UVs; free
leaf sectors use the leaf’s similarity mapping. An incomplete boundary is not
sent to the harmonic solver. With free boundary sectors, the Laplacian objective
can pull valid faces into nearly collinear configurations even while the global
95th-percentile distortion score improves. Boundary constraints prevent that
producer defect before packing or output validation. Unrelated existing defects
remain explicitly reported instead of triggering an incidental whole-patch
reconstruction merely because an unconstrained attachment failed.
