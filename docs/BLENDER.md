# Blender extension

`blender/cleanseams/` is the extension (id `cleanseams`, Blender 5.2+). It keeps the
reference add-on's settings (every `CleanSeamsSettings` property), panels (3D viewport and
UV editor sidebars) and operators (`cleanseams.run`, `cleanseams.restore`,
`cleanseams.pack`, `cleanseams.report`). Everything that computes lives in the native
engine `cleanseams_native`:

| Concern | Where |
|---|---|
| Geometry mapping, host edge/loop maps, restore fingerprint | `read_host_mesh`, `host_geometry_fingerprint` (native) |
| Profiles, budgets, percentages, view direction | `apply_ui_profile_json`, `resolve_ui_options_json` (native) |
| Packing settings (tiles as Python `int()`) | `resolve_ui_pack_options_json` (native; packing only) |
| Pins kept exactly, float32 storage of the written UVs | `Session.result_for_apply`, `check_stored_uv` (native) |
| Generation, refinement, unwrap policy, repair, contours, hygiene | `Session` worker (native pipeline) |
| Which charts Blender solves, planar n-gon maps, conformal retry | native host backend (session) |
| Status line, summary, refusal and repairs text, finite report JSON | `run_messages` (native) |
| Enhanced controls: kind, label, range, items, default, panel sections; their resolution | `option_registry_json` (`ui_enhanced_fields`, `ui_sections`), `resolve_ui_options_json` (native) |
| The panel's run and packing lines (why proposals were refused, budgets, symmetry, packer verification) | `run_messages` (`messages`), `packing_messages` (native) |
| RNA properties, drawing, bulk `foreach_get/foreach_set`, Blender's unwrap | Python glue (`props.py`, `ui.py`, `host.py`, `ops.py`) |

## Run

`cleanseams.run` reads the object's raw float32/int32 buffers on the main thread, maps them
natively, resolves the run options (policy `native_baseline`, the product default, named
in the run report), starts a `Session(backend="host")` and services its solve requests on the main
thread: each request is an isolated chart mesh (one vertex per cut sector) with a start map
and pins; Blender unwraps it with the RNA-discovered method (`MINIMUM_STRETCH` when
offered, else `CONFORMAL`), `fill_holes=False`, `correct_aspect=False` and the configured
iterations, and answers with the method it used. The native side validates every answer
(identity, revision, size; a non-finite map is ranked invalid and retried with CONFORMAL). The result is written back only when accepted and
`Session.result_for_apply` finds every pin as given: a pin the host solver returned within the
pipeline's own pin tolerance (numpy.allclose, as the 0.1.0 add-on accepted and applied it) is
written back exactly on its stored value, unless that would flip or collapse a face (
Blender's CONFORMAL unwrap returns pins a few float32 steps off); a pin beyond the tolerance is
refused and nothing is written. After writing, the glue reads
the layer back and `check_stored_uv` decides the float32 round trip exactly (Blender stores
float32 UVs, so the stored map must be the float32 rounding of the validated map, nothing
else). Both decisions are native; the glue only gathers the arrays. Nothing is written
before the result is applied; a write that fails midway is rolled back from the snapshot.
Packing resolves only the packing settings (`resolve_ui_pack_options_json`), so run
settings such as the current view or the island target never block it, and a malformed
tile list never blocks a run.

The run is an `ops.RunJob` (start, step, finish). `finish()` applies only when the
object, its mesh datablock, its UV layers, its seams/UVs/pins and the native geometry
fingerprint are unchanged (`Session.result_for_apply`); deletion, undo, file loads, UV
layer changes and same-count coordinate edits during a run end in "nothing applied".
Pinned UVs are compared as stored. Shared mesh data is updated for every user and the
status names them; library-linked or non-editable mesh data is refused before the run.
Unregistering and file loads cancel and join runs; sessions left open at exit are closed
by the module. Details: [SAFETY.md](SAFETY.md).

Before reading the object, the run asks the native options which
attribute names to transfer (`painted_attribute_names`: `cleanseams_required`,
`cleanseams_forbidden`, `cleanseams_allow`, `cleanseams_avoid`, or configured aliases). The list is empty unless the run's policy is `enhanced_v1`: a
`native_baseline` run transfers no painted attribute and excludes those attributes from its stale-state digest.
`host.HostData` reads each named attribute exactly as Blender stores it (domain, data type
and, for a scalar type, its values in host order), `read_host_mesh` maps it natively to core
edges, vertices, faces or corners, and the session receives it as
`constraints.attributes`. Domains, types, ranges, clamping, conflicts and their effect are
native (`resolve_preferences`, [ENHANCED_SEAMS.md](ENHANCED_SEAMS.md)); only `enhanced_v1`
reads them; the panel selects it with its seam policy control and names the four attributes.
The painted attributes are part of the run's state digest: painting during a run makes the
result stale. `cleanseams_lock` and `cleanseams_anchor` keep their existing reading.

## Enhanced controls and messages

The panel's **Seam policy** (after the profile) chooses *Native baseline* (the default, also for
older projects without a recorded seam policy) or *Enhanced* (`enhanced_v1`). *Enhanced settings* is a
collapsible box in the native sections (features with their switches and confidence
percentages, symmetry, distances, painted attribute names, search bounds), greyed out unless
the policy is *Enhanced*; [ENHANCED_SEAMS.md](ENHANCED_SEAMS.md) lists every control and its
option. `props.py` declares these properties at registration from `option_registry_json()`
(each `ui_enhanced_fields` entry's metadata and UI schema default) and `settings_json`
transfers them with the others; `ui.py` draws the registry's `ui_sections`. Without the
engine they are not declared, a file's saved values stay in its ID properties, and the
operators stay disabled. A project saved by the frozen 0.1.0 add-on retains supported
settings, with *Native baseline* and the native defaults for the new controls. An obsolete packing selection is refused
until Blender is explicitly selected.

After a run, the status line is the native summary (with an `enhanced: …` clause) or refusal,
and the box below it draws `run_messages`' lines exactly, one icon per level (info, warning,
error): the policy and outcome, why proposals were not taken, budgets, required relief slits,
painted attributes (including ignored ones and conflicts with their edge IDs), the solver, the
symmetry mode (strict, approximate report-only, asymmetric) and the distance policy, stops and
time budgets. Packing uses Blender's packer and appends a native line identifying its result
as not verified by the engine. Restore clears the lines. Restore after an enhanced run is the
same snapshot restore (below): the pre-run seams, UVs and pins on the identical geometry.

## Restore

A snapshot of seams, UVs and pins is kept per object with the native geometry fingerprint
(the source adapter's exact byte fingerprint of positions and index order). Restore writes
it back only onto identical geometry; any coordinate or topology edit refuses it.

Kept from the source add-on (measured in Blender, not changed here):

- Renaming the snapshot's UV layer is not a geometry change: restore puts the seams back,
  skips the UVs and pins of the missing layer and still reports "restored the previous
  seams and UVs".
- Snapshots are keyed by object name and live only while the extension is registered;
  unregistering clears them. A deleted object's snapshot stays in memory until then, and
  the fingerprint check guards any later object with the same name.
- Undo is independent: undoing a run returns the exact pre-run state, and the snapshot
  can still be restored onto it (Blender's background mode needs an explicit
  `ed.undo_push()` to start its undo system).

## Native module

`native_module.py` loads `cleanseams_native` from `CLEANSEAMS_NATIVE_DIR` or the package's
`native/` directory. A missing module or one built for another Python ABI is reported in
the panel and disables the operators.

See [Build and install](BUILD.md#python-binding) for the native module build.
