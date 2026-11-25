# Boundaries, host lifecycle and failure behavior

Everything that enters the engine from outside is checked before native work uses it, and
anything that fails does so as a structured `ContractError` (Python `ValueError`, CLI exit
2/3/4/5), never as a crash, a hang, an unbounded allocation or a write to the user's data.

## Untrusted input

| Boundary | Checks |
|---|---|
| JSON text (case files, CLI `--options/--result`, binding option/report strings) | `parse_untrusted_json` (io.hpp): nesting at most 128 levels, counted before parsing; syntax, invalid UTF-8 and number overflow become `invalid_input`. |
| Binary array descriptors | Dtype, exact `product(shape)*itemsize` (overflow-checked), offset+length within the file, relative path without `.`/`..`, backslashes, drive letters or control characters (a NUL-truncated name opened another file), no symlink escape, unresolvable paths as `io_error`, SHA-256, finite floats (an incoming UV map may hold NaN or infinities: the run then has no baseline; a case file stores such a map as a binary array), the memory budget. |
| Host arrays (bindings) | Exact dtype, rank, width, C-contiguity and alignment (`bindings/arrays.*`); counts, index ranges, polygon loop ranges, edge sides and finiteness natively (`read_host_mesh`). Unpinned UVs may be non-finite; a pinned loop's UV must be finite. |
| Session input (`Session::start`) | Synchronously, before the worker starts: mesh (snapshot), the run's table memory (`memory_limit`), UV count (non-finite UVs mean no baseline), pin range, finiteness and uniqueness (a duplicate used to overwrite the first pin silently), anchors, mask length and 0/1 values, pairs and mask not both, options. |
| Host answers | Identity (session, request, geometry revision, protocol version), kind, a structured error exactly when unsuccessful, a known unwrap method, one UV per scratch corner (a non-finite map is ranked invalid and retried); typed Python fields. Stale, foreign, late and duplicate answers are refused and the request stays open. |
| CLI | Missing or unreadable files and malformed numeric options are invalid input (exit 2), not internal errors (exit 8). |

## Memory budget

A run has one allowance: `SessionConfig::memory_budget_bytes` for a session,
`PipelineInput::memory_budget` for a direct run (2 GiB by default); the Blender host and
the command line use `recommended_memory_budget()`, half the physical memory and at
least 2 GiB. The
run is refused up front with `memory_limit` (field `memory_budget`, with the bytes needed)
when its simultaneous tables do not fit: `run_table_bytes` = mesh + cut surface +
quality geometry with one full report + seam network, each charged per element from
measured allocations (`mesh_table_bytes`, `cut_table_bytes`, `quality_table_bytes`,
`seam_table_bytes` in `mesh.hpp`; the cut and seam charges assume every edge cut). Each table build checks its own charge. Worker tasks inherit
the run allowance through `active_memory_budget()`. Chart solves keep separate
bounds within that allowance.

## Degenerate geometry

Zero-area and coincident faces, non-manifold edge fans, bow-tie vertices, inconsistent
winding, non-orientable strips, duplicate faces, subnormal and huge coordinates, far
translations and seeded random degradations either fail validation (repeated corners;
lengths or areas beyond float64) or run through the whole pipeline without a deadline and
end by validated counts.

## Applying results in Blender

`ops.RunJob` starts a native session, answers its solve requests on the main thread and,
in `finish()`, writes only if the object, its mesh datablock, its UV layers, its seams/UVs/
pins and the native geometry fingerprint are still those the run started from
(`HostData.live_fingerprint`, `Session.result_for_apply`). All arrays are staged before the
first write; a pinned UV is compared as stored (float32 read exactly as float64), never after
aligning maps, and a result that moved one is refused. A write that fails midway is rolled
back from the snapshot. Deletion, undo, file loads, UV layer changes, seam/UV edits and
same-count coordinate edits during a run all end in "nothing applied".

Written UVs are exactly the float32 rounding of the validated float64 map (at most half a
float32 ulp away); the engine's own validity is measured in float64, separately.

Shared mesh data is updated for every user, as in the source add-on, and the status names
the other objects. Library-linked or non-editable mesh data is refused before the run.

## Shutdown

Every binding `Session` is registered. `close()` cancels, wakes and joins the worker;
concurrent closes serialize. Unregistering the extension and loading a file cancel all runs;
`native_module.deactivate()` calls `close_all_sessions()`; the module's interpreter-exit hook
closes whatever is still open when Blender quits. No native worker outlives the module.
