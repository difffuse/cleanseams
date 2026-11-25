# Native sessions, host requests and bindings

`cleanseams::Session` (`include/cleanseams/session.hpp`, `src/engine/session.cpp`) runs
one native pipeline on one worker thread:

* `start(RunInput)` copies every buffer into an owned `MeshSnapshot`/`PipelineInput`
  synchronously; caller buffers may change or disappear immediately afterwards.
* `poll_progress()` drains stage events; each is also a cancellation checkpoint.
* With `SessionConfig::backend = "host"`, every solve becomes a typed `HostRequest`
  (`kind = "solve"`): the source BlenderBackend's isolated chart mesh (one vertex per
  cut sector), scratch-to-original corner map, planar n-gon start map, pins on scratch
  corners, method and iterations. Which charts are solved, planar n-gon maps and the
  conformal retry after a poor MINIMUM_STRETCH map remain native.
* `SessionConfig::host_unwrap_method` (Python `Session(host_unwrap_method=...)`)
  selects the primary generation solve: ANGLE_BASED, CONFORMAL or MINIMUM_STRETCH.
  It defaults to MINIMUM_STRETCH, preserving existing clients. It is independent
  of `fragment_unwrap_method`, whose default remains ANGLE_BASED. The host executes
  each request's method and reports the method actually used; it must not replace
  a native CONFORMAL retry with the primary selection. SMART_PROJECT is a region
  proposal operation, not an accepted primary unwrap method. Build metadata
  advertises this additive configuration as `features.host_method_selection`.
* `next_host_request()` hands a request out once; `submit_host_response()` accepts only
  the matching session/request/geometry/state revisions, chart keys, kind and protocol version with one
  UV per scratch corner — stale, duplicate, foreign, other-version and malformed
  answers are refused and the request stays open. A non-finite map is an answer, not a
  protocol error: as the source's `adapter.BlenderBackend`, it ranks invalid
  (Python tuple order, where a NaN decides), a minimum-stretch solve is retried with
  CONFORMAL, and the repair ladder follows.
* `request_cancel()` and `close()` wake a worker waiting for the host; `close()` joins
  it. Results carry the compatibility decision plus independent flags (uv/topology
  validity, budget, partial, cancelled, timed out), seams, map, exempt faces and the
  finite report.
* `start()` validates the UV count (non-finite UVs are kept: the run then has no
  baseline, as the source's `np.isfinite(uv_in).all()`), pin range/finiteness/uniqueness,
  anchors, edge-mask lengths and 0/1 values, pairs versus mask, and options before any
  worker exists, and the run's tables against `SessionConfig::memory_budget_bytes`
  (`check_run_memory`, a `memory_limit` refusal), which the run then uses for
  every table, solve and exposure field (`PipelineInput::memory_budget`). The Blender
  host passes `recommended_memory_budget()` (half the physical memory, at least 2 GiB).
  Answers also need a structured error exactly when unsuccessful and a
  known unwrap method. `result_for_apply(fingerprint)` returns the result only when the
  host's fingerprint of its live data equals the one the run started with (stale_state
  otherwise). Concurrent `close()` calls serialize; one joins.

The CPython module registers every `Session`: `close_all_sessions()` and
`live_sessions()` serve extension shutdown, and an interpreter-exit hook closes anything
still open (see [SAFETY.md](SAFETY.md)).

The worker never holds a Python object or calls into the host: every exchange goes
through the session queues. The CPython module (`bindings/module.cpp`) exposes
`Session` with exact dtype/rank/contiguity/finite checks, owned copies and owned NumPy
results; `ContractError` maps to a `ValueError` subclass. The GIL is released only
around pure native waits.

## Supplied layouts

The development `postprocess_layout` operation and owned chart result arrays are described
in [LAYOUT_POSTPROCESS.md](LAYOUT_POSTPROCESS.md). Check `repair_operations` and
`host_capabilities` before selecting optional algorithms; integration qualification remains
separate from capability availability.


## Optional fragment host operations

A positive `repair.host_proposal_seconds` enables fragment recovery only for
`postprocess_layout` / `repair_v1`. Session adds an explicit capability opt-in:

```python
session = cleanseams_native.Session(
    backend="superlu",  # global solver selection; fragment calls below use the host
    host_capabilities=["region_proposal"],
    fragment_unwrap_method="ANGLE_BASED",
)
```

The list defaults to empty, including for `backend="host"`; old host solve clients never
receive a new operation implicitly. Unknown/duplicate capabilities are rejected.
`fragment_unwrap_method` accepts ANGLE_BASED, CONFORMAL or MINIMUM_STRETCH, and defaults
to ANGLE_BASED. SMART_PROJECT is a projection, not an unwrap method.

The native worker selects a region and issues `kind="region_proposal"`. Its wire payload
is `HostRegionProposalRequest { mesh, scratch_to_original_corner }` under
`HostRequest.region_proposal`; `solve` and `pack` are null. The Python binding follows the
existing flattened solve-array format: `positions`, `face_offsets`, `corner_vertices`,
`scratch_to_original_corner`, plus `method="SMART_PROJECT"`, `kind` and `identity`.
Every array owns its memory. The host projects that scratch mesh and replies:

```python
reply = dict(identity=request["identity"], kind=request["kind"],
             success=True, method=request["method"], uv=proposed_uv)
session.submit_host_response(reply)
```

The host must execute Blender operations on its main thread. It does not choose regions,
derive cuts or accept replacements. Native code derives cuts and sends a separate
ordinary `kind="solve"` request whose scratch vertices are cut sectors. Its corner map
still identifies the original mesh corners. That request carries initial UVs, pins,
method and iterations in the existing solve format. The response echoes its identity,
kind and actual method. Returning `success=False` requires an error (a string in the
binding, structured `ErrorInfo` on the wire).

Each projection/unwrap pair shares a positive `state_revision`; the next projection
increments it. Every request also has its own monotonically increasing request ID and
the original geometry revision. Echo the entire identity, including chart keys. Wrong
identity, kind, method or shape does not consume the pending request. Nonfinite UVs are
consumed and cause optional group rejection, preserving the incumbent. Legacy solves
retain their existing nonfinite retry policy. Failed projection never triggers unwrap.

The queue hands out a request once. Cancellation takes priority over queued answers and
prevents another dispatch. Close cancels and joins; returned scratch arrays remain owned
by the caller. While an optional request waits, the worker checks its cumulative phase
and operation deadlines at intervals of at most 10 ms, subject to scheduling. Expiry
retires the request; a late response raises `stale_state`. The host should collect the
finished result after an operator returns late, rather than apply the orphaned proposal.
Earlier validated regions may survive timeout, while unfinished groups remain unchanged.

`recover_fragments.host_service_available`, `solver_backend`, `unwrap_methods`, proposal
and unwrap attempt counts make mixed native/host runs explicit. The top-level
`unwrap_calls` includes fragment unwraps and `host_proposal_calls` counts projections;
a global `backend="superlu"` label does not describe these optional host solves.
`before`/`after` count logical cut charts, consistent with returned chart membership;
explicit coincident cuts can count more charts than a coordinate-tolerance UV-island metric.

The host-request wire schema adds an optional `region_proposal` field; old solve/pack
documents without it still decode. New Session configuration fields also have defaults
and are optional in the schema. Protocol version remains 1 because capability opt-in
prevents dispatching this addition to unaware clients. Build metadata advertises both
`host_capabilities: region_proposal` and `repair_operations: recover_fragments`.
