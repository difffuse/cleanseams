# Fixed source mesh analysis (development candidate)

`MeshAnalysis analyze_mesh(const MeshInput&, uint64_t memory_budget)` exports the
native `MeshSnapshot` triangulation and original-input mappings without running seam
generation, an unwrap solver, postprocessing or packing. Its compiled Python entry
point is `cleanseams_native.analyze_mesh(mesh_payload, memory_budget=2147483648)`.
`build_info().features.mesh_analysis` advertises this helper. It does not import the
Python engine, SciPy or Blender modules.

The Python payload requires contiguous, aligned NumPy `positions[V,3]` (float64),
`face_offsets[F+1]` (uint64) and `corner_vertices[C]` (int64). Optional
`edge_vertices[H,2]`, `host_edge_indices[H]` and `host_loop_indices[C]` use int64.
These are geometry and host-identity arrays; unknown fields, UVs and solver settings
are rejected. Host edge/loop identities obey the existing [mesh identity contract](MESH_IDENTITY.md).
An absent or empty identity array uses the existing default mapping. No geometry is welded,
reordered or deleted. The binding copies inputs while holding the GIL, releases it for
native analysis, then returns owned, writable, contiguous arrays.

| Result | Type/shape | Meaning |
|---|---|---|
| `triangle_corners` | int64 `[T,3]` | Original polygon corner IDs in native fixed triangulation order |
| `triangle_faces` | int64 `[T]` | Original face ID of each triangle |
| `triangle_edges` | int64 `[T,3]` | Original core edge ID of each triangle side; -1 for an artificial polygon diagonal |
| `face_area_3d` | float64 `[F]` | Sum of unsigned 3D areas of the fixed triangles of each original face |
| `edge_vertices` | int64 `[E,2]` | Canonical sorted surface edges, retaining original vertex IDs |
| `corner_faces`, `corner_edges` | int64 `[C]` | Original face and canonical outgoing edge for every original corner |
| `host_edge_to_core` | int64 `[H]` | Supplied host edge row to canonical edge; -1 for a wire |
| `core_edge_to_host` | int64 `[E]` | Canonical edge to explicit host identity (or supplied row); -1 when absent |
| `corner_to_host_loop` | int64 `[C]` | Original corner to supplied host loop identity (or original corner index) |
| `face_degenerate_kind` | uint8 `[F]` | Existing native input-degeneracy classification |
| `triangulation_fallbacks` | uint64 scalar | Number of polygons for which native triangulation fell back to a fan |
| `geometry_key` | string | Native geometry-content key; distinct from a host restore fingerprint |

Zero-face inputs return correctly shaped empty arrays. Returned buffers survive input
mutation/deletion and do not share storage with later calls. Original vertex IDs can be
obtained by indexing the caller's `corner_vertices` with `triangle_corners`; the helper
never invents a replacement vertex/corner numbering.

`face_area_3d` deliberately uses triangle surface area, not the Newell/projected
polygon area in `MeshGeometry.face_area`. Those quantities differ on nonplanar faces.
Source UVs cannot select a different diagonal or hide a fold. Quad diagonal choice and
n-gon ear clipping are exactly the existing native policy. This is not a promise of
triangle-list equality with other triangulators, including convex n-gon fans.
Consumers must use one canonical triangle set consistently for area, distortion
and fold checks.

`face_degenerate_kind` uses the existing native codes: 0 ordinary, 1 zero area,
2 twisted/self-intersecting, 3 collapsed and 4 sliver. These are input-geometry
classifications under the native relative thresholds, not UV validity or topology
support guarantees.

Fallback triangles remain visible, with their fallback count and original degeneracy
classification. Exporting them does not certify that an unsupported/self-intersecting
polygon has a valid interior, and this API does not validate a UV map. Within-chart UV
validation and post-pack global overlap checks retain their separate contracts.

`mesh_analysis_bytes(TableCounts)` provides a conservative allocation allowance,
saturating at uint64 max on overflow. It includes native snapshot storage, one ingestion
copy and two owned analysis/transfer copies. The binding checks counts, identity-array
lengths and this allowance before copying input buffers. Native analysis also accounts
for any supplied `MeshInput.geometry_fingerprint` string when checking its budget.
Malformed data or insufficient allowance raises the existing structured native errors
(or a binding `ValueError` for buffer/type errors). This is an allowance for the
operation's geometry allocations, not a process-wide RSS limit or a charge for
caller-owned Python input arrays.
