# Mesh identity and host transfer

`MeshSnapshot` owns validated float64 positions, original polygon offsets and
corners, host edge rows, and host identities. It exposes two complementary
operations: `geometry_key()` is a 32-byte BLAKE2b content digest rendered as
64 lowercase hex characters; `same_geometry()` compares every stored geometry
array exactly, including the bits of signed zero. A cache must compare exact
geometry after a key hit. The supplied `MeshInput.geometry_fingerprint` string
does not establish native identity and is excluded from both operations.

The version-1 key hashes the bytes of `cleanseams.mesh.geometry.v1` followed by
one zero byte, then the following arrays in order. Each array begins with its
uint64 little-endian row count:

| Array | Encoding per row |
|---|---|
| positions | three float64 little-endian values |
| face offsets | uint64 little-endian |
| corner vertices | int64 little-endian |
| host edge endpoints | two int64 little-endian values |
| explicit host edge identities | int64 little-endian |
| explicit host loop identities | int64 little-endian |

Absent identity arrays and explicitly supplied arrays have different content
keys even when their mappings are equivalent. This conservative choice cannot
reuse an incompatible snapshot. The content key is distinct from the monotonic
`Revision` ordinal assigned by a persistent Session.

`HostMeshView` borrows six raw typed arrays only for a synchronous call. No view
is stored. `read_host_mesh()` gathers corners in polygon order into owned
`MeshInput`, retaining the corner-to-host-loop mapping and every host edge row,
including wires. `MeshSnapshot` then computes canonical core edge identities and
both mapping directions. The binding transfers complete contiguous arrays and
returns owned arrays; it rejects implicit dtype conversions, strided storage,
bad alignment and malformed shapes.

`read_host_attributes()` moves edge seam/lock flags and loop UV/pin arrays into
the snapshot's core ordering, extracts point anchors, and emits corner pin
constraints without merging corners that share a vertex. Wire edge flags do not
become surface seams. Missing masks become false; missing UVs remain absent.
Finite coordinates, mask values, lengths and loop identities are checked before
gathering. A pinned loop's UV must be finite; other UVs may be NaN or infinite, and the
run then measures no baseline. The binding accepts an optional `attributes` dictionary with boolean
NumPy masks and float64 UVs. The Blender integration discovers host attributes and handles layer and context operations.

`host_geometry_fingerprint()` reproduces the frozen adapter's restore token:
BLAKE2b with a 12-byte digest over vertex/edge/polygon/loop int64 counts, float32
position bytes, and int32 edge endpoints, loop vertex IDs, loop edge IDs,
polygon loop starts and totals. Order and signed-zero bits are preserved.
Encoding is explicitly little-endian, matching the recorded Blender host and
the initial supported architectures. This does not assert legacy fingerprint
compatibility with an untested big-endian Blender host.

The restore token always describes raw object-space Blender data, even when
the requested native positions are transformed into world space. Default
positions remain in object space. World-space positions follow the source's
upper 3-by-3 linear transform plus translation; there is no projective divide.
Scale diagnostics use the supplied `obj.scale`, not a decomposition of
`matrix_world`, preserving parenting/shear behavior. Input and transformed
coordinates must be finite; overflow is an explicit error.

The streaming hash implements the unkeyed sequential algorithm described in
[RFC 7693](https://www.rfc-editor.org/rfc/rfc7693.html).

The coarse [mesh-analysis API](MESH_ANALYSIS.md) returns these original/host mappings
alongside fixed source triangle corner IDs and triangle-based surface areas, without
creating a Session or invoking an unwrap backend.
