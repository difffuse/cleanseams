# Native contracts

The C++ library owns defaults, validation, UI conversions and versioned message
formats. `include/cleanseams/options.hpp` defines defaults;
`include/cleanseams/detail/option_fields.hpp` enumerates fields for serialization and
schema generation. `option_registry()` exposes UI metadata. The checked-in `schema/`
files describe the wire formats; `tools/export_schemas.py` exports them from a built
CLI for comparison.

Partial option patches fill omitted fields only. Explicit zeros remain explicit.
The core profile defaults to `custom`; the Blender UI defaults to `hard_surface`.
A core deadline of zero expires immediately; UI zero resolves to no deadline.
UI percentages convert to fractions once. Host context discovery remains in Blender;
option resolution and object-space direction conversion happen in the native core.
JSON null represents an unbounded `layout.max_defect`. NaN/infinity are rejected.

Policies are `native_baseline` (default), `legacy_010` (explicit compatibility mode),
`enhanced_v1` (optional seam planning) and `repair_v1` for supplied-layout processing.
Reports name the policy and operation. Compatibility settings remain supported for
existing callers; they are not interchangeable with current defaults.

Inputs must satisfy the schema and native range/finite checks. Negative deadlines,
invalid target-island bands, negative routing weights, empty variants, invalid native
enums and malformed arrays are rejected. Callers must use the returned report and
accepted/valid flags; successful transport alone does not establish valid output.

[Sessions](SESSIONS.md) describe ownership and host requests.
[Supplied layouts](LAYOUT_POSTPROCESS.md) describe repair operations, storage precision,
limits and original-corner mappings. [Mesh analysis](MESH_ANALYSIS.md) supplies fixed
triangulation and source mappings. [Enhanced seams](ENHANCED_SEAMS.md) describes the
optional policy. Native build metadata identifies supported operations and capabilities;
check it before selecting an operation rather than inferring support from a filename.
