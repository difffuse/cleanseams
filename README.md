# CleanSeams

CleanSeams provides seam generation, UV unwrapping, supplied-layout validation and
bounded recovery through a C++20 library, CLI and CPython binding. The optional
standalone Blender integration contains host glue only.

- [Build and install](docs/BUILD.md), [testing](docs/TESTING.md), [packaging](docs/PACKAGING.md).
- [Options and message contracts](docs/CONTRACTS.md), [sessions](docs/SESSIONS.md),
  [mesh analysis](docs/MESH_ANALYSIS.md), [supplied layouts](docs/LAYOUT_POSTPROCESS.md),
  [region recovery](docs/REGION_RECOVERY.md), [CLI and packing](docs/IO_PACKING.md).
- [Blender integration](docs/BLENDER.md) and [enhanced seam policy](docs/ENHANCED_SEAMS.md).

Use `find_package(CleanSeams CONFIG REQUIRED)` and link `CleanSeams::core`.
The [installed consumer](examples/installed/CMakeLists.txt) demonstrates the SDK.
The CPython module must match the target interpreter's ABI; Blender 5.2 uses CPython
3.13. There is no stable cross-version C++ or Python ABI promise.

The native version is `0.2.0`; the optional standalone Blender extension is `0.3.0`.
Identify exact builds by `cleanseams_cli --version --json`, source revision, native
source digest and archive hashes. Cross compilation does not qualify a native host.

Code is [GPL-3.0-or-later](COPYING). Preserve [dependency notices](LICENSES/README.md).
Run the CLI and C++ API smoke tests with `ctest --preset native-debug`.

## Upgrading from native 0.1 / extension 0.2

External packing integration has been removed. The CLI no longer offers `pack`,
`pack-request` or `pack-apply`, and the binding no longer exposes an external packer
client. Native seam generation and UV processing remain the library's responsibility;
packing belongs to the caller. The standalone Blender extension offers only Blender
packing. An obsolete saved packer selection cancels without modifying UVs; select
Blender explicitly to pack. The external executable setting is no longer accepted.

Generic process utilities are declared in `cleanseams/process.hpp`. Recompile C++
consumers against the matching SDK; the project does not promise a stable C++ ABI.
