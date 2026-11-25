# Build and install

Requirements: CMake 3.28+, a C++20 compiler, a C compiler and Ninja (or a supported
CMake generator). Dependencies are pinned by `deps.lock.json` and its checksum in
`cmake/deps.lock.sha256`. CMake verifies downloaded archives; `CS_DOWNLOAD_CACHE`
can reuse an existing archive cache. Do not modify extracted dependency sources.

## Library and CLI

Build the shared library and CLI, run the smoke tests and install the SDK.

```sh
cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --parallel 4
ctest --test-dir build/release --output-on-failure
cmake --install build/release --prefix /path/to/cleanseams
cmake -S examples/installed -B build/consumer -DCMAKE_PREFIX_PATH=/path/to/cleanseams
cmake --build build/consumer
```

The SDK installs public headers and locked nlohmann/json headers, the shared core,
CLI, CMake targets, notices and dependency metadata. Eigen, SuperLU and BLAS remain
private implementation dependencies. The exported targets are `CleanSeams::core`
and `CleanSeams::cli`. Build consumers with compatible compiler/runtime versions.
Fast-math and floating-point contraction remain disabled.

## Python binding

Add `-DCS_BUILD_BINDINGS=ON -DPython_EXECUTABLE=/path/to/target/python`.
Headers must match that interpreter; set `Python_INCLUDE_DIR` when necessary.
The binding embeds an isolated, hidden copy of the native core so another installed
addon cannot interpose a different core library. It is colocated with the CLI/core
under the installed `bin` directory. The installed CMake package exposes
`CleanSeams_HAS_BINDING`, `CleanSeams_PYTHON_DIR` and `CleanSeams_PYTHON_SOABI`.

Blender packaging requires its exact Python ABI and matching header provenance.
Use `tools/probe_blender.py --help` to inspect the actual target Blender. On macOS
use `CMAKE_BUILD_WITH_INSTALL_RPATH=ON` and `CMAKE_INSTALL_RPATH=@loader_path`;
on Linux use a literal `$ORIGIN`. Preserve compiler/runtime notices when bundling
runtime libraries. See [packaging](PACKAGING.md).

For the qualified macOS Blender 5.2 target, set
`CMAKE_OSX_DEPLOYMENT_TARGET=15.0` (or use the `macos-release` preset).
Do not inherit a newer host OS minor version: its rounded-up wheel tag can be
rejected by Blender even when the native module loads locally.

On Windows use a native x64 compiler and the target Python headers/import library;
CMake also supports an explicit cross-toolchain. Build `--config Release` for
multi-configuration generators. Windows, Linux and macOS are separate build and
acceptance targets. A successful cross-build is not a native Blender test pass.

The `native-debug`, `native-release`, `asan` and `tsan` presets enable the
[smoke tests](TESTING.md). Set `CS_BUILD_TESTS=OFF` when only the distributable
targets are needed.
