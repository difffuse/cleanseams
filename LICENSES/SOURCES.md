# Corresponding source and build recipe

The binaries in the CleanSeams artifacts (`libcleanseams_core.dylib`, `cleanseams_cli`,
`cleanseams_native.<SOABI>.so`) are built from:

1. **The CleanSeams project source**: `cleanseams-<version>-source.tar.gz`, written next to
   the binaries by `tools/package_extension.py --source-package`. It holds every tracked
   project file except generated build, distribution and execution outputs
   (`build/`, `dist/`, `execution/`). The CLI and C++ smoke checks are included.
   `SOURCE_PACKAGE.json` inside records the
   producer commit and whether its working tree was clean. When packaging from that Git
   checkout, `artifacts.json` records the same commit. An extracted archive has no Git
   history: the rebuilt artifact manifest reports `source_state: "unavailable"` with null
   revision/cleanliness fields. Retain the original receipt and archive hash separately;
   the receipt alone does not verify the contents of an extracted or modified tree.
2. **The third-party sources**: `cleanseams-<version>-third-party-sources.tar` holds the exact
   upstream archives the build verified, named `<name>-<sha256>.tar.gz`, and
   `deps.lock.json`:

   | Name | Version | Upstream commit | Archive SHA-256 | URL |
   |---|---|---|---|---|
   | eigen | 3.4.0 | 3147391d946bb4b6c68edd901f2add6ac1f31f8c | 8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72 | https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz |
   | superlu | 7.0.0 | 9c2665629b1fcc003726c87118dd1ccb895d94f8 | d7b91d4e0bb52644ca74c1a4dd466a694ddf1244a7bbf93cb453e8ca1f6527eb | https://codeload.github.com/xiaoyeli/superlu/tar.gz/refs/tags/v7.0.0 |
   | pybind11 | 3.0.1 | f5fbe867d2d26e4a0a9177a51f6e568868ad3dc8 | 741633da746b7c738bb71f1854f957b9da660bcd2dce68d71949037f0969d0ca | https://codeload.github.com/pybind/pybind11/tar.gz/refs/tags/v3.0.1 |
   | json | 3.12.0 | 55f93686c01528224f448c19128836e7df245f72 | 4b92eb0c06d10683f7447ce9406cb97cd4b453be18d7279320f7b2f025c10187 | https://codeload.github.com/nlohmann/json/tar.gz/refs/tags/v3.12.0 |

   The build extracts these archives unchanged; SuperLU's allocator hooks are injected with
   a forced include (`src/solver/superlu_hooks.h`), not by editing its files.
3. **The CPython 3.13.13 header set** the module is compiled against: the official
   `Python-3.13.13.tgz` (sha256 `f9cde7b0e2ec8165d7326e2a0f59ea2686ce9d0c617dbbb3d66a7e54d31b74b9`)
   `Include/` plus Blender's own `pyconfig.h`; see `host-header-provenance.json`. It is a
   build SDK, not shipped.

## Build recipe (macOS arm64)

The current component recipe targets macOS 15.0; runtime checks were performed on macOS
15.7.9 with Blender 5.2.1 LTS. This does not qualify the minimum OS or another platform.
Use CMake, Ninja, a C++20 compiler and Python 3.11+ for the tools. Recreate the target
Python headers using [BUILD.md](../docs/BUILD.md#python-binding) before configuring.

For an extracted source package, also extract the accompanying third-party archive
beside it. From inside `cleanseams-<version>-source/`, the cache below then resolves to
`../third-party-sources/`. CMake verifies the archive hashes against `deps.lock.json`.
For a development checkout, use its populated `build/downloads` cache instead.

```sh
# 1. Facts of the target Blender (Python version, SOABI, architecture, NumPy).
python3 tools/probe_blender.py --blender /Applications/Blender.app/Contents/MacOS/Blender \
  --output build/host.json
# 2. The release bindings build for that Python, relocatable (@loader_path).
cmake -S . -B build/blender-python-release -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCS_BUILD_BINDINGS=ON -DCS_BUILD_TESTS=OFF \
  -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=15.0 \
  -DCS_DOWNLOAD_CACHE="$PWD/../third-party-sources" \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_INSTALL_RPATH=@loader_path \
  -DPython_EXECUTABLE=/Applications/Blender.app/Contents/Resources/5.2/python/bin/python3.13 \
  -DPython_INCLUDE_DIR="$PWD/build/host-headers/python3.13"
cmake --build build/blender-python-release -j4
# 3. Binary artifacts (extension, wheel, CLI archive, artifacts.json).
python3 tools/package_extension.py --host build/host.json \
  --build build/blender-python-release --output dist

```

Maintainers generating new source/dependency packages from a Git checkout add
`--source-package` to step 3. That option requires tracked project files and verified
archives under `build/downloads`; it cannot reproduce source packages from an extracted
archive.

Compiler policy for every build: C++20, `-fno-fast-math -ffp-contract=off`, no
architecture-specific flags (`cmake/CompilerPolicy.cmake`). The packaging procedure
and supported artifact targets are described in `docs/PACKAGING.md`.

## Build recipe (Linux x86-64)

This candidate recipe targets Blender 5.2.1, CPython 3.13.13 and the header inventory in
`linux-host-header-provenance.json`. SDK and wheel checks ran under x86 emulation on Apple
arm64; this does not qualify native Linux hardware or the complete release corpus.
Use CMake 4.1+, a C++20 compiler, GNU Make (or Ninja), Git, GNU `readelf`, `patchelf`, and
Python 3.11+ with `tools/requirements-linux.txt` installed using `--require-hashes`.
The wheel is audited against `manylinux_2_35_x86_64`; a newer compiler/runtime requirement
that cannot satisfy that policy must fail rather than receive an inaccurate tag.

The Linux Blender runtime used here omits its headers. Recreate the official CPython
`Include/` directory as in the Mac recipe, but use the pinned **Linux** Blender SDK
configuration below. Fetching only Git metadata and the named blob avoids checking out the
large Blender library tree. These commands assume fresh `build/host-headers/python3.13`
and `build/blender-linux-sdk` directories:

```sh
mkdir -p build
curl -fL -o build/Python-3.13.13.tgz https://www.python.org/ftp/python/3.13.13/Python-3.13.13.tgz
printf '%s  %s\n' f9cde7b0e2ec8165d7326e2a0f59ea2686ce9d0c617dbbb3d66a7e54d31b74b9 build/Python-3.13.13.tgz | sha256sum -c -
tar -xzf build/Python-3.13.13.tgz -C build Python-3.13.13/Include
mkdir -p build/host-headers
cp -R build/Python-3.13.13/Include build/host-headers/python3.13
git init build/blender-linux-sdk
git -C build/blender-linux-sdk remote add origin https://projects.blender.org/blender/lib-linux_x64.git
git -C build/blender-linux-sdk config remote.origin.promisor true
git -C build/blender-linux-sdk config remote.origin.partialclonefilter blob:none
git -C build/blender-linux-sdk fetch --depth=1 --filter=blob:none origin 2b104b7a420a11f67fbad649cf3ea68156f8c9b9
git -C build/blender-linux-sdk show 2b104b7a420a11f67fbad649cf3ea68156f8c9b9:python/include/python3.13/pyconfig.h > build/host-headers/python3.13/pyconfig.h
printf '%s  %s\n' 97c991e1eed238ab1b61bda4813f72687f2e35bcf521f6c73d191633633b6f2a build/host-headers/python3.13/pyconfig.h | sha256sum -c -
```

Use the **actual paths** of the target Blender and its bundled interpreter throughout;
the packager checks the CMake cache against the probe. For an extracted source package,
extract the accompanying third-party source archive beside it as in the Mac recipe:

```sh
BLENDER=/opt/blender/blender
BPY=/opt/blender/5.2/python/bin/python3.13
python3 tools/probe_blender.py --blender "$BLENDER" --output build/host.json
cmake -S . -B build/linux-release -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCS_BUILD_BINDINGS=ON -DCS_BUILD_TESTS=OFF \
  -DCS_DOWNLOAD_CACHE="$PWD/../third-party-sources" \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON '-DCMAKE_INSTALL_RPATH=$ORIGIN' \
  -DPython_EXECUTABLE="$BPY" -DPython_INCLUDE_DIR="$PWD/build/host-headers/python3.13"
cmake --build build/linux-release --parallel 2
python3 tools/package_extension.py --blender "$BLENDER" --host build/host.json \
  --build build/linux-release --output dist/linux

```

An extracted source rebuild reports an unversioned native revision.
