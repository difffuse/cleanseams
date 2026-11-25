# Dependency provenance and notices

`../deps.lock.json` identifies exact upstream release tags, resolved commits, archive SHA-256 hashes and licenses. Sources are fetched only at configure and are excluded from Git/build products by default. `notices/` preserves upstream license texts; `../COPYING` preserves CleanSeams GPL-3.0-or-later licensing.

Upstream release references used to select the initial lock:
- Eigen 3.4.0: https://gitlab.com/libeigen/eigen/-/releases/3.4.0
- Serial SuperLU 7.0.0: https://github.com/xiaoyeli/superlu/releases/tag/v7.0.0
- pybind11 3.0.1: https://github.com/pybind/pybind11/releases/tag/v3.0.1
- nlohmann/json 3.12.0: https://github.com/nlohmann/json/releases/tag/v3.12.0
- Catch2 3.8.1: https://github.com/catchorg/Catch2/releases/tag/v3.8.1

Python SDK headers are part of the separately selected build environment. The local Blender SDK workaround uses the exact official CPython 3.13.13 release source plus that Blender installation's pyconfig.h; hashes are in LICENSES/host-header-provenance.json. See docs/BUILD.md for reproducible preparation. It is not a runtime dependency of the core or CLI.
