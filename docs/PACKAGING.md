# Packaging

Build from a clean distribution checkout. `tools/package_extension.py` packages an
optimized native build for a probed target Blender. Its outputs include the optional
standalone Blender extension, ABI-specific native wheel and CLI/core archive.
`--source-package` also emits the exact tracked source and verified dependency archives.
Preserve `COPYING`, `LICENSES/`, dependency sources and compiler/runtime notices.

```sh
python tools/package_extension.py --blender /path/to/blender \
  --build /path/to/native-build --output /path/to/new-dist --source-package
```

The packager validates architecture, ABI, loader paths, runtime floors, licenses and
archive membership. Current full extension/wheel inspection supports macOS ARM64 and
Linux x64. Windows SDK installation is available through CMake; Windows binary and
Blender acceptance must be performed on that target before claiming qualification.
Linux packaging requires the pinned tools in `tools/requirements-linux.txt` and
matching target header provenance.
