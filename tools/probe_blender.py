#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Read the facts a native module and extension archive must match from an actual Blender.

Runs the given Blender executable with `--factory-startup --background` and throwaway
BLENDER_USER_* directories (never the user's configuration), and records: Blender version
and build, the bundled Python's version, SOABI, extension suffixes, free-threading/debug
flags, architecture and platform, NumPy version and buffer/foreach capabilities, the
extension repositories and the `--command extension` subcommands this Blender offers,
plus the executable's architectures and minimum OS read from the bundle.

    python tools/probe_blender.py --blender /path/to/blender --output build/host.json

Nothing is guessed from the version number: every field comes from the running Blender or
its files. The same module is executed inside Blender (`--python probe_blender.py -- --inside`).
"""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

SCHEMA = "cleanseams.blender_host"

# Blender's platform identifiers (blender_manifest.toml `platforms`) by (sys.platform, machine).
BLENDER_PLATFORMS = {("darwin", "arm64"): "macos-arm64", ("darwin", "x86_64"): "macos-x64",
                     ("linux", "x86_64"): "linux-x64", ("win32", "AMD64"): "windows-x64",
                     ("win32", "ARM64"): "windows-arm64"}


# --------------------------------------------------------------------------- inside Blender
def _inside(output):
    import importlib.machinery
    import struct
    import sysconfig

    import bpy
    facts = {"bpy": {"version": list(bpy.app.version), "version_string": bpy.app.version_string,
                     "build_hash": bpy.app.build_hash.decode() if isinstance(bpy.app.build_hash, bytes) else str(bpy.app.build_hash),
                     "build_platform": bpy.app.build_platform.decode() if isinstance(bpy.app.build_platform, bytes) else str(bpy.app.build_platform),
                     "build_type": bpy.app.build_type.decode() if isinstance(bpy.app.build_type, bytes) else str(bpy.app.build_type),
                     "binary_path": bpy.app.binary_path, "background": bpy.app.background,
                     "factory_startup": bpy.app.factory_startup, "online_access": bpy.app.online_access}}
    config = sysconfig.get_config_vars()
    facts["python"] = {
        "version": sys.version.split()[0], "version_info": list(sys.version_info[:3]), "full": sys.version,
        "implementation": sys.implementation.name, "cache_tag": sys.implementation.cache_tag,
        "abiflags": getattr(sys, "abiflags", ""), "executable": sys.executable, "prefix": sys.prefix,
        "SOABI": config.get("SOABI"), "EXT_SUFFIX": config.get("EXT_SUFFIX"),
        "extension_suffixes": list(importlib.machinery.EXTENSION_SUFFIXES),
        "Py_GIL_DISABLED": bool(config.get("Py_GIL_DISABLED")), "Py_DEBUG": bool(config.get("Py_DEBUG")),
        "Py_ENABLE_SHARED": bool(config.get("Py_ENABLE_SHARED")), "LDLIBRARY": config.get("LDLIBRARY"),
        "MACOSX_DEPLOYMENT_TARGET": config.get("MACOSX_DEPLOYMENT_TARGET"),
        "sysconfig_platform": sysconfig.get_platform(), "include": sysconfig.get_paths().get("include"),
        "pointer_bits": struct.calcsize("P") * 8, "byteorder": sys.byteorder,
        "platform": sys.platform, "machine": platform.machine(),
        "mac_ver": platform.mac_ver()[0], "platform_string": platform.platform()}
    try:  # the wheel tags this interpreter accepts (pip's vendored packaging when present)
        try:
            from packaging import tags
        except ImportError:
            from pip._vendor.packaging import tags
        facts["python"]["supported_tags"] = [str(t) for t in tags.sys_tags()][:400]
    except Exception as exc:  # noqa: BLE001 - recorded, not fatal
        facts["python"]["supported_tags_error"] = repr(exc)
    import numpy as np
    np_facts = {"version": np.__version__, "file": np.__file__}
    try:
        from numpy._core import _multiarray_umath as mu
    except ImportError:  # NumPy 1.x
        from numpy.core import _multiarray_umath as mu
    np_facts["c_api_version"] = hex(mu._get_ndarray_c_version())
    # Buffer protocol and the dtypes the host glue moves through foreach_get/foreach_set.
    np_facts["buffer_formats"] = {name: memoryview(np.zeros(2, dtype=name)).format
                                  for name in ("float32", "float64", "int32", "int64", "bool")}
    mesh = bpy.data.meshes.new("cleanseams_probe")
    mesh.from_pydata([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], [], [(0, 1, 2, 3)])
    mesh.update()
    foreach = {}
    for label, collection, prop, dtype, width in (("co_float32", mesh.vertices, "co", np.float32, 3),
                                                  ("co_float64", mesh.vertices, "co", np.float64, 3),
                                                  ("loop_vertex_int32", mesh.loops, "vertex_index", np.int32, 1),
                                                  ("edge_vertices_int32", mesh.edges, "vertices", np.int32, 2)):
        values = np.empty(len(collection) * width, dtype=dtype)
        try:
            collection.foreach_get(prop, values)
            foreach[label] = bool(np.isfinite(values.astype(np.float64)).all())
        except Exception as exc:  # noqa: BLE001
            foreach[label] = repr(exc)
    uv = mesh.uv_layers.new(name="probe")
    values = np.empty(len(mesh.loops) * 2, dtype=np.float64)
    uv.uv.foreach_get("vector", values)
    foreach["uv_float64"] = True
    bpy.data.meshes.remove(mesh)
    np_facts["foreach_get"] = foreach
    facts["numpy"] = np_facts
    facts["extensions"] = {
        "user_extensions": bpy.utils.user_resource("EXTENSIONS"),
        "user_config": bpy.utils.user_resource("CONFIG"),
        "user_scripts": bpy.utils.user_resource("SCRIPTS"),
        "repos": [{"module": r.module, "directory": r.directory, "remote": bool(r.use_remote_url), "enabled": bool(r.enabled)}
                  for r in bpy.context.preferences.extensions.repos]}
    Path(output).write_text(json.dumps(facts, indent=1))


# -------------------------------------------------------------------------- outside Blender
def temp_user_env(root):
    """A copy of the environment whose Blender user directories all live under `root`."""
    root = Path(root).resolve()
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP", "CLEANSEAMS_NATIVE_DIR")}
    for name in ("CONFIG", "SCRIPTS", "DATAFILES", "EXTENSIONS"):
        path = root / name.lower()
        path.mkdir(parents=True, exist_ok=True)
        env["BLENDER_USER_" + name] = str(path)
    (root / "tmp").mkdir(parents=True, exist_ok=True)
    env["TMPDIR"] = str(root / "tmp")
    env["PYTHONNOUSERSITE"] = "1"
    return env


def check_user_dirs_isolated(facts, root):
    """Refuse to continue when Blender resolved a user directory outside `root`."""
    root = Path(root).resolve()
    for key in ("user_extensions", "user_config", "user_scripts"):
        path = Path(facts["extensions"][key]).resolve()
        if root not in path.parents and path != root:
            raise SystemExit(f"Blender user directory {path} is outside the temporary root {root}; refusing")


def blender_command(blender, env, args, timeout=600):
    return subprocess.run([str(blender), "--factory-startup", "--background", *args], env=env,
                          capture_output=True, text=True, timeout=timeout)


def bundle_facts(blender):
    """Architectures and bundle metadata of the executable (macOS bundle when present)."""
    exe = Path(blender).resolve()
    out = {"executable": str(exe), "bytes": exe.stat().st_size}
    if sys.platform == "darwin" and shutil.which("lipo"):
        archs = subprocess.run(["lipo", "-archs", str(exe)], capture_output=True, text=True)
        out["architectures"] = archs.stdout.split()
        plist = exe.parent.parent / "Info.plist"
        if plist.exists():
            import plistlib
            info = plistlib.loads(plist.read_bytes())
            out["bundle_version"] = info.get("CFBundleShortVersionString")
            out["minimum_os"] = info.get("LSMinimumSystemVersion")
    return out


def extension_subcommands(blender, env):
    """`--command extension` subcommands offered by this Blender (discovered, not assumed)."""
    run = blender_command(blender, env, ["--command", "extension", "--help"], timeout=300)
    match = re.search(r"\{([a-z,\-]+)\}", run.stdout)
    return {"exit": run.returncode, "subcommands": match.group(1).split(",") if match else []}


def probe(blender, user_root):
    user_root = Path(user_root)
    if user_root.exists():
        shutil.rmtree(user_root)
    env = temp_user_env(user_root)
    output = user_root / "probe.json"
    run = blender_command(blender, env, ["--python-exit-code", "1", "--python", str(Path(__file__).resolve()),
                                         "--", "--inside", str(output)])
    if run.returncode != 0 or not output.exists():
        raise SystemExit(f"Blender probe failed (exit {run.returncode}):\n{run.stdout[-4000:]}\n{run.stderr[-4000:]}")
    facts = json.loads(output.read_text())
    check_user_dirs_isolated(facts, user_root)
    py = facts["python"]
    facts["schema"] = SCHEMA
    facts["probed_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    facts["bundle"] = bundle_facts(blender)
    facts["extension_command"] = extension_subcommands(blender, env)
    facts["blender_platform"] = BLENDER_PLATFORMS.get((py["platform"], py["machine"]))
    facts["probe_host"] = {"platform": platform.platform(), "machine": platform.machine(),
                           "mac_ver": platform.mac_ver()[0]}
    facts["verified_by"] = ("Blender --factory-startup --background with BLENDER_USER_* under " + str(user_root.resolve()) +
                            "; bpy.app, sys, sysconfig, importlib.machinery, platform, numpy and foreach_get read in-process")
    return facts


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--" in argv:  # running inside Blender: arguments follow `--`
        inner = argv[argv.index("--") + 1:]
        if inner and inner[0] == "--inside":
            _inside(inner[1])
            return 0
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--blender", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--user-root", type=Path, help="throwaway Blender user directory (default: next to --output)")
    args = parser.parse_args(argv)
    user_root = args.user_root or args.output.parent / "blender-probe-user"
    facts = probe(args.blender, user_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(facts, indent=1) + "\n")
    py = facts["python"]
    print(f"Blender {facts['bpy']['version_string']} Python {py['version']} {py['SOABI']} {py['machine']} "
          f"({facts['blender_platform']}) NumPy {facts['numpy']['version']} -> {args.output}")
    return 0


if __name__ == "__main__":
    if "bpy" in sys.modules:  # executed by Blender's --python
        main(sys.argv[sys.argv.index("--"):] if "--" in sys.argv else [])
    else:
        raise SystemExit(main())
