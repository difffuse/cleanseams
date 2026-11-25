#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build CleanSeams candidate artifacts for one probed Blender platform/ABI.

    python tools/package_extension.py --host build/host.json \
        --build build/blender-python-release --output dist [--source-package]

`--host` is the probe of an actual Blender (tools/probe_blender.py, or pass `--blender` to
probe here). `--build` is a CMake build of this project with CS_BUILD_BINDINGS=ON for that
Blender's Python, a loader-relative rpath (CMAKE_BUILD_WITH_INSTALL_RPATH=ON,
CMAKE_INSTALL_RPATH=@loader_path on Mac or $ORIGIN on Linux) and an optimized build type.
Linux also requires the retained target header receipt and the pinned auditwheel tools.
Installed smoke and exact-source release verification are separate. Produced in --output:

* cleanseams-<version>-<platform>.zip - the Blender extension: the host glue
  (blender/cleanseams/*.py), the manifest, native/ (the self-contained cleanseams_native
  binding for Blender's SOABI), COPYING and LICENSES/. Built and validated by Blender's own
  `--command extension build/validate`.
* cleanseams_native-<version>-<python>-<abi>-<platform>.whl - the native module as a wheel for
  the same interpreter ABI (tag checked against the tags that interpreter accepts).
* cleanseams_cli-<version>-<platform>.tar.gz - the command-line engine and core library.
* cleanseams-<version>-source.tar.gz and cleanseams-<version>-third-party-sources.tar
  (--source-package): the corresponding source (tracked project files without the frozen
  private fixture data or generated build outputs) and the verified dependency archives.
* artifacts.json - hashes, sizes, members, platform/ABI and binary facts of every artifact.

Refused rather than guessed: a module built for another Python, another architecture, a
non-relocatable rpath, a dependency on libpython or on build-tree paths, a minimum OS newer
than Blender's, a platform other than the probed one, or a debug build. The frozen Python
test data, tests and build tools never enter a runtime artifact.
"""
import argparse
import base64
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe_blender  # noqa: E402

SCHEMA = "cleanseams.artifacts"
# Platforms enabled for the complete artifact pipeline. Standalone inspection of another
# binary format does not enable its packaging or advertise it in the manifest.
INSPECTABLE_PLATFORMS = {"macos-arm64", "linux-x64"}
GLUE_FILES = ("__init__.py", "compat.py", "host.py", "native_module.py", "ops.py", "props.py", "ui.py")
# System libraries a shipped macOS binary may link besides libcleanseams_core.
MACOS_SYSTEM_LIBS = ("/usr/lib/libc++.1.dylib", "/usr/lib/libSystem.B.dylib")
CORE_INSTALL_NAME = "@rpath/libcleanseams_core.dylib"
LINUX_SYSTEM_LIBS = {"libc.so.6", "libm.so.6", "libstdc++.so.6", "libgcc_s.so.1",
                     "libpthread.so.0", "libdl.so.2", "librt.so.1", "ld-linux-x86-64.so.2"}
LINUX_CORE_SONAME = "libcleanseams_core.so"
LINUX_WHEEL_POLICY = "manylinux_2_35_x86_64"
# Deterministic archive timestamps (SOURCE_DATE_EPOCH when given).
EPOCH = int(os.environ.get("SOURCE_DATE_EPOCH", "315532800"))  # 1980-01-01


class PackagingError(RuntimeError):
    pass


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(argv, **kw):
    return subprocess.run([str(a) for a in argv], capture_output=True, text=True, **kw)


def read_manifest(path):
    import tomllib
    return tomllib.loads(Path(path).read_text())


def cmake_cache(build):
    cache = {}
    path = Path(build) / "CMakeCache.txt"
    if not path.exists():
        raise PackagingError(f"{build} is not a configured CMake build (no CMakeCache.txt)")
    for line in path.read_text().splitlines():
        m = re.match(r"^([A-Za-z_0-9]+):[A-Z_]+=(.*)$", line)
        if m:
            cache[m.group(1)] = m.group(2)
    return cache


# ------------------------------------------------------------------ binary inspection
def macho_facts(path):
    """Architectures, linked libraries, rpaths, minimum OS and signature of a Mach-O file."""
    path = Path(path)
    archs = run(["lipo", "-archs", path]).stdout.split()
    libs_out = run(["otool", "-L", path]).stdout.splitlines()[1:]
    libs = [line.strip().split(" (compatibility")[0] for line in libs_out if line.strip()]
    load = run(["otool", "-l", path]).stdout
    rpaths = re.findall(r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (\S+)", load)
    install_name = None
    m = re.search(r"cmd LC_ID_DYLIB\n\s+cmdsize \d+\n\s+name (\S+)", load)
    if m:
        install_name = m.group(1)
        libs = [lib for lib in libs if lib != install_name]
    minos = re.findall(r"cmd LC_BUILD_VERSION\n(?:.*\n){0,3}?\s+minos (\S+)", load)
    signature = run(["codesign", "--verify", "--strict", path])
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path), "architectures": archs,
            "linked": libs, "rpaths": rpaths, "install_name": install_name, "minos": minos[0] if minos else None,
            "signature_valid": signature.returncode == 0}


def version_tuple(text):
    return tuple(int(x) for x in str(text).split("."))


def elf_facts(path):
    """Inspect a dynamic ELF with GNU readelf; no binary code is executed.

    Numeric symbol-version maxima are recorded alongside the complete requirements.
    They are not, by themselves, a manylinux policy or host-compatibility certificate.
    """
    path = Path(path).resolve()
    try:
        result = run(["readelf", "--wide", "--file-header", "--program-headers", "--dynamic",
                      "--version-info", "--dyn-syms", path], env=dict(os.environ, LC_ALL="C"), timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PackagingError(f"cannot inspect ELF {path.name}: {exc}") from exc
    if result.returncode:
        raise PackagingError(f"readelf failed for {path.name}: {result.stderr[-2000:]}")
    output = result.stdout
    if any(section not in output for section in ("ELF Header:", "Program Headers:", "Symbol table '.dynsym' contains")) \
            or not any(section in output for section in ("Version symbols section", "No version information found")):
        raise PackagingError(f"incomplete ELF inspection for {path.name}: header, symbols or versions missing")
    def header(label):
        match = re.search(r"^\s*" + re.escape(label) + r":\s*(\S[^\n]*)$", output, re.MULTILINE)
        if not match:
            raise PackagingError(f"incomplete ELF inspection for {path.name}: missing {label}")
        return match.group(1).strip()
    machine = header("Machine")
    arch = {"Advanced Micro Devices X86-64": "x86_64", "AArch64": "aarch64", "Intel 80386": "i386"}.get(machine, machine)
    elf_class, data, elf_type = header("Class"), header("Data"), header("Type").split()[0]
    if "Dynamic section at offset" not in output:
        raise PackagingError(f"{path.name}: no inspected dynamic ELF section")
    linked = re.findall(r"\(NEEDED\)[^\n]*\[([^\]]+)\]", output)
    names = re.findall(r"\(SONAME\)[^\n]*\[([^\]]+)\]", output)
    if len(names) > 1:
        raise PackagingError(f"{path.name}: multiple ELF SONAME entries")
    paths = re.findall(r"\((?:RPATH|RUNPATH)\)[^\n]*\[([^\]]*)\]", output)
    interpreter = re.findall(r"\[Requesting program interpreter: ([^\]]+)\]", output)
    if len(interpreter) > 1:
        raise PackagingError(f"{path.name}: multiple ELF interpreters")
    exports = set()
    for line in output.splitlines():
        fields = line.split()
        if len(fields) >= 8 and fields[0].endswith(":") and fields[0][:-1].isdigit() \
                and fields[4] in ("GLOBAL", "WEAK", "UNIQUE") and fields[5] in ("DEFAULT", "PROTECTED") \
                and fields[6] != "UND":
            exports.add(fields[7].split("@")[0])
    requirements = {}
    library = None
    for line in output.partition("Version needs section")[2].splitlines():
        if line.startswith("Version ") and " section " in line:
            break
        match = re.search(r"\bFile:\s+(\S+)", line)
        if match:
            library = match.group(1)
            requirements.setdefault(library, set())
        match = re.search(r"\bName:\s+(\S+)", line)
        if match:
            if library is None:
                raise PackagingError(f"{path.name}: symbol version has no dependency")
            requirements[library].add(match.group(1))
    floors = {}
    for versions in requirements.values():
        for version in versions:
            match = re.fullmatch(r"(GLIBC|GLIBCXX|CXXABI|GCC)_([0-9]+(?:\.[0-9]+)*)", version)
            if match:
                family, number = match.groups()
                if family not in floors or version_tuple(number) > version_tuple(floors[family]):
                    floors[family] = number
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path), "format": "ELF",
            "architectures": [arch], "elf_class": elf_class, "byte_order": data, "elf_type": elf_type,
            "linked": linked, "rpaths": [part for value in paths for part in value.split(":")],
            "install_name": names[0] if names else None, "interpreter": interpreter[0] if interpreter else None,
            "exports": sorted(exports), "required_symbol_versions": {k: sorted(v) for k, v in sorted(requirements.items())},
            "symbol_version_floors": floors}


def verify_elf_binaries(files, host):
    """Verify the Linux x64 binary layout; staged ABI/import tests remain separate."""
    if host.get("blender_platform") != "linux-x64" or host["python"].get("machine") != "x86_64":
        raise PackagingError("ELF inspection supports only the probed Linux x86-64 target")
    facts = {key: elf_facts(files[key]) for key in ("core", "module", "cli")}
    problems = []
    for key, info in facts.items():
        name = info["file"]
        if info["architectures"] != ["x86_64"] or info["elf_class"] != "ELF64" \
                or info["byte_order"] != "2's complement, little endian":
            problems.append(f"{name}: wrong ELF architecture/class/byte order")
        if info["elf_type"] not in (("EXEC", "DYN") if key == "cli" else ("DYN",)):
            problems.append(f"{name}: invalid ELF type {info['elf_type']}")
        allowed = LINUX_SYSTEM_LIBS | ({LINUX_CORE_SONAME} if key == "cli" else set())
        bad = [lib for lib in info["linked"] if lib not in allowed]
        if bad:
            problems.append(f"{name}: unapproved dynamic dependencies {bad}")
        if key == "cli" and LINUX_CORE_SONAME not in info["linked"]:
            problems.append(f"{name}: does not link the shipped {LINUX_CORE_SONAME}")
        if info["rpaths"] != ["$ORIGIN"] and not (key == "core" and not info["rpaths"]):
            problems.append(f"{name}: non-relocatable or unexpected rpaths {info['rpaths']}")
        expected_interpreter = "/lib64/ld-linux-x86-64.so.2" if key == "cli" else None
        if info["interpreter"] != expected_interpreter:
            problems.append(f"{name}: unexpected ELF interpreter {info['interpreter']}")
        if any("GLIBC_PRIVATE" in versions for versions in info["required_symbol_versions"].values()):
            problems.append(f"{name}: depends on private glibc ABI")
    if facts["core"]["install_name"] != LINUX_CORE_SONAME:
        problems.append(f"core SONAME {facts['core']['install_name']} != {LINUX_CORE_SONAME}")
    private = [name for name in facts["core"]["exports"] if not name.startswith("_Z") and name not in ("_init", "_fini")]
    if private:
        problems.append(f"private C exports from core: {private}")
    if "PyInit_cleanseams_native" not in facts["module"]["exports"]:
        problems.append("module does not export PyInit_cleanseams_native")
    if problems:
        raise PackagingError("; ".join(problems))
    return facts


def check_binary(facts, arch, max_os, allowed_links, rpaths_expected, problems):
    name = facts["file"]
    if facts["architectures"] != [arch]:
        problems.append(f"{name}: architectures {facts['architectures']} != [{arch}]")
    bad = [lib for lib in facts["linked"] if lib not in allowed_links]
    if bad:
        problems.append(f"{name}: links {bad} (allowed: {sorted(allowed_links)})")
    if any("python" in lib.lower() for lib in facts["linked"]):
        problems.append(f"{name}: links a Python library")
    if facts["rpaths"] != rpaths_expected:
        problems.append(f"{name}: rpaths {facts['rpaths']} != {rpaths_expected} (build with "
                        "CMAKE_BUILD_WITH_INSTALL_RPATH=ON CMAKE_INSTALL_RPATH=@loader_path)")
    if facts["minos"] is None or version_tuple(facts["minos"]) > version_tuple(max_os):
        problems.append(f"{name}: minimum macOS {facts['minos']} is newer than the verified host {max_os}")
    if not facts["signature_valid"]:
        problems.append(f"{name}: code signature does not verify")


def wheel_platform_tag(platform_name, minos, arch):
    if platform_name == "linux-x64" and arch == "x86_64":
        # Intermediate tag only. qualify_linux_wheel must audit before publication.
        return "linux_x86_64"
    if platform_name.startswith("macos"):
        major, minor = (version_tuple(minos) + (0,))[:2]
        # From macOS 11 wheel tags carry only the major version: a binary needing 13.3 is
        # tagged 14_0 so that no interpreter on 13.0-13.2 accepts it.
        if major >= 11:
            return f"macosx_{major + (1 if minor else 0)}_0_{arch}"
        return f"macosx_{major}_{minor}_{arch}"
    raise PackagingError(f"no wheel platform tag rule for {platform_name}")


# ------------------------------------------------------------------ archives
def zip_deterministic(path, members):
    """members: list of (archive name, bytes, mode)."""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data, mode in sorted(members):
            info = zipfile.ZipInfo(name, date_time=time.gmtime(EPOCH)[:6])
            info.external_attr = (mode & 0xFFFF) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)


def tar_deterministic(path, members, compress=True):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for name, data, mode in sorted(members):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = EPOCH
            info.mode = mode
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    raw = buffer.getvalue()
    if compress:
        with open(path, "wb") as f, gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=EPOCH) as gz:
            gz.write(raw)
    else:
        Path(path).write_bytes(raw)


def archive_members(path):
    """(name, sha256, bytes) of every member of a zip, wheel or tar(.gz) artifact."""
    out = []
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                data = z.read(info)
                out.append({"name": info.filename, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
    else:
        with tarfile.open(path) as tar:
            for info in tar.getmembers():
                if info.isfile():
                    data = tar.extractfile(info).read()
                    out.append({"name": info.name, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
    return sorted(out, key=lambda m: m["name"])


def license_members(prefix):
    """COPYING plus LICENSES/** and the dependency lock/provenance, under `prefix`."""
    if not (ROOT / "LICENSES/host-header-provenance.json").is_file():
        raise PackagingError("missing build header provenance: LICENSES/host-header-provenance.json")
    members = [(prefix + "COPYING", (ROOT / "COPYING").read_bytes(), 0o644)]
    for path in sorted((ROOT / "LICENSES").rglob("*")):
        if path.is_file():
            members.append((prefix + "LICENSES/" + path.relative_to(ROOT / "LICENSES").as_posix(), path.read_bytes(), 0o644))
    members.append((prefix + "LICENSES/deps.lock.json", (ROOT / "deps.lock.json").read_bytes(), 0o644))
    return members


# ------------------------------------------------------------------ steps
def load_host(args, work):
    if args.host:
        host = json.loads(Path(args.host).read_text())
    else:
        if not args.blender:
            raise PackagingError("give --host (tools/probe_blender.py output) or --blender to probe")
        host = probe_blender.probe(args.blender, work / "probe-user")
        (work / "host.json").write_text(json.dumps(host, indent=1) + "\n")
    if host.get("schema") != probe_blender.SCHEMA:
        raise PackagingError("host facts are not a tools/probe_blender.py record")
    return host


def inspect_build(build, host, allow_debug):
    cache = cmake_cache(build)
    py = host["python"]
    problems = []
    if cache.get("CS_BUILD_BINDINGS", "").upper() not in ("ON", "TRUE", "1"):
        problems.append("build has CS_BUILD_BINDINGS off")
    if cache.get("CMAKE_BUILD_TYPE") not in ("Release", "RelWithDebInfo") and not allow_debug:
        problems.append(f"build type {cache.get('CMAKE_BUILD_TYPE')!r} is not an optimized release build")
    built_for = cache.get("Python_EXECUTABLE") or cache.get("_Python_EXECUTABLE")
    if not built_for or Path(built_for).resolve() != Path(py["executable"]).resolve():
        problems.append(f"module built for {built_for}, not Blender's Python {py['executable']}")
    module = Path(build) / ("cleanseams_native" + py["EXT_SUFFIX"])
    core_names = {"macos-arm64": "libcleanseams_core.dylib", "macos-x64": "libcleanseams_core.dylib",
                  "linux-x64": LINUX_CORE_SONAME}
    if host.get("blender_platform") not in core_names:
        raise PackagingError(f"build inspection does not support {host.get('blender_platform')}")
    core = Path(build) / core_names[host["blender_platform"]]
    cli = Path(build) / "cleanseams_cli"
    for path in (module, core, cli):
        if not path.exists():
            problems.append(f"missing {path} (a module named for SOABI {py['SOABI']} is required)")
    if problems:
        raise PackagingError("; ".join(problems))
    return cache, {"module": module, "core": core, "cli": cli}


def verify_host_headers(cache, host, receipt_path):
    """Match the retained target receipt to build headers and the actual interpreter."""
    receipt_path = Path(receipt_path)
    try:
        receipt = json.loads(receipt_path.read_text())
    except (OSError, ValueError) as exc:
        raise PackagingError(f"cannot read target header receipt: {exc}") from exc
    py = host["python"]
    version = ".".join(map(str, py["version_info"][:3]))
    if receipt.get("schema") != "cleanseams.host_headers" or type(receipt.get("schema_version")) is not int \
            or receipt["schema_version"] != 1:
        raise PackagingError("invalid target header receipt schema")
    if receipt.get("target") != host["blender_platform"] or receipt.get("version") != version \
            or receipt.get("blender_build_hash") != host["bpy"]["build_hash"]:
        raise PackagingError("header receipt target/version/Blender build differs from the probed host")
    headers = receipt.get("headers")
    if not isinstance(headers, dict) or not {"Python.h", "pyconfig.h", "patchlevel.h"}.issubset(headers):
        raise PackagingError("header receipt has no complete header inventory")
    include = cache.get("Python_INCLUDE_DIR") or cache.get("_Python_INCLUDE_DIR")
    if not include or not Path(include).is_dir():
        raise PackagingError(f"build Python header directory unavailable: {include}")
    directory = Path(include).resolve()
    for name, digest in headers.items():
        if not isinstance(name, str) or Path(name).is_absolute() or ".." in Path(name).parts \
                or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise PackagingError(f"invalid header receipt entry: {name}")
    actual = {p.relative_to(directory).as_posix(): p for p in directory.rglob("*") if p.is_file()}
    if set(actual) != set(headers):
        raise PackagingError("build header inventory differs from the target receipt")
    for name, path in actual.items():
        if path.is_symlink() or not path.resolve().is_relative_to(directory) or sha256_file(path) != headers[name]:
            raise PackagingError(f"build header changed or resolves outside the SDK: {name}")
    if receipt.get("pyconfig_sha256") != headers["pyconfig.h"]:
        raise PackagingError("pyconfig header hash differs from its receipt")
    patchlevel = (directory / "patchlevel.h").read_text()
    declared = re.search(r'^\s*#\s*define\s+PY_VERSION\s+"([^"]+)"', patchlevel, re.MULTILINE)
    if not declared or declared.group(1) != version:
        raise PackagingError("build header Python version differs from the probed interpreter")
    code = '''import json, platform, sys, sysconfig
with open(sys.argv[1]) as stream:
    headers = sysconfig.parse_config_h(stream)
runtime = sysconfig.get_config_vars()
print(json.dumps({'version': list(sys.version_info[:3]), 'machine': platform.machine(),
                  'soabi': runtime.get('SOABI'), 'macros': len(headers),
                  'differences': {k: [v, runtime.get(k)] for k, v in headers.items() if runtime.get(k) != v}}))
'''
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME") and not k.startswith(("LD_", "DYLD_"))}
    try:
        result = run([py["executable"], "-I", "-c", code, directory / "pyconfig.h"], env=env, timeout=60)
        if result.returncode:
            raise PackagingError(f"target header configuration probe failed: {result.stderr[-2000:]}")
        configuration = json.loads(result.stdout)
    except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
        raise PackagingError(f"cannot verify target header configuration: {exc}") from exc
    if configuration.get("version") != py["version_info"][:3] or configuration.get("machine") != py["machine"] \
            or configuration.get("soabi") != py["SOABI"] or not configuration.get("macros") \
            or configuration.get("differences"):
        raise PackagingError(f"build header configuration differs from target interpreter: {configuration}")
    return {"receipt": receipt_path.name, "receipt_sha256": sha256_file(receipt_path),
            "target": host["blender_platform"], "python_version": version,
            "include_directory": str(directory), "headers_verified": len(headers),
            "configuration_macros_verified": configuration["macros"]}


def verify_binaries(files, host):
    """Inspect the three binaries using their platform's loader/ABI rules.

    macOS additionally checks signatures and deployment floors against the probed host.
    ELF records complete symbol requirements; staged runtime checks remain mandatory.
    """
    if host.get("blender_platform", "").startswith("linux"):
        return verify_elf_binaries(files, host)
    arch = host["python"]["machine"]
    max_os = host["python"]["mac_ver"] or host["probe_host"]["mac_ver"]
    facts = {key: macho_facts(path) for key, path in files.items()}
    problems = []
    system = set(MACOS_SYSTEM_LIBS)
    # The core links nothing of ours; an @loader_path rpath on it is inert.
    check_binary(facts["core"], arch, max_os, system, facts["core"]["rpaths"] if facts["core"]["rpaths"] in ([], ["@loader_path"]) else [], problems)
    if facts["core"]["install_name"] != CORE_INSTALL_NAME:
        problems.append(f"core install name {facts['core']['install_name']} != {CORE_INSTALL_NAME}")
    check_binary(facts["module"], arch, max_os, system, ["@loader_path"], problems)
    check_binary(facts["cli"], arch, max_os, system | {CORE_INSTALL_NAME}, ["@loader_path"], problems)
    exported = run(["nm", "-gU", files["module"]]).stdout
    if "_PyInit_cleanseams_native" not in exported:
        problems.append("module does not export PyInit_cleanseams_native")
    if problems:
        raise PackagingError("; ".join(problems))
    return facts


def staged_cli_check(stage, py_facts):
    """The relocated CLI runs from the staged copy and loads the staged core library."""
    stage = Path(stage).resolve()
    if sys.platform not in ("darwin", "linux"):
        raise PackagingError(f"staged runtime inspection does not support {sys.platform}")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DYLD_", "LD_"))}
    trace_variable = "LD_DEBUG" if sys.platform == "linux" else "DYLD_PRINT_LIBRARIES"
    env[trace_variable] = "libs" if sys.platform == "linux" else "1"
    result = run([stage / "cleanseams_cli", "--version", "--json"], env=env, timeout=60)
    if result.returncode != 0:
        raise PackagingError(f"staged cleanseams_cli --version failed: {result.stderr[-2000:]}")
    if sys.platform == "linux":
        # glibc reports actual initialized images here, not just its search candidates.
        initialized = re.findall(r"^\s*\d+:\s+calling init:\s+(.+)$", result.stderr, re.MULTILINE)
        loaded = [p.strip() for p in initialized if Path(p.strip()).name == LINUX_CORE_SONAME]
    else:
        loaded = [line.split("> ")[-1].strip() for line in result.stderr.splitlines() if "libcleanseams_core" in line]
    if not loaded or any(Path(p).resolve().parent != stage.resolve() for p in loaded):
        raise PackagingError(f"staged CLI loaded the core from {loaded}, not from {stage}")
    version = json.loads(result.stdout)
    # The module imports in Blender's own interpreter (isolated mode) from the staged copy.
    # Ask the process's loader for its image list. Hardened macOS interpreters ignore
    # DYLD_*; Linux dl_iterate_phdr also preserves paths with spaces and Unicode.
    code = '''import ctypes, json, sys
sys.path.insert(0, sys.argv[1])
import cleanseams_native as m
d = ctypes.CDLL(None)
if sys.platform == "linux":
    class ImageInfo(ctypes.Structure):
        _fields_ = [("address", ctypes.c_void_p), ("name", ctypes.c_char_p)]
    callback_type = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ImageInfo), ctypes.c_size_t, ctypes.c_void_p)
    images = []
    @callback_type
    def collect(info, size, data):
        name = info.contents.name
        if name:
            images.append(name.decode(sys.getfilesystemencoding(), errors="surrogateescape"))
        return 0
    d.dl_iterate_phdr.argtypes = [callback_type, ctypes.c_void_p]
    d.dl_iterate_phdr.restype = ctypes.c_int
    if d.dl_iterate_phdr(collect, None) != 0:
        raise RuntimeError("native loader image enumeration did not complete")
else:
    d._dyld_image_count.restype = ctypes.c_uint32
    d._dyld_get_image_name.restype = ctypes.c_char_p
    d._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    images = [d._dyld_get_image_name(i).decode() for i in range(d._dyld_image_count())]
print(json.dumps({'file': m.__file__, 'has_session': hasattr(m, 'Session'),
                  'core': [p for p in images if 'libcleanseams_core' in p]}))
'''
    env_py = {k: v for k, v in env.items() if k not in ("PYTHONPATH", "PYTHONHOME", trace_variable)}
    imported = run([py_facts["executable"], "-I", "-c", code, stage], env=env_py, timeout=120)
    if imported.returncode != 0:
        raise PackagingError(f"staged module does not import in Blender's Python: {imported.stderr[-2000:]}")
    module_import = json.loads(imported.stdout.splitlines()[-1])
    expected_module = stage / ("cleanseams_native" + py_facts["EXT_SUFFIX"])
    actual_module = Path(module_import["file"]).resolve()
    if actual_module != expected_module or not module_import.get("has_session"):
        raise PackagingError(f"native module {actual_module} is not from {expected_module} or has no Session API")
    if module_import["core"]:
        raise PackagingError(f"staged binding resolved a shared core: {module_import['core']}")
    return {"cli_version": version, "module_import": module_import, "core_loaded_from_stage": True,
            "binding_core_is_private": True}


def blender_extension(blender, host, stage_dir, output, work):
    """Validate and build the extension with this Blender's own extension tooling."""
    env = probe_blender.temp_user_env(work / "blender-build-user")
    commands = host["extension_command"]["subcommands"]
    for needed in ("build", "validate"):
        if needed not in commands:
            raise PackagingError(f"this Blender has no `--command extension {needed}`")
    log = {}
    validate = probe_blender.blender_command(blender, env, ["--command", "extension", "validate", str(stage_dir)])
    log["validate_source"] = {"exit": validate.returncode, "output": (validate.stdout + validate.stderr)[-3000:]}
    if validate.returncode != 0:
        raise PackagingError(f"extension validate failed: {log['validate_source']['output']}")
    out_dir = work / "blender-build-output"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    build = probe_blender.blender_command(blender, env, ["--command", "extension", "build", "--source-dir", str(stage_dir),
                                                         "--output-dir", str(out_dir), "--split-platforms"])
    log["build"] = {"exit": build.returncode, "output": (build.stdout + build.stderr)[-3000:]}
    if build.returncode != 0:
        raise PackagingError(f"extension build failed: {log['build']['output']}")
    built = sorted(out_dir.glob("*.zip"))
    platform_suffix = host["blender_platform"].replace("-", "_")
    matching = [p for p in built if p.stem.endswith(platform_suffix)]
    if len(built) != 1 or len(matching) != 1:
        raise PackagingError(f"expected one {platform_suffix} archive from extension build, got {[p.name for p in built]}")
    final = output / matching[0].name
    shutil.copyfile(matching[0], final)
    validate = probe_blender.blender_command(blender, env, ["--command", "extension", "validate", str(final)])
    log["validate_archive"] = {"exit": validate.returncode, "output": (validate.stdout + validate.stderr)[-3000:]}
    if validate.returncode != 0:
        raise PackagingError(f"built archive does not validate: {log['validate_archive']['output']}")
    return final, log


def target_manifest(template, host):
    target = host.get("blender_platform")
    if target not in INSPECTABLE_PLATFORMS:
        raise PackagingError(f"artifact assembly does not support {target}")
    return dict(template, platforms=[target])


def stage_extension(stage, files, host, manifest):
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    source = ROOT / "blender/cleanseams"
    present = sorted(p.name for p in source.iterdir() if p.is_file() and p.suffix == ".py")
    if tuple(present) != tuple(sorted(GLUE_FILES)):
        raise PackagingError(f"unexpected host glue file set {present}; update GLUE_FILES and the audit together")
    for name in GLUE_FILES + ("blender_manifest.toml",):
        shutil.copyfile(source / name, stage / name)
    manifest_path = stage / "blender_manifest.toml"
    if read_manifest(manifest_path).get("platforms") != manifest.get("platforms"):
        content, count = re.subn(r"(?m)^platforms[ \t]*=[ \t]*\[[^\]]*\]",
                                "platforms = " + json.dumps(manifest["platforms"]), manifest_path.read_text())
        if count != 1:
            raise PackagingError("source manifest must contain exactly one platform list")
        manifest_path.write_text(content)
    if read_manifest(manifest_path) != manifest:
        raise PackagingError("staged manifest differs from the selected artifact metadata")
    native = stage / "native"
    native.mkdir()
    shutil.copy2(files["module"], native / files["module"].name)
    for name, data, _mode in license_members(""):
        target = stage / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    platforms = manifest.get("platforms")
    if platforms != [host["blender_platform"]]:
        raise PackagingError(f"manifest platforms {platforms} must list exactly the probed platform "
                             f"{host['blender_platform']} (only verified platforms are advertised)")


def build_wheel(output, files, host, version, binary_facts):
    py = host["python"]
    cp = "cp" + "".join(str(x) for x in py["version_info"][:2])
    abi = cp + ("t" if py["Py_GIL_DISABLED"] else "") + ("d" if py["Py_DEBUG"] else "")
    plat = wheel_platform_tag(host["blender_platform"], (binary_facts.get("module") or {}).get("minos"), py["machine"])
    tag = f"{cp}-{abi}-{plat}"
    supported = py.get("supported_tags")
    if supported is not None and tag not in supported:
        raise PackagingError(f"wheel tag {tag} is not accepted by Blender's Python ({supported[:4]} ...)")
    dist = f"cleanseams_native-{version}.dist-info"
    members = [(files["module"].name, files["module"].read_bytes(), 0o755)]
    licenses = []
    for name, data, mode in license_members(""):
        members.append((f"{dist}/licenses/{name}", data, mode))
        licenses.append(name)
    metadata = ["Metadata-Version: 2.4", "Name: cleanseams_native", f"Version: {version}",
                "Summary: CleanSeams native seam and UV engine (self-contained Python binding)",
                "License-Expression: GPL-3.0-or-later", f"Requires-Python: =={py['version_info'][0]}.{py['version_info'][1]}.*",
                "Requires-Dist: numpy"]
    metadata += [f"License-File: {name}" for name in licenses]
    members.append((f"{dist}/METADATA", ("\n".join(metadata) + "\n").encode(), 0o644))
    wheel = ["Wheel-Version: 1.0", "Generator: cleanseams tools/package_extension.py", "Root-Is-Purelib: false", f"Tag: {tag}"]
    members.append((f"{dist}/WHEEL", ("\n".join(wheel) + "\n").encode(), 0o644))
    record = []
    for name, data, _mode in members:
        digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
        record.append(f"{name},sha256={digest},{len(data)}")
    record.append(f"{dist}/RECORD,,")
    members.append((f"{dist}/RECORD", ("\n".join(sorted(record)) + "\n").encode(), 0o644))
    path = output / f"cleanseams_native-{version}-{tag}.whl"
    zip_deterministic(path, members)
    return path, tag


def qualify_linux_wheel(wheel, output, host):
    """Require auditwheel policy/ISA checks; retain the already relocatable native bytes."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-m", "auditwheel", "repair", "--only-plat", "--plat", LINUX_WHEEL_POLICY,
               "--wheel-dir", str(output), str(wheel)]
    try:
        version = run([sys.executable, "-m", "auditwheel", "--version"], timeout=30)
        result = run(command, timeout=180)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PackagingError(f"auditwheel unavailable or timed out: {exc}") from exc
    if version.returncode or result.returncode:
        raise PackagingError(f"auditwheel refused {LINUX_WHEEL_POLICY}: {(result.stdout + result.stderr)[-4000:]}")
    wheels = sorted(output.glob("*.whl"))
    tag = wheel.name.rsplit("-", 1)[0].rsplit("-", 2)[1:]
    tag = "-".join([*tag, LINUX_WHEEL_POLICY])
    if len(wheels) != 1 or not wheels[0].name.endswith("-" + tag + ".whl"):
        raise PackagingError(f"auditwheel emitted unexpected artifacts: {[p.name for p in wheels]}")
    if tag not in host["python"].get("supported_tags", []):
        raise PackagingError(f"audited wheel tag {tag} is not accepted by the probed Blender interpreter")
    def native_members(path):
        with zipfile.ZipFile(path) as archive:
            return {n: hashlib.sha256(archive.read(n)).hexdigest() for n in archive.namelist()
                    if not n.endswith("/") and archive.read(n)[:4] == b"\x7fELF"}
    if not native_members(wheel) or native_members(wheel) != native_members(wheels[0]):
        raise PackagingError("auditwheel changed native members; review/requalify those binaries before packaging")
    with zipfile.ZipFile(wheels[0]) as archive:
        metadata = [n for n in archive.namelist() if n.endswith(".dist-info/WHEEL")]
        if len(metadata) != 1 or f"Tag: {tag}" not in archive.read(metadata[0]).decode().splitlines():
            raise PackagingError("auditwheel output metadata does not match its filename tag")
    return wheels[0], tag, {"policy": LINUX_WHEEL_POLICY, "tool": version.stdout.strip(), "command": command,
                           "exit": result.returncode, "output": (result.stdout + result.stderr)[-6000:],
                           "native_members_unchanged": True}


def build_cli_archive(output, files, host, version):
    top = f"cleanseams_cli-{version}-{host['blender_platform']}/"
    loader_path = "$ORIGIN" if host["blender_platform"] == "linux-x64" else "@loader_path"
    readme = (f"CleanSeams command-line engine {version} for {host['blender_platform']}.\n"
              f"Run ./cleanseams_cli --version --json. {files['core'].name} must stay next to the executable\n"
              f"({loader_path}). Licensed GPL-3.0-or-later (COPYING); third-party notices in LICENSES/.\n")
    members = [(top + "cleanseams_cli", files["cli"].read_bytes(), 0o755),
               (top + files["core"].name, files["core"].read_bytes(), 0o755),
               (top + "README.txt", readme.encode(), 0o644)]
    members += license_members(top)
    path = output / f"cleanseams_cli-{version}-{host['blender_platform']}.tar.gz"
    tar_deterministic(path, members)
    return path


# Generated output directories never enter the source archive.
SOURCE_EXCLUDE = ("build/", "dist/", "execution/")


def git_output(*args):
    """Require a successful Git query; empty output is not proof of a clean tree."""
    try:
        result = run(["git", "-C", ROOT, *args])
    except OSError as exc:
        raise PackagingError(f"source Git query unavailable: {exc}") from exc
    if result.returncode:
        raise PackagingError(f"source Git query failed: {' '.join(args)}")
    return result.stdout


def product_tree_changes():
    """Uncommitted changes to files that enter the product (evidence logs excluded)."""
    prefix = git_output("rev-parse", "--show-prefix").rstrip("\n")
    status = git_output("status", "--porcelain=v1", "--untracked-files=all", "-z", "--", ".")
    changed = set()
    entries = iter(status.split("\0"))
    for entry in entries:
        if not entry:
            continue
        paths = [entry[3:]]
        # Porcelain -z puts the destination first, then the original path without
        # a status prefix. Either side may change a retained product file.
        if "R" in entry[:2] or "C" in entry[:2]:
            original = next(entries, "")
            if not original:
                raise PackagingError("incomplete Git rename/copy status")
            paths.append(original)
        for path in paths:
            if path.startswith(prefix):
                rel = path[len(prefix):]
                if not rel.startswith(SOURCE_EXCLUDE):
                    changed.add(rel)
    return sorted(changed)


def tracked_product_files():
    """An enclosing repository does not supply history for an untracked extraction."""
    names = [n for n in git_output("ls-files", "-z", "--", ".").split("\0")
             if n and not n.startswith(SOURCE_EXCLUDE)]
    if not names:
        raise PackagingError("source Git query found no tracked product files")
    return names


def product_source_state():
    """Report Git facts or explicit unknowns, never infer them from an archive receipt."""
    try:
        commit = git_output("rev-parse", "--verify", "HEAD").strip()
        tracked_product_files()
        dirty = product_tree_changes()
    except PackagingError:
        return {"source_state": "unavailable", "commit": None, "working_tree_dirty": None,
                "uncommitted_product_files": None}
    return {"source_state": "git", "commit": commit, "working_tree_dirty": bool(dirty),
            "uncommitted_product_files": dirty}


def build_source_package(output, version, *, dependency_cache=None):
    names = tracked_product_files()
    dirty = product_tree_changes()
    commit = git_output("rev-parse", "--verify", "HEAD").strip()
    top = f"cleanseams-{version}-source/"
    members = []
    for name in names:
        path = ROOT / name
        if path.is_file():
            members.append((top + name, path.read_bytes(), 0o755 if os.access(path, os.X_OK) else 0o644))
    note = {"commit": commit, "working_tree_clean": not dirty, "excluded_prefixes": list(SOURCE_EXCLUDE),
            "dependencies": json.loads((ROOT / "deps.lock.json").read_text())["dependencies"]}
    members.append((top + "SOURCE_PACKAGE.json", (json.dumps(note, indent=1) + "\n").encode(), 0o644))
    source = output / f"cleanseams-{version}-source.tar.gz"
    tar_deterministic(source, members)
    lock = json.loads((ROOT / "deps.lock.json").read_text())
    cache = Path(dependency_cache) if dependency_cache is not None else ROOT / "build/downloads"
    deps = [(f"third-party-sources/deps.lock.json", (ROOT / "deps.lock.json").read_bytes(), 0o644)]
    for dep in lock["dependencies"]:
        archive = cache / f"{dep['name']}-{dep['sha256']}.tar.gz"
        if not archive.exists() or sha256_file(archive) != dep["sha256"]:
            raise PackagingError(f"verified dependency archive missing or changed: {archive}")
        deps.append((f"third-party-sources/{archive.name}", archive.read_bytes(), 0o644))
    third = output / f"cleanseams-{version}-third-party-sources.tar"
    tar_deterministic(third, deps, compress=False)
    return source, third, {"commit": commit, "working_tree_clean": not dirty, "files": len(members)}


def record(path, kind, output, **extra):
    return {"kind": kind, "path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path),
            "members": archive_members(path) if kind != "manifest" else [], **extra}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--host", type=Path, help="tools/probe_blender.py output for the target Blender")
    parser.add_argument("--blender", type=Path, help="Blender executable (probed when --host is absent; runs extension build/validate)")
    parser.add_argument("--build", type=Path, required=True, help="release bindings build for that Blender's Python")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, help="staging directory (default: <output>/.work)")
    parser.add_argument("--source-package", action="store_true", help="also write the source and third-party source packages")
    parser.add_argument("--allow-debug", action="store_true", help="accept a Debug build (never for release artifacts)")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    work = (args.work or output / ".work").resolve()
    work.mkdir(parents=True, exist_ok=True)
    if output.exists():
        for old in list(output.glob("*.zip")) + list(output.glob("*.whl")) + list(output.glob("*.tar*")) + [output / "artifacts.json"]:
            if old.exists():
                old.unlink()
    output.mkdir(parents=True, exist_ok=True)
    try:
        host = load_host(args, work)
        blender = args.blender or Path(host["bpy"]["binary_path"])
        expected_system = "linux" if host["blender_platform"] == "linux-x64" else "darwin"
        if sys.platform != expected_system or host["blender_platform"] not in INSPECTABLE_PLATFORMS:
            raise PackagingError(f"packaging for {host.get('blender_platform')} on {sys.platform} is not implemented or verified")
        manifest = target_manifest(read_manifest(ROOT / "blender/cleanseams/blender_manifest.toml"), host)
        version = manifest["version"]
        cache, files = inspect_build(args.build, host, args.allow_debug)
        headers = None
        if host["blender_platform"] == "linux-x64":
            headers = verify_host_headers(cache, host, ROOT / "LICENSES/linux-host-header-provenance.json")
        binary_facts = verify_binaries(files, host)
        cli_stage = work / "stage-native"
        if cli_stage.exists():
            shutil.rmtree(cli_stage)
        cli_stage.mkdir(parents=True)
        for path in files.values():
            shutil.copy2(path, cli_stage / path.name)
        staged = staged_cli_check(cli_stage, host["python"])
        ext_stage = work / "stage-extension"
        stage_extension(ext_stage, files, host, manifest)
        extension, blender_log = blender_extension(blender, host, ext_stage, output, work)
        wheel_audit = None
        if host["blender_platform"] == "linux-x64":
            wheel_input = work / "wheel-input"
            wheel_input.mkdir(exist_ok=True)
            provisional, _ = build_wheel(wheel_input, files, host, version, binary_facts)
            wheel, tag, wheel_audit = qualify_linux_wheel(provisional, work / "wheel-audited", host)
            shipped_wheel = output / wheel.name
            shutil.copy2(wheel, shipped_wheel)
            wheel = shipped_wheel
            os_floor = {"kind": "elf-symbol-versions",
                        "requirements": {k: v["required_symbol_versions"] for k, v in binary_facts.items()},
                        "verified_on": host["python"]["platform_string"],
                        "note": "symbol requirements and staged checks; installed smoke and native host qualification are separate"}
        else:
            wheel, tag = build_wheel(output, files, host, version, binary_facts)
            os_floor = {"artifact": max((f["minos"] for f in binary_facts.values()), key=version_tuple),
                        "blender": host["bundle"].get("minimum_os"), "verified_on": host["python"]["mac_ver"],
                        "note": "the artifacts need the newer of these floors; only verified_on was tested"}
        cli = build_cli_archive(output, files, host, version)
        platform_record = {"platform": host["blender_platform"], "python": host["python"]["version"],
                           "soabi": host["python"]["SOABI"], "ext_suffix": host["python"]["EXT_SUFFIX"],
                           "architecture": host["python"]["machine"]}
        artifacts = [record(extension, "blender_extension", output, **platform_record, blender_min=manifest["blender_version_min"]),
                     record(wheel, "native_wheel", output, **platform_record, wheel_tag=tag),
                     record(cli, "cli_archive", output, platform=host["blender_platform"], architecture=host["python"]["machine"])]
        source_info = None
        if args.source_package:
            source, third, source_info = build_source_package(output, version, dependency_cache=cache.get("CS_DOWNLOAD_CACHE"))
            artifacts += [record(source, "source_package", output), record(third, "third_party_sources", output)]
    except PackagingError as exc:
        print(f"package_extension: {exc}", file=sys.stderr)
        return 2
    summary = {"schema": SCHEMA, "version": 1, "product_version": version,
               "engine": staged["cli_version"], **product_source_state(),
               "host": {"blender": host["bpy"]["version_string"], "build_hash": host["bpy"]["build_hash"],
                        "executable": host["bundle"]["executable"], "python": host["python"]["version"],
                        "python_executable": host["python"]["executable"],
                        "soabi": host["python"]["SOABI"], "ext_suffix": host["python"]["EXT_SUFFIX"],
                        "architecture": host["python"]["machine"], "platform": host["blender_platform"],
                        "minimum_os": host["bundle"].get("minimum_os"), "numpy": host["numpy"]["version"],
                        "numpy_c_api": host["numpy"]["c_api_version"], "free_threaded": host["python"]["Py_GIL_DISABLED"]},
               "build": {"directory": str(Path(args.build).resolve()), "type": cache.get("CMAKE_BUILD_TYPE"),
                         "compiler": staged["cli_version"].get("compiler"),
                         "deployment_target": cache.get("CMAKE_OSX_DEPLOYMENT_TARGET"),
                         "deps_lock_sha256": sha256_file(ROOT / "deps.lock.json"), "header_provenance": headers},
               "binaries": binary_facts, "staged_checks": staged, "blender_tooling": blender_log,
               "os_floor": os_floor, "wheel_audit": wheel_audit,
               "manifest": {"id": manifest["id"], "version": version, "blender_version_min": manifest["blender_version_min"],
                            "platforms": manifest.get("platforms"), "license": manifest.get("license")},
               "source_package": source_info, "artifacts": artifacts}
    (output / "artifacts.json").write_text(json.dumps(summary, indent=1) + "\n")
    for a in artifacts:
        print(f"{a['kind']:20s} {a['path']}  {a['bytes']} bytes  sha256 {a['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
