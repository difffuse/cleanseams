# SPDX-License-Identifier: GPL-3.0-or-later
"""Locate and load the native engine (cleanseams_native) for this Blender's Python.

There is no Python fallback: when the module is missing or built for another Python ABI
the extension still registers, reports the diagnostic and disables its operators.
Search order: the directory in CLEANSEAMS_NATIVE_DIR (development and tests), then the
package's own `native/` directory (the shipped extension).
"""

import importlib.machinery
import importlib.util
import os
from pathlib import Path
import sys

MODULE = "cleanseams_native"
# Every attribute the add-on calls unconditionally,
# so a stale module is refused here rather than failing later with an AttributeError.
REQUIRED = ("Session", "read_host_mesh", "host_geometry_fingerprint", "resolve_ui_options_json",
            "apply_ui_profile_json", "run_messages", "packing_messages",
            "painted_attribute_names", "state_digest", "check_stored_uv", "ContractError", "close_all_sessions",
            "live_sessions", "recommended_memory_budget", "option_registry_json")

_module = None
_diagnostic = "native engine not loaded"


class NativeUnavailable(RuntimeError):
    """The native engine cannot be used in this Blender; the message says why."""


def search_dirs():
    dirs = []
    configured = os.environ.get("CLEANSEAMS_NATIVE_DIR")
    if configured:
        dirs.append(Path(configured))
    dirs.append(Path(__file__).resolve().parent / "native")
    return dirs


def _candidates(directory, name):
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.iterdir() if p.name.startswith(name + ".") and p.suffix in (".so", ".pyd"))


def load(name=MODULE):
    """Import the engine from the first search directory holding it. Raises NativeUnavailable
    with the reason: missing module, or modules built for another ABI than this Python's."""
    suffixes = importlib.machinery.EXTENSION_SUFFIXES
    seen = []
    for directory in search_dirs():
        found = _candidates(directory, name)
        seen += found
        for path in found:
            if not any(path.name == name + suffix for suffix in suffixes):
                continue
            # Keep imports owned by the addon package. The compiled binding has
            # module-local Python types and an embedded private core; both are
            # required, since a Python alias alone cannot isolate a shared ABI.
            fullname = f"{__package__ or __name__}.{name}"
            module = sys.modules.get(fullname)
            if module is not None:
                if Path(getattr(module, "__file__", "")).resolve() != path.resolve():
                    raise NativeUnavailable(f"{fullname} already belongs to another native path; restart Blender")
            else:
                try:
                    spec = importlib.util.spec_from_file_location(fullname, path)
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                except (ImportError, OSError) as exc:
                    raise NativeUnavailable(f"{path.name} could not be loaded: {exc}") from exc
            missing = [symbol for symbol in REQUIRED if not hasattr(module, symbol)]
            if missing:
                raise NativeUnavailable(f"{path.name} lacks {', '.join(missing)}; rebuild the native engine")
            sys.modules[fullname] = module
            return module
    if seen:
        tags = ", ".join(p.name for p in seen)
        raise NativeUnavailable(f"native engine built for another Python ({tags}); this Blender needs "
                                f"{name}{suffixes[0]}")
    raise NativeUnavailable(f"native engine {name} not found in " + ", ".join(str(d) for d in search_dirs()))


def activate():
    """Load once at registration; remember the diagnostic instead of raising."""
    global _module, _diagnostic
    try:
        _module = load()
        _diagnostic = ""
    except NativeUnavailable as exc:
        _module = None
        _diagnostic = str(exc)
    return _module


def deactivate():
    """Close every session of the engine (cancel, wake and join its worker) before the
    extension lets go of the module: no native worker outlives its owner. The module also
    closes them itself at interpreter exit."""
    global _module, _diagnostic
    module, _module = _module, None
    if module is not None and hasattr(module, "close_all_sessions"):
        module.close_all_sessions()
    _diagnostic = "native engine not loaded"


def get():
    """The loaded engine, or NativeUnavailable with the registration diagnostic."""
    if _module is None:
        raise NativeUnavailable(_diagnostic)
    return _module


def available():
    return _module is not None


def diagnostic():
    return _diagnostic
