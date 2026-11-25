# SPDX-License-Identifier: GPL-3.0-or-later
"""cleanseams: seams that follow the quad flow, straight or smooth UV contours, and nothing
applied that does not pass validation.

This package is the Blender host glue of the native engine (cleanseams_native): settings,
panels, operators, bulk array transfer and Blender's unwrap on the isolated chart meshes the
engine requests. Every computation is native; there is no Python engine and no fallback.
"""

__version__ = "0.2.0"


def register():  # pragma: no cover - Blender only
    from . import native_module, props, ops, ui
    native_module.activate()  # a missing engine is reported, never replaced
    props.register()
    ops.register()
    ui.register()


def unregister():  # pragma: no cover - Blender only
    from . import native_module, props, ops, ui
    ui.unregister()
    ops.unregister()
    props.unregister()
    native_module.deactivate()
