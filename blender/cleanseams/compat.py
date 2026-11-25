# SPDX-License-Identifier: GPL-3.0-or-later
"""Blender version differences, resolved through RNA rather than assumed."""

import bpy


def unwrap_methods():
    props = bpy.ops.uv.unwrap.get_rna_type().properties
    return [item.identifier for item in props["method"].enum_items]


def unwrap_method(requested=None):
    """The requested method when this Blender offers it, otherwise the best interior solver
    it has: minimum stretch (SLIM) from 4.3, conformal (LSCM) before."""
    items = unwrap_methods()
    if requested in items:
        return requested
    if "MINIMUM_STRETCH" in items:
        return "MINIMUM_STRETCH"
    return "CONFORMAL" if "CONFORMAL" in items else items[0]


def unwrap_kwargs(method):
    props = bpy.ops.uv.unwrap.get_rna_type().properties
    kw = {"method": method, "margin": 0.001}
    if "correct_aspect" in props:
        kw["correct_aspect"] = False
    if "fill_holes" in props:
        # The core has already opened every chart into a disk. Artificial hole
        # filling would introduce triangles outside that validated cut surface.
        kw["fill_holes"] = False
    return kw


def seam_attribute(me, create=False):
    attr = me.attributes.get("uv_seam")
    if attr is None and create:
        attr = me.attributes.new(name="uv_seam", type="BOOLEAN", domain="EDGE")
    return attr
