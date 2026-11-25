# SPDX-License-Identifier: GPL-3.0-or-later
"""Scene settings of the add-on (RNA only). Profiles, budgets and percentages are resolved
by the native engine (`resolve_ui_options_json`, `apply_ui_profile_json`). The enhanced
controls (T41) are declared from the native option registry at registration: their kind,
label, description, range, enum items, default and panel section are native; this module
only turns them into RNA properties and transfers their values."""

import json

import bpy
from bpy.props import (BoolProperty, EnumProperty, FloatProperty, IntProperty, StringProperty,
                       PointerProperty, FloatVectorProperty)

from . import native_module

# The native UiSettings fields of the frozen add-on, in the order of the extension's properties.
UI_FIELDS = ("profile", "generate", "protect_existing", "feature_angle", "refine_paths", "contour_mode", "budget",
             "min_island", "solver_iterations", "min_island_faces", "target_enabled", "target_min", "target_max",
             "deadline_s", "hide_from", "hide_direction", "show_advanced", "w_turn", "w_flow", "w_drift",
             "corridor_rings", "smooth_radius", "max_displacement", "min_fusion", "pack_after", "keep_orientation",
             "pack_width", "pack_height", "pack_padding", "pack_border", "pack_tiles",
             "last_status", "last_report")

# T41: the native UiSettings fields declared from the registry at registration (none without
# the engine), and the native panel sections that draw them.
NATIVE_FIELDS = []
SECTIONS = []
PROPERTY_TYPES = {"BoolProperty": BoolProperty, "EnumProperty": EnumProperty, "FloatProperty": FloatProperty,
                  "IntProperty": IntProperty, "StringProperty": StringProperty, "FloatVectorProperty": FloatVectorProperty}
PROPERTY_KEYWORDS = ("name", "description", "min", "max", "subtype", "size")
# Host-only fields: never sent back by a profile update.
HOST_FIELDS = ("profile", "last_status", "last_report", "last_messages")


def native_declarations():
    """RNA properties of the native UI fields this file does not spell out, exactly as the
    native registry describes them, and the registry's panel sections."""
    registry = json.loads(native_module.get().option_registry_json())
    metadata, schema = registry["ui_metadata"], registry["ui"]["properties"]
    declared = {}
    for name in registry["ui_enhanced_fields"]:
        info = metadata[name]
        keywords = {key: info[key] for key in PROPERTY_KEYWORDS if key in info}
        if "items" in info:
            keywords["items"] = [tuple(item) for item in info["items"]]
        keywords["default"] = schema[name]["default"]
        declared[name] = PROPERTY_TYPES[info["property_type"]](**keywords)
    return declared, registry["ui_sections"]


def field_names():
    """Every UiSettings field this registration declares: the frozen add-on's, then the native ones."""
    return UI_FIELDS + tuple(NATIVE_FIELDS)


def settings_json(st):
    """The property values as the native UiSettings JSON (transport only)."""
    out = {}
    for name in field_names():
        value = getattr(st, name)
        out[name] = list(value) if name == "hide_direction" else value
    return json.dumps(out)


def _profile_update(self, context):
    if not native_module.available():
        return
    applied = json.loads(native_module.get().apply_ui_profile_json(settings_json(self), self.profile))
    for name in field_names():
        if name in HOST_FIELDS:
            continue
        value = applied[name]
        current = list(getattr(self, name)) if name == "hide_direction" else getattr(self, name)
        if current != value:
            setattr(self, name, value)


class CleanSeamsSettings(bpy.types.PropertyGroup):
    profile: EnumProperty(
        name="Profile",
        items=[('auto', 'Auto', 'Choose from geometric panels, revolved parts and curved surfaces'),
               ("hard_surface", "Hard surface", "Seams on hard edges, straight sides where "
                "the topology is a grid"),
               ("organic", "Organic", "No feature loops, smoothed contours that keep the "
                "curves"),
               ("custom", "Custom", "Whatever is set below")],
        default="hard_surface", update=_profile_update)
    generate: BoolProperty(
        name="Generate seams", default=True,
        description="Trace the quad layout and merge patches under the distortion budget. "
                    "Off: only refine the existing seams and clean the contours")
    protect_existing: BoolProperty(
        name="Protect existing seams", default=True,
        description="Seams already on the mesh are kept as locked cuts")
    feature_angle: FloatProperty(
        name="Feature angle", default=60.0, min=0.0, max=180.0,
        description="Edges sharper than this (degrees) become seams and are never merged "
                    "away; 0 disables")
    refine_paths: BoolProperty(
        name="Refine seam paths", default=True,
        description="Re-route open seam chains inside a corridor so they stop hopping between "
                    "quad rows; anchors, junctions and locked seams are kept")
    contour_mode: EnumProperty(
        name="Contours",
        items=[("auto", "Auto", "Straighten arcs that are edge loops and already nearly "
                "straight, smooth the rest"),
               ("panels", "Panels", "Straighten every arc between anchors"),
               ("organic", "Organic", "Smooth short oscillations, keep the overall curve"),
               ("off", "Off", "Leave the unwrap's contours as they are")],
        default="auto")
    budget: EnumProperty(
        name="Distortion budget",
        items=[("tight", "Tight", "P95 anisotropy 1.3, area 1.4"),
               ("normal", "Normal", "P95 anisotropy 1.6, area 1.8"),
               ("loose", "Loose", "P95 anisotropy 2.0, area 2.4")],
        default="normal")
    min_island: FloatProperty(
        name="Absorb islands below %", default=0.2, min=0.0, max=10.0,
        description="Islands smaller than this share of the surface are merged into a "
                    "neighbour under a relaxed budget")
    solver_iterations: IntProperty(
        name="Solver iterations", default=10, min=1, max=50,
        description="Minimum-stretch iterations of Blender's unwrap")
    min_island_faces: IntProperty(
        name="Minimum island faces", default=8, min=0,
        description="Absorb smaller islands into a validated neighbour, except whole "
                    "components and boundaries protected by the artist")
    target_enabled: BoolProperty(name='Target island count', default=False,
                                 description='Search up to five distortion budgets; an unreachable target is reported')
    target_min: IntProperty(name='Minimum', default=3, min=1)
    target_max: IntProperty(name='Maximum', default=8, min=1)
    deadline_s: FloatProperty(
        name="Time budget (seconds)", default=0.0, min=0.0,
        description="Return a validated partial layout when time runs out; 0 is unlimited. "
                    "An active solve finishes before the budget is checked again")
    hide_from: EnumProperty(
        name="Hide seams from",
        items=[('none', 'All views equally', 'Use overall surface visibility'),
               ('view', 'Current view', 'Prefer seams away from the current 3D view'),
               ('custom', 'Custom direction', 'Direction toward the viewer in object coordinates')]
              + [(axis, axis + ' object axis', 'Keep seams away from the ' + axis + ' side')
                 for axis in ('+X', '-X', '+Y', '-Y', '+Z', '-Z')], default='none')
    hide_direction: FloatVectorProperty(
        name="Viewing direction", size=3, default=(0, -1, 0), subtype='DIRECTION',
        description="Direction from the object toward the viewer, in object coordinates")
    # advanced
    show_advanced: BoolProperty(name="Advanced", default=False)
    w_turn: FloatProperty(name="Turn weight", default=4.0, min=0.0, max=20.0)
    w_flow: FloatProperty(name="Flow weight", default=2.0, min=0.0, max=20.0)
    w_drift: FloatProperty(name="Drift weight", default=3.0, min=0.0, max=20.0)
    corridor_rings: IntProperty(name="Corridor rings", default=2, min=1, max=6)
    smooth_radius: FloatProperty(name="Smooth radius", default=3.0, min=0.5, max=12.0,
                                 description="In local edge lengths")
    max_displacement: FloatProperty(name="Max displacement", default=1.5, min=0.1, max=6.0,
                                    description="Band for smoothed contours, in edge lengths")
    min_fusion: FloatProperty(name="Min fusion", default=0.15, min=0.0, max=0.9,
                              description="A merge must remove at least this share of the "
                                          "smaller island's outline")
    pack_after: EnumProperty(
        name="Then pack with",
        items=[("none", "Nothing", "Leave the islands where the unwrap put them"),
               ("blender", "Blender", "bpy.ops.uv.pack_islands")],
        default="none")
    keep_orientation: BoolProperty(name="Keep island orientation", default=True)
    pack_width: IntProperty(name="Texture width", default=2048, min=1)
    pack_height: IntProperty(name="Texture height", default=2048, min=1)
    pack_padding: FloatProperty(name="Island gap (texels)", default=4, min=0)
    pack_border: FloatProperty(name="Tile border (texels)", default=2, min=0)
    pack_tiles: StringProperty(name="UDIM tiles", default="1001",
                               description="Comma-separated UDIM tiles, for example 1001,1002")
    # status
    last_status: StringProperty(name="Last status", default="")
    last_report: StringProperty(name="Last report", default="")


def settings(context):
    return context.scene.cleanseams


CLASSES = (CleanSeamsSettings,)


def register():
    # T41: the enhanced controls come from the native registry. Without the engine they are not
    # declared (their values saved in a file stay in its ID properties) and nothing runs.
    NATIVE_FIELDS.clear()
    SECTIONS.clear()
    if native_module.available():
        declared, sections = native_declarations()
        CleanSeamsSettings.__annotations__.update(declared)
        NATIVE_FIELDS.extend(declared)
        SECTIONS.extend(sections)
    for c in CLASSES:
        bpy.utils.register_class(c)
    bpy.types.Scene.cleanseams = PointerProperty(type=CleanSeamsSettings)


def unregister():
    del bpy.types.Scene.cleanseams
    for c in reversed(CLASSES):
        bpy.utils.unregister_class(c)
    for name in NATIVE_FIELDS:
        CleanSeamsSettings.__annotations__.pop(name, None)
    NATIVE_FIELDS.clear()
    SECTIONS.clear()
