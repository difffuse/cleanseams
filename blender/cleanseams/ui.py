# SPDX-License-Identifier: GPL-3.0-or-later
"""Sidebar panel in the 3D viewport and the UV editor (drawing only). The enhanced controls
and the run messages are drawn as the native registry and run_diagnostics give them."""

import json

import bpy

from . import native_module, props

# Blender icons of the native message levels.
LEVEL_ICONS = {"info": "INFO", "warning": "ERROR", "error": "CANCEL"}


def draw_enhanced(layout, st):
    """The native sections of the enhanced controls: sections without a required policy (the
    seam policy) at the top, the others in a collapsible box, greyed out unless their policy
    is selected."""
    boxed = [section for section in props.SECTIONS if "requires" in section]
    for section in props.SECTIONS:
        if "requires" not in section:
            for name in section["fields"]:
                layout.prop(st, name)
    if not boxed:
        return
    box = layout.box()
    box.prop(st, "show_enhanced", icon="TRIA_DOWN" if st.show_enhanced else "TRIA_RIGHT", emboss=False)
    if not st.show_enhanced:
        return
    for section in boxed:
        col = box.column(align=True)
        col.active = st.policy == section["requires"]
        col.label(text=section["label"])
        for name in section["fields"]:
            col.prop(st, name)


def run_messages(st):
    """The native lines of the last run (level, topic, text), as run_messages returned them."""
    text = getattr(st, "last_messages", "")
    if not text:
        return []
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return []


def draw_settings(layout, st):
    layout.prop(st, "profile")
    draw_enhanced(layout, st)
    col = layout.column(align=True)
    col.prop(st, "generate")
    row = col.row()
    row.active = st.generate
    row.prop(st, "protect_existing")
    row = col.row()
    row.active = st.generate
    row.prop(st, "feature_angle")
    col.prop(st, "refine_paths")
    col.prop(st, "contour_mode")
    col.prop(st, "budget")
    col.prop(st, "min_island")
    col.prop(st, "min_island_faces")
    row = col.row()
    row.active = st.generate
    row.prop(st, 'target_enabled')
    if st.target_enabled:
        row = col.row(align=True)
        row.active = st.generate
        row.prop(st, 'target_min')
        row.prop(st, 'target_max')
    col.prop(st, "deadline_s")
    col.prop(st, "hide_from")
    if st.hide_from == 'custom':
        col.prop(st, "hide_direction")
    layout.separator()
    layout.prop(st, "pack_after")
    if st.pack_after != "none":
        layout.prop(st, "keep_orientation")
    box = layout.box()
    box.prop(st, "show_advanced", icon="TRIA_DOWN" if st.show_advanced else "TRIA_RIGHT",
             emboss=False)
    if st.show_advanced:
        c = box.column(align=True)
        c.prop(st, "w_turn")
        c.prop(st, "w_flow")
        c.prop(st, "w_drift")
        c.prop(st, "corridor_rings")
        c.prop(st, "smooth_radius")
        c.prop(st, "max_displacement")
        c.prop(st, "min_fusion")
        c.prop(st, "solver_iterations")


def draw_actions(layout, context, st):
    if not native_module.available():
        box = layout.box()
        box.alert = True
        for line in _wrap("Native engine unavailable: " + native_module.diagnostic(), 38):
            box.label(text=line)
    row = layout.row(align=True)
    row.scale_y = 1.4
    row.operator("cleanseams.run", icon="UV")
    row = layout.row(align=True)
    row.operator("cleanseams.restore", icon="LOOP_BACK")
    row.operator("cleanseams.pack", icon="PACKAGE")
    row.operator("cleanseams.report", icon="TEXT")
    if st.last_status:
        box = layout.box()
        for line in _wrap(st.last_status, 38):
            box.label(text=line)
        # T41: every native line of the run in full, with the icon of its level.
        for message in run_messages(st):
            icon = LEVEL_ICONS.get(message["level"], "INFO")
            for index, line in enumerate(_wrap(message["text"], 38, None)):
                box.label(text=line, icon=icon if index == 0 else "BLANK1")


def _wrap(text, width, limit=8):
    out, cur = [], ""
    for word in text.split(" "):
        if len(cur) + len(word) + 1 > width and cur:
            out.append(cur)
            cur = word
        else:
            cur = (cur + " " + word).strip()
    if cur:
        out.append(cur)
    return out[:limit]


class CLEANSEAMS_PT_view3d(bpy.types.Panel):
    bl_label = "Clean Seams"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Clean Seams"

    def draw(self, context):
        st = props.settings(context)
        draw_settings(self.layout, st)
        self.layout.separator()
        draw_actions(self.layout, context, st)


class CLEANSEAMS_PT_uv(bpy.types.Panel):
    bl_label = "Clean Seams"
    bl_space_type = "IMAGE_EDITOR"
    bl_region_type = "UI"
    bl_category = "Clean Seams"

    def draw(self, context):
        st = props.settings(context)
        draw_settings(self.layout, st)
        self.layout.separator()
        draw_actions(self.layout, context, st)


CLASSES = (CLEANSEAMS_PT_view3d, CLEANSEAMS_PT_uv)


def register():
    for c in CLASSES:
        bpy.utils.register_class(c)


def unregister():
    for c in reversed(CLASSES):
        bpy.utils.unregister_class(c)
