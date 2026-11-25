# SPDX-License-Identifier: GPL-3.0-or-later
"""Host glue: bulk reads and writes of Blender arrays on the main thread, restore
snapshots, and the service that answers the native engine's solve requests.

No domain decision is taken here. Geometry mapping, fingerprints, options, chart
selection, validation and summaries are native (cleanseams_native); this module only
moves arrays between Blender and the engine and runs Blender's own unwrap operator on
the isolated chart mesh the engine asks for.
"""

import threading

import bmesh
import bpy
import numpy as np

from . import compat, native_module

# Evidence for tests: every host service call and the thread it ran on.
HOST_CALLS = {"total": 0, "off_main_thread": 0}


class StaleState(RuntimeError):
    """The object, its mesh or its UV layers changed while the engine worked; nothing of
    the result may be written."""


def alive(struct):
    """Whether a Blender object/mesh reference still points at live data (deletion, undo
    and file loads free the underlying data and leave the Python wrapper dangling)."""
    try:
        struct.name
        return True
    except ReferenceError:
        return False


def uv_state(me):
    """UV layer names in order and the active layer: what a run's result is written into."""
    active = me.uv_layers.active.name if me.uv_layers.active else None
    return tuple(layer.name for layer in me.uv_layers), active


def state_digest(me, painted_names=()):
    """Digest of what a run reads besides geometry and replaces: seams, locks, anchors, the
    painted constraint attributes it transfers, the active UV layer's values and its pins.
    A change while the engine works makes the result stale (it was computed from, and would
    overwrite, other data). The arrays are gathered here; the digest is native. The painted
    names are the native painted_attribute_names of the run's options: none unless its
    policy is enhanced_v1, so other runs digest exactly the Gate A parts."""
    parts = []
    for name, domain, count in (("uv_seam", "EDGE", len(me.edges)), ("cleanseams_lock", "EDGE", len(me.edges)),
                                ("cleanseams_anchor", "POINT", len(me.vertices))):
        parts.append((name, _bool_attribute(me, name, domain, count)))
    for name in painted_names:
        painted = _painted_attribute(me, name)
        parts.append((f"painted:{name}:{painted['domain']}:{painted['data_type']}" if painted else f"painted:{name}",
                      painted["values"] if painted else None))
    layer = me.uv_layers.active
    if layer is not None:
        uv = np.empty(len(me.loops) * 2, dtype=np.float32)
        layer.uv.foreach_get("vector", uv)
        parts.append((layer.name, uv))
        if len(layer.pin) == len(me.loops):
            pins = np.zeros(len(me.loops), dtype=bool)
            layer.pin.foreach_get("value", pins)
            parts.append(("pin", pins))
    return native_module.get().state_digest(parts)


def mesh_users(me):
    """Objects of the current file that use this mesh datablock."""
    return [o.name for o in bpy.data.objects if o.type == "MESH" and o.data == me]  # bpy.data keeps name order


def library_refusal(me):
    """A reason when the mesh cannot be written in this file (linked library data or a
    non-editable override), else None. Writing it would be lost or rejected by Blender."""
    if me.library is not None:
        return f"mesh '{me.name}' is linked from library '{me.library.filepath}'; make it local first"
    if not me.is_editable:
        return f"mesh '{me.name}' is not editable (library override); make it local or editable first"
    return None


def _record_host_call():
    HOST_CALLS["total"] += 1
    if threading.get_ident() != threading.main_thread().ident:
        HOST_CALLS["off_main_thread"] += 1


# ------------------------------------------------------------------- reading
def raw_arrays(me):
    """Blender's float32/int32 buffers in raw host order (what the fingerprint hashes)."""
    nv, ne, nl, nf = len(me.vertices), len(me.edges), len(me.loops), len(me.polygons)
    co = np.empty(nv * 3, dtype=np.float32)
    me.vertices.foreach_get("co", co)
    edges = np.empty(ne * 2, dtype=np.int32)
    me.edges.foreach_get("vertices", edges)
    loop_vertices = np.empty(nl, dtype=np.int32)
    me.loops.foreach_get("vertex_index", loop_vertices)
    loop_edges = np.empty(nl, dtype=np.int32)
    me.loops.foreach_get("edge_index", loop_edges)
    starts = np.empty(nf, dtype=np.int32)
    me.polygons.foreach_get("loop_start", starts)
    totals = np.empty(nf, dtype=np.int32)
    me.polygons.foreach_get("loop_total", totals)
    return co.reshape(-1, 3), edges.reshape(-1, 2), loop_vertices, loop_edges, starts, totals


def fingerprint(me):
    """The exact geometry and index-order fingerprint (native)."""
    return native_module.get().host_geometry_fingerprint(*raw_arrays(me))


def _bool_attribute(me, name, domain, count):
    attr = me.attributes.get(name)
    if attr is None or attr.domain != domain or attr.data_type != "BOOLEAN":
        return None
    values = np.zeros(count, dtype=bool)
    attr.data.foreach_get("value", values)
    return values


# Blender attribute types with one scalar "value" per element, read as float64.
SCALAR_TYPES = {"BOOLEAN": bool, "FLOAT": np.float32, "INT": np.int32, "INT8": np.int8}


def _painted_attribute(me, name):
    """A painted attribute exactly as Blender stores it: its domain, data type and, for a
    scalar type, its values in host order. Domains, types, ranges, clamping and conflicts are
    interpreted natively (T29); nothing is converted, filtered or guessed here."""
    attr = me.attributes.get(name)
    if attr is None:
        return None
    values = np.zeros(0, dtype=np.float64)
    dtype = SCALAR_TYPES.get(attr.data_type)
    if dtype is not None:
        raw = np.zeros(len(attr.data), dtype=dtype)
        attr.data.foreach_get("value", raw)
        values = np.ascontiguousarray(raw, dtype=np.float64)
    return {"name": name, "domain": attr.domain, "data_type": attr.data_type, "values": values}


class HostData:
    """One object's mesh as the native engine maps it, plus the host index maps."""

    def __init__(self, obj, painted_names=()):
        native = native_module.get()
        me = obj.data
        self.obj = obj
        # The identity the result belongs to, rechecked immediately before writing.
        self.object_pointer, self.mesh_pointer = obj.as_pointer(), me.as_pointer()
        self.uv_state = uv_state(me)
        # The attribute names come from the native options (documented names or aliases).
        self.painted_names = tuple(painted_names)
        self.state = state_digest(me, self.painted_names)
        raw = raw_arrays(me)
        nv, ne, nl = len(me.vertices), len(me.edges), len(me.loops)
        attributes = {}
        seams = _bool_attribute(me, "uv_seam", "EDGE", ne)
        if seams is not None:
            attributes["seams"] = seams
        locked = _bool_attribute(me, "cleanseams_lock", "EDGE", ne)
        if locked is not None:
            attributes["locked"] = locked
        anchors = _bool_attribute(me, "cleanseams_anchor", "POINT", nv)
        if anchors is not None:
            attributes["anchors"] = anchors
        painted = [p for p in (_painted_attribute(me, name) for name in self.painted_names) if p is not None]
        if painted:
            attributes["painted"] = painted
        self.uv_layer = me.uv_layers.active.name if me.uv_layers.active else None
        if self.uv_layer:
            layer = me.uv_layers[self.uv_layer]
            uv = np.empty(nl * 2, dtype=np.float64)
            layer.uv.foreach_get("vector", uv)
            attributes["uv"] = uv.reshape(-1, 2)
            if len(layer.pin) == nl:
                pins = np.zeros(nl, dtype=bool)
                layer.pin.foreach_get("value", pins)
                attributes["pins"] = pins
        matrix = np.ascontiguousarray([[float(v) for v in row] for row in obj.matrix_world], dtype=np.float64)
        scale = np.ascontiguousarray([float(v) for v in obj.scale], dtype=np.float64)
        self.mapped = native.read_host_mesh(*raw, matrix, scale, memory_budget=native.recommended_memory_budget(),
                                            attributes=attributes)
        self.n_edges, self.n_loops = ne, nl

    @property
    def nonuniform_scale(self):
        return bool(self.mapped["nonuniform_scale"])

    def run_input(self, options_json, protect_existing):
        m = self.mapped
        seams = m["seams"].astype(bool)
        locked = m["locked"].astype(bool)
        if protect_existing:
            locked = locked | seams
        data = {"positions": m["positions"], "face_offsets": m["face_offsets"],
                "corner_vertices": m["corner_vertices"], "seam_mask": seams, "lock_mask": locked,
                "pin_corners": m["pin_corners"], "pin_uv": m["pin_uv"],
                "anchors": np.asarray(m["anchors"], dtype=np.int64), "options": options_json,
                "geometry_fingerprint": m["geometry_fingerprint"]}
        if m["uv"] is not None:
            data["uv"] = m["uv"]
        if m["painted"]:
            data["painted_attributes"] = m["painted"]
        return data

    def live_fingerprint(self):
        """The native fingerprint of the object's current data, after checking that the
        object, its mesh datablock and its UV layers are still the ones the run read.
        Raises StaleState otherwise; nothing is written by this call."""
        obj = self.obj
        if not alive(obj):
            raise StaleState("the object was deleted, or replaced by undo or a file load, while CleanSeams ran")
        if obj.as_pointer() != self.object_pointer or obj.type != "MESH" or not alive(obj.data):
            raise StaleState("the object changed while CleanSeams ran")
        me = obj.data
        if me.as_pointer() != self.mesh_pointer:
            raise StaleState("the object's mesh data was replaced while CleanSeams ran")
        if uv_state(me) != self.uv_state:
            raise StaleState("the UV layers changed while CleanSeams ran")
        if state_digest(me, self.painted_names) != self.state:
            raise StaleState("the seams, UVs or pins changed while CleanSeams ran")
        refusal = library_refusal(me)
        if refusal:
            raise StaleState(refusal)
        return fingerprint(me)

    def apply(self, seams_core, uv_core):
        """Write seams and UVs to the object; the native side checks what Blender stored.
        Every array is staged before the first write; the caller has just rechecked the
        native fingerprint and pins (Session.result_for_apply, or the native packing apply)
        and the live identity is checked again."""
        m = self.mapped
        seams_bl = np.zeros(self.n_edges, dtype=bool)
        to_host = m["core_edge_to_host"]
        valid = to_host >= 0
        seams_bl[to_host[valid]] = np.asarray(seams_core, dtype=bool)[valid]
        uv_bl = pins_bl = None
        if uv_core is not None:
            loops = m["host_loop_indices"]
            uv_bl = np.zeros((self.n_loops, 2), dtype=np.float64)
            uv_bl[loops] = uv_core
            pins_bl = np.zeros(self.n_loops, dtype=bool)
            pins_bl[loops] = np.asarray(m["pins"], dtype=bool)
        if self.live_fingerprint() != m["geometry_fingerprint"]:
            raise StaleState("the geometry changed while CleanSeams ran")
        me = self.obj.data
        write_seams(me, seams_bl)
        if uv_bl is not None:
            if not self.uv_layer:
                me.uv_layers.new(name="UVMap")
                self.uv_layer = me.uv_layers.active.name
            write_uv(me, self.uv_layer, uv_bl)
            write_pins(self.obj, self.uv_layer, pins_bl)
            # Blender stores UVs in float32: the native check accepts exactly the float32
            # rounding of the validated map and nothing else (raises; the caller rolls back).
            back = np.empty(self.n_loops * 2, dtype=np.float32)
            me.uv_layers[self.uv_layer].uv.foreach_get("vector", back)
            native_module.get().check_stored_uv(np.ascontiguousarray(uv_core, dtype=np.float64),
                                                np.asarray(loops, dtype=np.int64), back.reshape(-1, 2))
        me.update()


# ------------------------------------------------------------------- writing
def write_seams(me, seams_bl):
    attr = compat.seam_attribute(me, create=True)
    attr.data.foreach_set("value", np.asarray(seams_bl, dtype=bool))


def write_uv(me, layer_name, uv_bl):
    layer = me.uv_layers[layer_name]
    layer.uv.foreach_set("vector", np.asarray(uv_bl, dtype=np.float64).ravel())


def write_pins(obj, layer_name, pins_bl):
    """Pins through bmesh (object mode): creates the pin attribute when missing."""
    me = obj.data
    pins_bl = np.asarray(pins_bl, dtype=bool)
    layer = me.uv_layers[layer_name]
    if len(layer.pin) == len(me.loops):
        layer.pin.foreach_set("value", pins_bl)
        return
    if not pins_bl.any():
        return
    bm = bmesh.new()
    bm.from_mesh(me)
    uvl = bm.loops.layers.uv.get(layer_name)
    bm.faces.ensure_lookup_table()
    for f in bm.faces:
        for loop in f.loops:
            loop[uvl].pin_uv = bool(pins_bl[loop.index])
    bm.to_mesh(me)
    bm.free()


# ----------------------------------------------------------------- snapshots
class Snapshot:
    """Seams, UVs and pins before a run, keyed by the exact native geometry fingerprint."""

    def __init__(self, obj):
        me = obj.data
        self.obj_name = obj.name
        self.fingerprint = fingerprint(me)
        ne, nl = len(me.edges), len(me.loops)
        self.seams = np.zeros(ne, dtype=bool)
        attr = compat.seam_attribute(me)
        if attr is not None:
            attr.data.foreach_get("value", self.seams)
        self.uv_layer = me.uv_layers.active.name if me.uv_layers.active else None
        self.uv = None
        self.pins = None
        if self.uv_layer:
            layer = me.uv_layers[self.uv_layer]
            uv = np.zeros(nl * 2, dtype=np.float64)
            layer.uv.foreach_get("vector", uv)
            self.uv = uv.reshape(-1, 2)
            self.pins = np.zeros(nl, dtype=bool)
            if len(layer.pin) == nl:
                layer.pin.foreach_get("value", self.pins)

    def restore(self, obj):
        """Put the state back only onto the exact geometry it was taken from."""
        if not alive(obj) or obj.type != "MESH" or library_refusal(obj.data):
            return False
        me = obj.data
        if fingerprint(me) != self.fingerprint:
            return False
        write_seams(me, self.seams)
        if self.uv_layer and self.uv is not None and self.uv_layer in me.uv_layers:
            write_uv(me, self.uv_layer, self.uv)
            write_pins(obj, self.uv_layer, self.pins)
        elif self.uv_layer is None:
            for layer in list(me.uv_layers):
                me.uv_layers.remove(layer)
        me.update()
        return True


class ModeGuard:
    """Restore object mode, active object and selection when done."""

    def __init__(self, obj):
        self.obj = obj
        self.mode = obj.mode
        active = bpy.context.view_layer.objects.active
        # Pointers, not wrappers: objects may be deleted, or replaced by undo or a file
        # load, before the guard exits.
        self.active = active.as_pointer() if active is not None else None
        self.selected = {o.as_pointer() for o in bpy.context.view_layer.objects if o is not None and o.select_get()}

    def __enter__(self):
        bpy.context.view_layer.objects.active = self.obj
        if self.obj.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        return self

    def __exit__(self, *exc):
        if alive(self.obj) and self.obj.name in bpy.context.view_layer.objects:
            bpy.context.view_layer.objects.active = self.obj
            if self.obj.mode != "OBJECT":
                bpy.ops.object.mode_set(mode="OBJECT")
            if self.mode != "OBJECT":
                try:
                    bpy.ops.object.mode_set(mode=self.mode)
                except RuntimeError:
                    pass
        active = None
        for o in bpy.context.view_layer.objects:
            if o is not None:
                pointer = o.as_pointer()
                o.select_set(pointer in self.selected)
                if pointer == self.active:
                    active = o
        if active is not None:
            bpy.context.view_layer.objects.active = active
        return False


# ------------------------------------------------------------ solve service
class SolveService:
    """Answers the engine's solve requests with Blender's unwrap on an isolated scratch
    mesh (the source BlenderBackend._native_unwrap), always on the main thread."""

    def __init__(self, iterations):
        self.iterations = int(iterations)
        self._scratch = None

    def close(self):
        scratch, self._scratch = self._scratch, None
        # A file load or undo may already have freed the scratch object with its file.
        if scratch is not None and alive(scratch):
            mesh = scratch.data
            bpy.data.objects.remove(scratch, do_unlink=True)
            if mesh is not None and alive(mesh) and mesh.users == 0:
                bpy.data.meshes.remove(mesh)
            bpy.context.view_layer.update()

    def answer(self, request):
        _record_host_call()
        identity = request["identity"]
        try:
            uv, method = self._unwrap(request)
        except Exception as exc:  # noqa: BLE001 - reported to the engine as a typed failure
            return {"identity": identity, "kind": request["kind"], "success": False, "error": str(exc)}
        return {"identity": identity, "kind": request["kind"], "success": True, "uv": uv, "method": method}

    def _unwrap(self, request):
        positions = request["positions"]
        offsets = request["face_offsets"].astype(np.int64)
        corners = request["corner_vertices"]
        n_corners = len(corners)
        method = compat.unwrap_method(request["method"])
        if self._scratch is not None and not alive(self._scratch):
            self._scratch = None  # freed by undo or a file load: build a new one
        if self._scratch is None:
            me = bpy.data.meshes.new("CleanSeams solver mesh")
            self._scratch = bpy.data.objects.new("CleanSeams solver", me)
            self._scratch.hide_render = True
            bpy.context.scene.collection.objects.link(self._scratch)
            bpy.context.view_layer.update()
        me = self._scratch.data
        me.clear_geometry()
        me.vertices.add(len(positions))
        me.vertices.foreach_set("co", np.ascontiguousarray(positions, dtype=np.float64).ravel())
        me.loops.add(n_corners)
        me.loops.foreach_set("vertex_index", np.asarray(corners, dtype=np.int32))
        me.polygons.add(len(offsets) - 1)
        me.polygons.foreach_set("loop_start", offsets[:-1].astype(np.int32))
        me.polygons.foreach_set("loop_total", np.diff(offsets).astype(np.int32))
        me.update(calc_edges=True)
        if not me.uv_layers:
            me.uv_layers.new(name="UVMap")
        initial = request["initial_uv"] if request["initial_uv"] is not None else np.zeros((n_corners, 2))
        pinned = {int(c): request["pin_uv"][i] for i, c in enumerate(request["pin_corners"])}
        selection_mode = tuple(bpy.context.tool_settings.mesh_select_mode)
        with ModeGuard(self._scratch):
            for obj in bpy.context.view_layer.objects:
                obj.select_set(obj == self._scratch)
            bpy.ops.object.mode_set(mode="EDIT")
            try:
                bpy.context.tool_settings.mesh_select_mode = (False, False, True)
                bpy.ops.mesh.select_all(action="SELECT")
                bm = bmesh.from_edit_mesh(me)
                bm.faces.ensure_lookup_table()
                layer = bm.loops.layers.uv.verify()
                for face, start in zip(bm.faces, offsets[:-1]):
                    for offset, loop in enumerate(face.loops):
                        corner = int(start) + offset
                        loop[layer].uv = pinned.get(corner, initial[corner])
                        loop[layer].pin_uv = corner in pinned
                bmesh.update_edit_mesh(me)
                kw = compat.unwrap_kwargs(method)
                if method == "MINIMUM_STRETCH":
                    kw["iterations"] = int(request["iterations"])
                bpy.ops.uv.unwrap(**kw)
            finally:
                bpy.ops.object.mode_set(mode="OBJECT")
                bpy.context.tool_settings.mesh_select_mode = selection_mode
        values = np.zeros(n_corners * 2)
        me.uv_layers.active.uv.foreach_get("vector", values)
        return values.reshape(-1, 2), method
