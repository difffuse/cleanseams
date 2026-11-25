# SPDX-License-Identifier: GPL-3.0-or-later
"""Operators: run the native UV stage on the active object, restore the previous state,
pack the islands with Blender, show the report. The engine runs in a native worker; this
module only transports arrays and answers its solve requests on the main thread."""

import json
import time

import bpy

from . import host, native_module, props

# restore snapshots per object name, valid only for the exact geometry they were taken on
SNAPSHOTS = {}
LAST = {}
# T41: the native run_diagnostics lines of the last run per object (JSON, drawn by the panel)
MESSAGES = {}
# Operator report types of the native message levels.
REPORT_TYPES = {"info": "INFO", "warning": "WARNING", "error": "ERROR"}


def resolve_options(st, context=None, obj=None):
    """The pipeline options for these settings (native resolution of budgets, profiles,
    percentages and the view direction)."""
    native = native_module.get()
    direction = matrix = None
    if st.hide_from == "view":
        from mathutils import Vector
        context = context or bpy.context
        obj = obj or context.active_object
        region = getattr(context, "region_data", None)
        if region is None and context.screen is not None:
            region = next((a.spaces.active.region_3d for a in context.screen.areas if a.type == "VIEW_3D"), None)
        if region is None or obj is None:
            raise ValueError("Current view requires a 3D viewport; choose an object axis or custom direction")
        direction = list(region.view_rotation @ Vector((0, 0, 1)))
        linear = obj.matrix_world.to_3x3()
        matrix = [linear[row][col] for row in range(3) for col in range(3)]
    return native.resolve_ui_options_json(props.settings_json(st), view_direction=direction, object_to_world=matrix)


# Runs in progress (their sessions are cancelled and closed on unregister or file load).
JOBS = set()


class Refusal(ValueError):
    """A precondition the user can fix (no faces, library-linked mesh): reported as a
    warning, nothing started and nothing written."""


class EngineFailure(RuntimeError):
    """The native engine ended the run with an error rather than a result (where the
    source add-on raised): reported as an error with the engine's message, nothing
    written."""


# Evidence for tests: stale results refused before any write.
STALE = {"refused": 0}


class RunJob:
    """One native run on one object, stepped on the main thread: start, step until the
    worker is done (answering its solve requests), then finish. Nothing is written to the
    object before finish(), and finish() writes only when the object, its mesh, its UV
    layers and the native geometry fingerprint are still the ones the run started from."""

    def __init__(self, context, obj, st, report_fn=None):
        self.native = native_module.get()
        self.obj, self.report_fn = obj, report_fn
        self.started = time.time()
        self.guard = host.ModeGuard(obj)
        self.session = self.service = None
        self.guard.__enter__()
        try:
            me = obj.data
            if len(me.polygons) == 0:
                raise Refusal("mesh has no faces")
            refusal = host.library_refusal(me)
            if refusal:
                raise Refusal(refusal)
            self.users = host.mesh_users(me)
            self.snap = host.Snapshot(obj)
            options = resolve_options(st, context, obj)
            # The painted attributes the options name are transferred raw (T29). The native
            # list is empty unless the policy is enhanced_v1 (GB M2): other runs transfer none,
            # and their stale-state digest covers what it covered at Gate A.
            self.data = host.HostData(obj, self.native.painted_attribute_names(options))
            self.service = host.SolveService(st.solver_iterations)
            # the frozen add-on had no memory limit: half the physical memory (Gate A F3)
            self.session = self.native.Session(backend="host", memory_budget=self.native.recommended_memory_budget())
            self.session.start(self.data.run_input(options, st.protect_existing))
        except BaseException:
            self.close()
            raise
        JOBS.add(self)

    def step(self):
        """Drain progress and answer at most one solve request; True once the run ended."""
        for event in self.session.poll_progress():
            if self.report_fn:
                self.report_fn(event)
        request = self.session.next_host_request()
        if request is not None:
            self.session.submit_host_response(self.service.answer(request))
            return False
        return self.session.ready()

    def cancel(self):
        self.session.request_cancel()

    def close(self):
        JOBS.discard(self)
        if self.session is not None:
            self.session.close()
        if self.service is not None:
            self.service.close()
        guard, self.guard = self.guard, None
        if guard is not None:
            guard.__exit__(None, None, None)

    def finish(self):
        """Apply an accepted result to the unchanged object. Returns (ok, message)."""
        try:
            try:
                # The native fingerprint of the live data, taken right before writing.
                current = self.data.live_fingerprint()
                result = self.session.result_for_apply(current)
            except (host.StaleState, self.native.ContractError) as exc:
                if isinstance(exc, self.native.ContractError) and not str(exc).startswith("apply.geometry_fingerprint"):
                    raise
                STALE["refused"] += 1
                return False, f"The object changed while CleanSeams ran; nothing applied ({exc})"
            obj = self.obj
            if result["diagnostics"] and not result["cancelled"]:
                # The operator reports "Clean Seams failed ...: <message>" as the source's did
                # for an exception, instead of the bare error code as a refusal reason.
                raise EngineFailure(result["diagnostics"][0][2])
            report = json.loads(result["report_json"])
            report["blender_time"] = time.time() - self.started
            if len(self.users) > 1:
                report["shared_mesh_users"] = self.users
            messages = self.native.run_messages(json.dumps(report), result["accepted"], result["reason"],
                                                result["invalid_charts"], self.data.nonuniform_scale)
            LAST[obj.name] = json.loads(messages["report_json"])
            MESSAGES[obj.name] = messages["messages_json"]
            if result["cancelled"]:
                return False, "cancelled; kept the previous state"
            if not result["accepted"]:
                return False, messages["status"]
            try:
                self.data.apply(result["seams"], result["uv"])
            except host.StaleState as exc:
                STALE["refused"] += 1
                return False, f"The object changed while CleanSeams ran; nothing applied ({exc})"
            except Exception:  # noqa: BLE001 - a partial write is rolled back
                self.snap.restore(obj)
                raise
            SNAPSHOTS[obj.name] = self.snap
            status = messages["status"]
            if len(self.users) > 1:
                # Shared mesh data: every user shows the new seams and UVs (as the source
                # add-on did); say so rather than changing other objects silently.
                others = [name for name in self.users if name != obj.name]
                status += f"; mesh data is shared, also changed on: {', '.join(others)}"
            return True, status
        finally:
            self.close()


def execute_on_object(context, obj, st, report_fn=None, poll_s=0.0005):
    """Runs the native pipeline on one object and applies the result. Returns (ok, message)."""
    job = RunJob(context, obj, st, report_fn)
    try:
        while not job.step():
            time.sleep(poll_s)
    except BaseException:
        job.close()
        raise
    return job.finish()


def cancel_all_jobs():
    """Cancel and close every run in progress (unregister, file load)."""
    for job in list(JOBS):
        try:
            job.cancel()
        finally:
            job.close()


class CLEANSEAMS_OT_run(bpy.types.Operator):
    bl_idname = "cleanseams.run"
    bl_label = "Clean Seams"
    bl_description = ("Generate or refine seams along the quad flow, unwrap, clean the UV "
                      "contours and apply only what passes validation")
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return native_module.available() and obj is not None and obj.type == "MESH" and len(obj.data.polygons) > 0

    def execute(self, context):
        st = props.settings(context)
        obj = context.active_object
        MESSAGES.pop(obj.name, None)
        st.last_messages = ""
        try:
            ok, msg = execute_on_object(context, obj, st,
                                        report_fn=lambda e: print(f"[cleanseams] {e['stage']} {e['done']}/{e['total']}", flush=True))
        except Refusal as exc:
            st.last_status = f"refused: {exc}"
            self.report({"WARNING"}, f"Clean Seams did not run: {exc}")
            return {"CANCELLED"}
        except Exception as exc:  # noqa: BLE001
            st.last_status = f"error: {exc}"
            self.report({"ERROR"}, f"Clean Seams failed and restored the object: {exc}")
            return {"CANCELLED"}
        st.last_status = msg
        st.last_report = json.dumps(LAST.get(obj.name, {}), allow_nan=False)
        st.last_messages = MESSAGES.get(obj.name, "")
        if not ok:
            self.report({"WARNING"}, msg)
            return {"CANCELLED"}
        self.report({"INFO"}, msg)
        if st.pack_after != "none":
            bpy.ops.cleanseams.pack()
        return {"FINISHED"}


class CLEANSEAMS_OT_restore(bpy.types.Operator):
    bl_idname = "cleanseams.restore"
    bl_label = "Restore previous"
    bl_description = "Put back the seams, UVs and pins the object had before the last run"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return native_module.available() and obj is not None and obj.name in SNAPSHOTS

    def execute(self, context):
        obj = context.active_object
        snap = SNAPSHOTS.get(obj.name)
        with host.ModeGuard(obj):
            ok = snap.restore(obj)
        if not ok:
            self.report({"WARNING"}, "The mesh changed since the run; nothing restored")
            return {"CANCELLED"}
        del SNAPSHOTS[obj.name]
        st = props.settings(context)
        st.last_status = "restored the previous seams and UVs"
        st.last_messages = ""  # the run's lines no longer describe the object
        return {"FINISHED"}


class CLEANSEAMS_OT_pack(bpy.types.Operator):
    bl_idname = "cleanseams.pack"
    bl_label = "Pack islands"
    bl_description = "Pack the validated islands with Blender"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == "MESH"

    def execute(self, context):
        st = props.settings(context)
        obj = context.active_object
        choice = st.pack_after
        if choice not in {"none", "blender"}:
            self.report({"WARNING"}, "Packing choice is no longer available; select Blender to pack")
            return {"CANCELLED"}
        with host.ModeGuard(obj):
            bpy.ops.object.mode_set(mode="EDIT")
            bpy.ops.mesh.select_all(action="SELECT")
            try:
                bpy.ops.uv.pack_islands(rotate=not st.keep_orientation, margin=0.003)
            finally:
                bpy.ops.object.mode_set(mode="OBJECT")
        if native_module.available():
            # T41: a host packer the engine does not verify is named as such (native line),
            # after the run's own lines; the run's status line is kept.
            packed = native_module.get().packing_messages("", "blender", st.last_messages)
            st.last_messages = packed["messages_json"]
            self.report({REPORT_TYPES[packed["level"]]}, packed["status"])
        return {"FINISHED"}


class CLEANSEAMS_OT_report(bpy.types.Operator):
    bl_idname = "cleanseams.report"
    bl_label = "Show report"
    bl_description = "Write the last report as JSON into a text block"

    def execute(self, context):
        st = props.settings(context)
        name = "cleanseams_report.json"
        text = bpy.data.texts.get(name) or bpy.data.texts.new(name)
        text.clear()
        try:
            text.write(json.dumps(json.loads(st.last_report or "{}"), indent=2))
        except json.JSONDecodeError:
            text.write(st.last_report)
        self.report({"INFO"}, f"Report written to text block {name}")
        return {"FINISHED"}


CLASSES = (CLEANSEAMS_OT_run, CLEANSEAMS_OT_restore, CLEANSEAMS_OT_pack, CLEANSEAMS_OT_report)


@bpy.app.handlers.persistent
def _before_file_load(*_):
    """A file load frees every object of the old file: stop runs on them, and drop restore
    snapshots, which describe objects of the old file only."""
    cancel_all_jobs()
    SNAPSHOTS.clear()


def register():
    for c in CLASSES:
        bpy.utils.register_class(c)
    if _before_file_load not in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.append(_before_file_load)


def unregister():
    if _before_file_load in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.remove(_before_file_load)
    cancel_all_jobs()
    SNAPSHOTS.clear()
    LAST.clear()
    MESSAGES.clear()
    for c in reversed(CLASSES):
        bpy.utils.unregister_class(c)
