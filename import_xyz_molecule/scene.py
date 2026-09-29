"""Scene helpers: viewport framing, default-cube cleanup, camera fitting,
and the one-click figure render.
"""

import json
import math
import os

import bpy

from . import builder, framing, nodes

FIGURE_CAMERA_NAME = "XYZ_Figure_Camera"
FIGURE_SUN_NAME = "XYZ_Figure_Sun"


# ---------------------------------------------------------------------------
# Viewport
# ---------------------------------------------------------------------------

def remove_default_cube():
    """Delete the untouched startup cube (and only that).

    Conservative checks: a mesh object named 'Cube' with exactly 8
    vertices, no modifiers, no children, no animation data.
    """
    obj = bpy.data.objects.get("Cube")
    if (obj is None or obj.type != 'MESH' or obj.data is None
            or len(obj.data.vertices) != 8 or obj.modifiers
            or obj.children or obj.animation_data):
        return False
    mesh = obj.data
    bpy.data.objects.remove(obj)
    if mesh.users == 0:
        bpy.data.meshes.remove(mesh)
    return True


def frame_object_in_viewports(context, obj):
    """Zoom every 3D viewport in the current screen onto ``obj``."""
    try:
        for other in context.selected_objects:
            if other != obj:
                other.select_set(False)
        obj.select_set(True)
        context.view_layer.objects.active = obj

        window = context.window
        if window is None:
            return False
        framed = False
        for area in window.screen.areas:
            if area.type != 'VIEW_3D':
                continue
            for region in area.regions:
                if region.type == 'WINDOW':
                    with context.temp_override(
                        window=window, area=area, region=region
                    ):
                        bpy.ops.view3d.view_selected()
                    framed = True
        return framed
    except Exception:
        return False  # headless / unusual context: framing is best-effort


# ---------------------------------------------------------------------------
# What to frame
# ---------------------------------------------------------------------------

def _read_point_attr(mesh, name, default, n, vector=False):
    attr = mesh.attributes.get(name)
    if attr is None or attr.domain != 'POINT' or len(attr.data) != n:
        return [default] * n
    if vector:
        flat = [0.0] * (3 * n)
        attr.data.foreach_get("vector", flat)
        return [tuple(flat[3 * i:3 * i + 3]) for i in range(n)]
    values = [0.0] * n if attr.data_type == 'FLOAT' else [0] * n
    attr.data.foreach_get("value", values)
    return values


def current_local_positions(context, obj):
    """Atom positions (object space) at the current frame, including
    trajectory shape keys, without evaluating the Geometry Nodes tree."""
    mesh = obj.data
    n = len(mesh.vertices)
    base = [0.0] * (3 * n)
    mesh.vertices.foreach_get("co", base)
    basis = [tuple(base[3 * i:3 * i + 3]) for i in range(n)]
    key = mesh.shape_keys
    if key is None or len(key.key_blocks) < 2:
        return basis

    def coords(block):
        flat = [0.0] * (3 * n)
        block.data.foreach_get("co", flat)
        return [tuple(flat[3 * i:3 * i + 3]) for i in range(n)]

    try:
        key_eval = key.evaluated_get(context.evaluated_depsgraph_get())
    except Exception:
        key_eval = key
    blocks = list(key.key_blocks)
    if not key.use_relative:
        return framing.interpolate_absolute(
            [b.frame for b in blocks], [coords(b) for b in blocks],
            key_eval.eval_time,
        )
    values = [b.value for b in key_eval.key_blocks][1:]
    return framing.mix_relative(basis, [coords(b) for b in blocks[1:]],
                                values)


def _slab_elements(obj, context):
    try:
        settings = json.loads(obj.get(builder.PROP_SETTINGS, "{}"))
    except (TypeError, ValueError):
        settings = {}
    slab = settings.get("slab_elements")
    if not slab:
        addon = context.preferences.addons.get(__package__)
        prefs = addon.preferences if addon else None
        raw = getattr(prefs, "default_slab_elements", "Au,Cu,Pt,Ga,Si")
        slab = [s.strip() for s in raw.split(",") if s.strip()]
    return set(slab)


def visible_atoms_world(context, obj, frame_target='WHOLE'):
    """World-space centers and drawn radii of the atoms that will show in
    the render: hidden elements are skipped, ghost atoms only count when
    shown, and ``frame_target='ADSORBATE'`` skips the slab elements
    (falling back to everything when nothing else is left).
    """
    mesh = obj.data
    n = len(mesh.vertices)
    positions = current_local_positions(context, obj)
    hidden = _read_point_attr(mesh, "ElementHidden", 0, n)
    ghost = _read_point_attr(mesh, "IsGhost", 0, n)
    index = _read_point_attr(mesh, "IndexOfElement", 0, n)
    elements = builder.get_molecule_elements(obj)

    display = getattr(obj, "xyz_mol_display", None)
    mod = nodes.find_molecule_modifier(obj)

    def setting(prop, socket, default):
        if display is not None:
            return getattr(display, prop)
        value = nodes.get_gn_input(mod, socket) if mod else None
        return default if value is None else value

    atom_scale = setting("atom_scale", "Atom Scale", 0.7)
    bond_radius = setting("bond_radius", "Bond Radius", 0.1)
    jitter = setting("scale_jitter", "Scale Jitter", 0.0)
    show_ghosts = setting("show_ghosts", "Show Ghost Atoms", True)
    style = obj.get("xyz_style", 'BALL_AND_STICK')
    attr_name, factor = framing.drawn_radius_factor(
        style, atom_scale, bond_radius, jitter
    )
    base_radii = (_read_point_attr(mesh, attr_name, 0.0, n)
                  if attr_name else [1.0] * n)

    keep = [i for i in range(n)
            if not hidden[i] and (show_ghosts or not ghost[i])]
    if frame_target == 'ADSORBATE':
        slab = _slab_elements(obj, context)
        adsorbate = [
            i for i in keep
            if not (0 <= index[i] < len(elements)
                    and elements[index[i]] in slab)
        ]
        if adsorbate:
            keep = adsorbate
    if not keep:
        keep = list(range(n))

    mw = obj.matrix_world
    scale = max(abs(v) for v in mw.to_scale())
    points = [tuple(mw @ _vec(positions[i])) for i in keep]
    radii = [base_radii[i] * factor * scale for i in keep]
    return points, radii


def _vec(p):
    from mathutils import Vector
    return Vector(p)


# ---------------------------------------------------------------------------
# Figure camera and light (temporary: removed after the render)
# ---------------------------------------------------------------------------

def make_topdown_camera(context, frame):
    cam_data = bpy.data.cameras.new(FIGURE_CAMERA_NAME)
    cam_data.type = 'ORTHO'
    cam_data.ortho_scale = frame["ortho_scale"]
    cam_data.sensor_fit = 'AUTO'
    cam_data.clip_start = frame["clip_start"]
    cam_data.clip_end = frame["clip_end"]
    cam_obj = bpy.data.objects.new(FIGURE_CAMERA_NAME, cam_data)
    cx, cy = frame["center"]
    cam_obj.location = (cx, cy, frame["camera_z"])
    cam_obj.rotation_euler = (0.0, 0.0, 0.0)   # looks down -Z, +Y is up
    context.scene.collection.objects.link(cam_obj)
    return cam_obj


def make_sun(context):
    sun_data = bpy.data.lights.new(FIGURE_SUN_NAME, 'SUN')
    sun_data.energy = 3.0
    sun_data.angle = math.radians(15.0)  # soft shadows
    sun_obj = bpy.data.objects.new(FIGURE_SUN_NAME, sun_data)
    sun_obj.rotation_euler = (math.radians(35), 0.0, math.radians(35))
    context.scene.collection.objects.link(sun_obj)
    return sun_obj


def _remove_object_and_data(obj):
    data = obj.data
    bpy.data.objects.remove(obj)
    if data is not None and data.users == 0:
        if isinstance(data, bpy.types.Camera):
            bpy.data.cameras.remove(data)
        elif isinstance(data, bpy.types.Light):
            bpy.data.lights.remove(data)


# ---------------------------------------------------------------------------
# Figure rendering
# ---------------------------------------------------------------------------

def _pick_engine(scene, high_quality):
    if high_quality:
        candidates = ['CYCLES']
    else:
        candidates = ['BLENDER_EEVEE_NEXT', 'BLENDER_EEVEE', 'CYCLES']
    for engine in candidates:
        try:
            scene.render.engine = engine
            return engine
        except TypeError:
            continue
    return scene.render.engine


_EEVEE_AO_SETTINGS = (("use_gtao", True),        # EEVEE legacy AO
                      ("use_raytracing", True),  # EEVEE-Next raytraced AO
                      ("gtao_distance", 1.0))


def _snapshot(scene):
    """Everything render_figure changes, so it can be put back."""
    render = scene.render
    view = scene.view_settings
    saved = {
        "camera": scene.camera,
        "engine": render.engine,
        "resolution_x": render.resolution_x,
        "resolution_y": render.resolution_y,
        "resolution_percentage": render.resolution_percentage,
        "film_transparent": render.film_transparent,
        "filepath": render.filepath,
        "file_format": render.image_settings.file_format,
        "color_mode": render.image_settings.color_mode,
        "view_transform": view.view_transform,
        "look": view.look,
        "eevee": {},
        "cycles_samples": None,
    }
    eevee = getattr(scene, "eevee", None)
    for attr, _value in _EEVEE_AO_SETTINGS:
        if eevee is not None and hasattr(eevee, attr):
            saved["eevee"][attr] = getattr(eevee, attr)
    cycles = getattr(scene, "cycles", None)
    if cycles is not None and hasattr(cycles, "samples"):
        saved["cycles_samples"] = cycles.samples
    return saved


def _restore(scene, saved):
    render = scene.render
    view = scene.view_settings
    for step in (
        lambda: setattr(scene, "camera", saved["camera"]),
        lambda: setattr(render, "engine", saved["engine"]),
        lambda: setattr(render, "resolution_x", saved["resolution_x"]),
        lambda: setattr(render, "resolution_y", saved["resolution_y"]),
        lambda: setattr(render, "resolution_percentage",
                        saved["resolution_percentage"]),
        lambda: setattr(render, "film_transparent",
                        saved["film_transparent"]),
        lambda: setattr(render, "filepath", saved["filepath"]),
        lambda: setattr(render.image_settings, "file_format",
                        saved["file_format"]),
        lambda: setattr(render.image_settings, "color_mode",
                        saved["color_mode"]),
        lambda: setattr(view, "view_transform", saved["view_transform"]),
        lambda: setattr(view, "look", saved["look"]),
    ):
        try:
            step()
        except Exception:
            pass   # one failed restore must not block the others
    eevee = getattr(scene, "eevee", None)
    for attr, value in saved["eevee"].items():
        try:
            setattr(eevee, attr, value)
        except Exception:
            pass
    if saved["cycles_samples"] is not None:
        try:
            scene.cycles.samples = saved["cycles_samples"]
        except Exception:
            pass


def _composite_over_white(filepath):
    """Flatten the saved transparent PNG onto a white background.

    Runs as a post-process on the written file, so it needs no
    compositor and works identically on Blender 4.x and 5.x.
    """
    import numpy as np  # bundled with Blender

    img = bpy.data.images.load(filepath)
    try:
        n = img.size[0] * img.size[1] * 4
        pixels = np.empty(n, dtype=np.float32)
        img.pixels.foreach_get(pixels)
        pixels = pixels.reshape(-1, 4)
        alpha = pixels[:, 3:4]
        pixels[:, 0:3] = pixels[:, 0:3] * alpha + (1.0 - alpha)
        pixels[:, 3] = 1.0
        img.pixels.foreach_set(pixels.reshape(-1))
        img.filepath_raw = filepath
        img.file_format = 'PNG'
        img.save()
    finally:
        bpy.data.images.remove(img)


def figure_output_path(obj, suffix="_figure"):
    source = obj.get(builder.PROP_SOURCE, "") if obj else ""
    if source:
        base, _ext = os.path.splitext(source)
        return f"{base}{suffix}.png"
    directory = bpy.path.abspath("//") or os.path.expanduser("~")
    name = obj.name if obj else "molecule"
    return os.path.join(directory, f"{name}{suffix}.png")


def render_figure(context, obj, filepath=None, resolution=1600,
                  background='TRANSPARENT', high_quality=False,
                  samples=64, frame_target='WHOLE'):
    """One-click figure: top-down orthographic view of the molecule,
    framed on its visible atoms, AO, transparent or white background,
    PNG saved next to the source file.

    ``resolution`` is the longer image side; the shorter side follows the
    structure's footprint. ``frame_target``: 'WHOLE' or 'ADSORBATE'.

    The camera and sun exist only for the render, and every scene/render
    setting touched here is restored afterwards.

    Returns the absolute output path.
    """
    scene = context.scene
    render = scene.render

    points, radii = visible_atoms_world(context, obj, frame_target)
    frame = framing.topdown_frame(points, radii, long_side=resolution)

    if filepath is None:
        filepath = figure_output_path(obj)

    saved = _snapshot(scene)
    camera = sun = None
    try:
        camera = make_topdown_camera(context, frame)
        sun = make_sun(context)
        scene.camera = camera

        engine = _pick_engine(scene, high_quality)
        if engine == 'CYCLES':
            try:
                scene.cycles.samples = samples
            except AttributeError:
                pass
        else:
            eevee = getattr(scene, "eevee", None)
            for attr, value in _EEVEE_AO_SETTINGS:
                if eevee is not None and hasattr(eevee, attr):
                    try:
                        setattr(eevee, attr, value)
                    except Exception:
                        pass

        # Element colors are specified in sRGB; the default AgX / Filmic
        # transforms desaturate them. Standard shows them as defined.
        try:
            scene.view_settings.view_transform = 'Standard'
            scene.view_settings.look = 'None'
        except TypeError:
            pass

        res_x, res_y = frame["resolution"]
        render.resolution_x = res_x
        render.resolution_y = res_y
        render.resolution_percentage = 100
        render.film_transparent = True
        render.image_settings.file_format = 'PNG'
        render.image_settings.color_mode = 'RGBA'
        render.filepath = filepath

        bpy.ops.render.render(write_still=True)
    finally:
        _restore(scene, saved)
        for temp in (camera, sun):
            if temp is not None:
                try:
                    _remove_object_and_data(temp)
                except Exception:
                    pass

    if background == 'WHITE':
        _composite_over_white(filepath)
    return filepath
