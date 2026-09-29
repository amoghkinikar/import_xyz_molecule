"""Headless smoke test inside real Blender.

Run from the repository root:

    blender -b --factory-startup --python-exit-code 1 -P tests/blender_smoke_test.py

On Windows, use the full path to blender.exe, e.g.
    "C:\\Program Files\\Blender Foundation\\Blender 5.1\\blender.exe" -b ...

Prints one PASS/FAIL line per check and exits with code 1 on any
failure. Paste the whole output back.
"""

import os
import struct
import sys
import tempfile
import traceback

import bpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import import_xyz_molecule as addon  # noqa: E402
from import_xyz_molecule import builder, nodes, scene  # noqa: E402

RESULTS = []


def check(name, condition, detail=""):
    RESULTS.append(bool(condition))
    print(f"{'PASS' if condition else 'FAIL'}  {name}"
          + (f"  [{detail}]" if detail else ""))


def write(tmp, name, text):
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def png_size(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    return struct.unpack(">II", head[16:24])


def import_xyz(path, **kw):
    before = set(bpy.data.objects)
    result = bpy.ops.import_mesh.xyz_molecule(filepath=path, **kw)
    new = [o for o in bpy.data.objects
           if o not in before and builder.is_molecule(o)]
    return result, (new[0] if new else None)


def evaluated_counts(obj):
    """(atom instances, bond polygons) of the evaluated molecule."""
    dg = bpy.context.evaluated_depsgraph_get()
    dg.update()
    instances = sum(1 for inst in dg.object_instances
                    if inst.is_instance and inst.parent
                    and inst.parent.original == obj)
    ev = obj.evaluated_get(dg)
    mesh = ev.to_mesh()
    polys = len(mesh.polygons) if mesh else 0
    ev.to_mesh_clear()
    return instances, polys


def main():
    print("Blender", bpy.app.version_string)
    addon.register()
    tmp = tempfile.mkdtemp(prefix="xyzmol_")

    # ---- 1. static molecule --------------------------------------------
    static = write(tmp, "co2.xyz",
                   "3\nCO2\nO -1.16 0 0\nC 0 0 0\nO 1.16 0 0\n")
    res, obj = import_xyz(static)
    check("static import finished", res == {'FINISHED'} and obj, res)
    mod = nodes.find_molecule_modifier(obj)
    check("static: Dynamic Bonds off",
          nodes.get_gn_input(mod, "Dynamic Bonds") in (False, 0))
    atoms, polys = evaluated_counts(obj)
    check("static: 3 atom instances", atoms == 3, atoms)
    check("static: bonds drawn", polys > 0, polys)
    check("static: BondCutoff edge attribute",
          obj.data.attributes.get("BondCutoff") is not None)

    # ---- 2. trajectory: C-C bond breaks at frame 3 ---------------------
    traj_text = ""
    for k, d in enumerate((1.40, 1.55, 3.00)):
        traj_text += f"2\nframe {k}\nC1 0 0 0\nC2 {d} 0 0\n"
    traj = write(tmp, "break.xyz", traj_text)
    bpy.context.scene.frame_start = 1
    bpy.context.scene.frame_end = 2
    res, obj = import_xyz(traj, frames='ALL')
    check("trajectory import finished", res == {'FINISHED'} and obj, res)
    key = obj.data.shape_keys
    check("trajectory: 3 shape keys",
          key is not None and len(key.key_blocks) == 3,
          None if key is None else len(key.key_blocks))
    check("trajectory: key frames increasing (absolute mode)",
          key is not None and not key.use_relative,
          None if key is None else [b.frame for b in key.key_blocks])
    check("trajectory: scene end frame extended to 3",
          bpy.context.scene.frame_end == 3, bpy.context.scene.frame_end)
    mod = nodes.find_molecule_modifier(obj)
    check("trajectory: Dynamic Bonds on",
          nodes.get_gn_input(mod, "Dynamic Bonds") in (True, 1))
    check("trajectory: labels mapped to carbon",
          builder.get_molecule_elements(obj) == ["C"],
          builder.get_molecule_elements(obj))
    for frame, dist, bonded in ((1, 1.40, True), (2, 1.55, True),
                                (3, 3.00, False)):
        bpy.context.scene.frame_set(frame)
        pos = scene.current_local_positions(bpy.context, obj)
        measured = abs(pos[1][0] - pos[0][0])
        check(f"trajectory frame {frame}: C-C distance {dist}",
              abs(measured - dist) < 1e-4, round(measured, 5))
        atoms, polys = evaluated_counts(obj)
        check(f"trajectory frame {frame}: bond "
              f"{'shown' if bonded else 'hidden'}",
              (polys > 0) == bonded, f"polys={polys}, atoms={atoms}")
    bpy.context.scene.frame_set(1)

    # ---- 3. periodic trajectory with changing cell ----------------------
    per_text = ""
    for k, a in enumerate((5.0, 5.2)):
        per_text += (f'2\nLattice="{a} 0 0 0 20 0 0 0 20" '
                     'Properties=species:S:1:pos:R:3 pbc="T F F"\n'
                     f"C 0.5 10 10\nC {a - 0.9} 10 10\n")
    per = write(tmp, "chain.extxyz", per_text)
    res, obj = import_xyz(per, frames='ALL')
    check("periodic trajectory import finished",
          res == {'FINISHED'} and obj, res)
    ghosts = obj.data.attributes.get("IsGhost")
    n_ghost = 0
    if ghosts is not None:
        vals = [0] * len(ghosts.data)
        ghosts.data.foreach_get("value", vals)
        n_ghost = sum(vals)
    check("periodic trajectory: ghost atoms present", n_ghost > 0, n_ghost)
    cell = next((c for c in obj.children if "_cell" in c.name), None)
    check("periodic trajectory: cell box animated",
          cell is not None and cell.data.shape_keys is not None)

    # ---- 4. Last frame only ------------------------------------------------
    res, obj_last = import_xyz(traj, frames='LAST')
    check("Frames=Last: no shape keys",
          obj_last is not None and obj_last.data.shape_keys is None)
    pos = scene.current_local_positions(bpy.context, obj_last)
    check("Frames=Last: last geometry",
          abs(abs(pos[1][0] - pos[0][0]) - 3.0) < 1e-4)

    # ---- 5. top-down figure render ---------------------------------------
    ribbon = "6\nribbon\n" + "".join(
        f"C {1.4 * i:.3f} {0.7 * (i % 2):.3f} 0\n" for i in range(6))
    rib = write(tmp, "ribbon.xyz", ribbon)
    res, obj = import_xyz(rib)
    sc = bpy.context.scene
    user_cam_data = bpy.data.cameras.new("UserCam")
    user_cam = bpy.data.objects.new("UserCam", user_cam_data)
    sc.collection.objects.link(user_cam)
    sc.camera = user_cam
    before = {
        "camera": sc.camera, "engine": sc.render.engine,
        "res": (sc.render.resolution_x, sc.render.resolution_y),
        "view": sc.view_settings.view_transform,
        "film": sc.render.film_transparent,
        "objects": len(bpy.data.objects),
        "cameras": len(bpy.data.cameras), "lights": len(bpy.data.lights),
    }
    out = os.path.join(tmp, "fig.png")
    try:
        path = scene.render_figure(bpy.context, obj, filepath=out,
                                   resolution=400)
        ok = os.path.isfile(path)
        check("figure rendered", ok, path)
        if ok:
            w, h = png_size(path)
            check("figure: wide image for a ribbon", w == 400 and h < 300,
                  f"{w}x{h}")
    except Exception as exc:
        traceback.print_exc()
        check("figure rendered", False, exc)
    after = {
        "camera": sc.camera, "engine": sc.render.engine,
        "res": (sc.render.resolution_x, sc.render.resolution_y),
        "view": sc.view_settings.view_transform,
        "film": sc.render.film_transparent,
        "objects": len(bpy.data.objects),
        "cameras": len(bpy.data.cameras), "lights": len(bpy.data.lights),
    }
    for field in before:
        check(f"figure: {field} restored", before[field] == after[field],
              f"{before[field]} -> {after[field]}")

    # ---- 6. unregister ---------------------------------------------------
    try:
        addon.unregister()
        check("unregister", True)
    except Exception as exc:
        traceback.print_exc()
        check("unregister", False, exc)

    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed} passed, {failed} failed")
    return failed


if __name__ == "__main__":
    try:
        failures = main()
    except Exception:
        traceback.print_exc()
        failures = 1
    sys.exit(1 if failures else 0)
