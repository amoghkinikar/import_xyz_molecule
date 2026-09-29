"""
XYZ Molecule Importer
=====================

A Blender add-on that imports .xyz / extended-xyz molecular geometry files
(and generic point CSVs) as meshes with atoms instanced via Geometry
Nodes, colored per element, and bonds built as skinned edges.

Inspired by the design of Molecular Nodes
(https://github.com/BradyAJohnston/MolecularNodes).

Features
--------
* File > Import > XYZ Molecule (.xyz) and Points CSV (.csv).
* Auto-center at the origin, auto-frame the viewport, and remove the
  untouched default cube, so the molecule is visible immediately.
* One-click "Render Figure": top-down orthographic view framed on the
  visible atoms (or only the adsorbate), ambient occlusion, transparent
  or white background, PNG saved next to the source file.
* Multi-frame .xyz (trajectories, geometry optimisations) imported as an
  animation, with bonds that form and break as the atoms move.
* Display presets: Ball-and-Stick, Space-Filling (van der Waals radii),
  Licorice/Sticks, Wireframe.
* Extended XYZ support: Lattice (unit cell box + periodic bonds across
  cell boundaries via ghost atoms) and extra columns (forces, charges,
  spins) imported as point attributes.
* Build Supercell: repeat a periodic structure NxMxL along its lattice
  vectors, with bonds throughout the interior and the cell redrawn as a
  grid of the original cells.
* Color by any scalar attribute with viridis/plasma/coolwarm colormaps.
* Element tools: hide hydrogens, per-element visibility, select by
  element.
* Reload button: re-import the same file with all current settings.
* Instance swapping: instance any scene object for a chosen element
  (or declare it up front with an "object" key in the element JSON).
* Per-instance variation (random rotation, scale jitter, seed) and
  alignment of instances to a per-point vector attribute.

Installation
------------
Blender 5.1+: Edit > Preferences > Get Extensions > (dropdown, top
right) Install from Disk... > select the .zip.
"""

bl_info = {
    "name": "XYZ Molecule Importer",
    "author": "Amogh Kinikar",
    "version": (2, 2, 0),
    "blender": (5, 1, 0),
    "location": "File > Import > XYZ Molecule (.xyz)  |  "
                "3D View > Sidebar > XYZ Molecule",
    "description": "Import .xyz/extxyz structures and point CSVs as "
                   "Geometry Nodes molecules",
    "doc_url": "https://github.com/amoghkinikar/import_xyz_molecule",
    "category": "Import-Export",
}

# Support add-on reloading: re-import submodules on repeated registration.
if "bpy" in locals():
    import importlib
    for _mod_name in ("elements", "parsing", "bonds", "trajectory",
                      "framing", "materials", "nodes", "builder", "scene",
                      "properties", "operators", "ui"):
        if _mod_name in locals():
            importlib.reload(locals()[_mod_name])

try:
    import bpy
except ImportError:
    # Outside Blender (e.g. unit tests importing the pure-Python
    # submodules); register()/unregister() are unavailable there.
    bpy = None

if bpy is not None:
    from . import (  # noqa: E402,F401
        bonds, builder, elements, framing, materials, nodes, operators,
        parsing, properties, scene, trajectory, ui,
    )

    def register():
        properties.register()
        for cls in operators.classes:
            bpy.utils.register_class(cls)
        for cls in ui.classes:
            bpy.utils.register_class(cls)
        bpy.types.TOPBAR_MT_file_import.append(ui.menu_func_import)

    def unregister():
        bpy.types.TOPBAR_MT_file_import.remove(ui.menu_func_import)
        for cls in reversed(ui.classes):
            bpy.utils.unregister_class(cls)
        for cls in reversed(operators.classes):
            bpy.utils.unregister_class(cls)
        properties.unregister()


if __name__ == "__main__":
    register()
