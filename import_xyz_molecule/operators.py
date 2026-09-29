"""Operators: import (.xyz / .csv), reload, display presets, element
tools, instance swapping, color-by-attribute, figure rendering, and the
edit-mode attribute stamp.
"""

import json
import os

import bmesh
import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    StringProperty,
)
from bpy.types import Operator
from bpy_extras.io_utils import ImportHelper

from . import builder, materials, nodes, parsing, scene
from .elements import (
    load_user_elements,
    merge_elements,
    normalize_symbol,
    pseudo_element_color,
)
from .properties import (
    active_molecule,
    element_enum_items,
    seed_display_settings,
    sync_display_to_modifier,
)


def get_prefs(context):
    addon = context.preferences.addons.get(__package__)
    return addon.preferences if addon else None


# ---------------------------------------------------------------------------
# Shared import machinery
# ---------------------------------------------------------------------------

class ImportSettingsMixin:
    """Properties and post-import behavior shared by the XYZ and CSV
    importers."""

    user_elements_path: StringProperty(
        name="Element Overrides (JSON)",
        description="Optional JSON overriding element radii/colors, and "
                    'optionally template objects: {"Au": {"radius": 150, '
                    '"color": "FFD700", "object": "MyMesh"}}',
        subtype='FILE_PATH',
        default="",
    )
    atom_scale: FloatProperty(
        name="Atom Scale",
        description="Sphere radius as a fraction of the element radius",
        default=0.7, min=0.0,
    )
    bond_scale: FloatProperty(
        name="Bond Cutoff Scale",
        description="Two atoms are bonded when their distance is below the "
                    "sum of their covalent radii times this factor",
        default=1.1, min=0.5, max=2.0,
    )
    bond_radius: FloatProperty(
        name="Bond Radius",
        description="Radius of the bond cylinders (Angstrom)",
        default=0.1, min=0.0,
    )
    world_scale: FloatProperty(
        name="World Scale",
        description="Blender units per Angstrom (object scale)",
        default=0.01, min=0.0001, precision=4,
    )
    center: BoolProperty(
        name="Center at Origin",
        description="Move the molecule's bounding-box center to the world "
                    "origin (the original offset is stored on the object)",
        default=True,
    )
    frame_view: BoolProperty(
        name="Frame in Viewport",
        description="Zoom the 3D viewport onto the imported molecule",
        default=True,
    )
    remove_cube: BoolProperty(
        name="Remove Default Cube",
        description="Delete the untouched startup cube so the molecule "
                    "is not hidden next to it",
        default=True,
    )
    setup_world: BoolProperty(
        name="Set Up Studio World",
        description="Create and assign the blackbody/highlight world "
                    "lighting",
        default=True,
    )

    def base_settings(self):
        return {
            "user_elements_path": self.user_elements_path,
            "atom_scale": self.atom_scale,
            "bond_scale": self.bond_scale,
            "bond_radius": self.bond_radius,
            "world_scale": self.world_scale,
            "center": self.center,
        }

    def draw_common_geometry(self, layout):
        box = layout.box()
        box.label(text="Geometry", icon='MESH_UVSPHERE')
        box.prop(self, "atom_scale")
        box.prop(self, "bond_scale")
        box.prop(self, "bond_radius")
        box.prop(self, "world_scale")

    def draw_common_scene(self, layout):
        box = layout.box()
        box.label(text="Scene", icon='SCENE_DATA')
        box.prop(self, "center")
        box.prop(self, "frame_view")
        box.prop(self, "remove_cube")
        box.prop(self, "setup_world")

    def finish_import(self, context, obj, info, warnings=()):
        seed_display_settings(
            obj,
            atom_scale=self.atom_scale,
            bond_radius=self.bond_radius,
            show_atoms=True,
            show_bonds=True,
            show_ghosts=True,
        )
        if self.remove_cube:
            scene.remove_default_cube()
        if self.setup_world:
            materials.make_world(context)
        if self.frame_view:
            scene.frame_object_in_viewports(context, obj)
        for warning in warnings:
            self.report({'WARNING'}, warning)
        if info.get("frame_end_extended"):
            self.report({'INFO'}, "Scene end frame extended to "
                                  f"{info['frame_end_extended']}")

        parts = [f"{info['n_real']} atoms", f"{info['n_bonds']} bonds"]
        if info.get("n_frames", 1) > 1:
            parts.append(f"{info['n_frames']} frames")
        if info["n_ghosts"]:
            parts.append(f"{info['n_ghosts']} periodic ghost atoms")
        if info["extra_attributes"]:
            parts.append(
                "attributes: " + ", ".join(info["extra_attributes"])
            )
        self.report(
            {'INFO'},
            f"Imported {', '.join(parts)} from "
            f"{os.path.basename(self.filepath)}"
        )


# ---------------------------------------------------------------------------
# XYZ import
# ---------------------------------------------------------------------------

class IMPORT_MESH_OT_xyz_molecule(Operator, ImportSettingsMixin,
                                  ImportHelper):
    """Import an .xyz or extended-xyz file as a Geometry Nodes molecule"""
    bl_idname = "import_mesh.xyz_molecule"
    bl_label = "Import XYZ Molecule"
    bl_options = {'REGISTER', 'UNDO'}

    filename_ext = ".xyz"
    filter_glob: StringProperty(
        default="*.xyz;*.extxyz", options={'HIDDEN'}
    )

    periodic_bonds: BoolProperty(
        name="Periodic Bonds",
        description="When the file has a Lattice (extended xyz), add bonds "
                    "across cell boundaries via ghost atoms",
        default=True,
    )
    show_unit_cell: BoolProperty(
        name="Unit Cell Box",
        description="Draw the lattice as a wireframe box (extended xyz)",
        default=True,
    )
    scale_slab: BoolProperty(
        name="Scale Slab Atoms",
        description="Draw slab elements with larger spheres. This is a "
                    "display setting only: bond detection always uses the "
                    "true covalent radii",
        default=False,
    )
    slab_elements: StringProperty(
        name="Slab Elements",
        description="Comma-separated list of slab element symbols",
        default="Au,Cu,Pt,Ga,Si",
    )
    slab_scale: FloatProperty(
        name="Slab Scale",
        description="Scale factor applied to slab atom radii",
        default=2.0, min=0.0,
    )
    frames: EnumProperty(
        name="Frames",
        description="What to import from a multi-frame file (trajectory, "
                    "geometry optimisation)",
        items=[
            ('ALL', "All (animated)",
             "One shape key per frame, one frame per scene frame; bonds "
             "form and break as the atoms move"),
            ('FIRST', "First", "Only the first frame"),
            ('LAST', "Last",
             "Only the last frame (e.g. the relaxed geometry of a "
             "geometry optimisation)"),
        ],
        default='ALL',
    )
    frame_step: IntProperty(
        name="Frame Step",
        description="Import every Nth frame (long MD runs)",
        default=1, min=1,
    )
    unwrap: BoolProperty(
        name="Unwrap Periodic Jumps",
        description="For animated periodic files, undo atoms jumping "
                    "across the cell between frames (minimum image), so "
                    "they do not fly across the cell in the animation",
        default=True,
    )

    def invoke(self, context, event):
        prefs = get_prefs(context)
        if prefs is not None:
            if not self.user_elements_path:
                self.user_elements_path = prefs.user_elements_path
            self.slab_elements = prefs.default_slab_elements
        return super().invoke(context, event)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        box = layout.box()
        box.label(text="Element Data", icon='PRESET')
        box.prop(self, "user_elements_path")

        self.draw_common_geometry(layout)

        box = layout.box()
        box.label(text="Periodic (extended xyz)", icon='MESH_GRID')
        box.prop(self, "periodic_bonds")
        box.prop(self, "show_unit_cell")

        box = layout.box()
        box.label(text="Multi-frame files", icon='RENDER_ANIMATION')
        box.prop(self, "frames")
        sub = box.column()
        sub.enabled = self.frames == 'ALL'
        sub.prop(self, "frame_step")
        sub.prop(self, "unwrap")

        box = layout.box()
        box.label(text="Slab", icon='MOD_LATTICE')
        box.prop(self, "scale_slab")
        sub = box.column()
        sub.enabled = self.scale_slab
        sub.prop(self, "slab_elements")
        sub.prop(self, "slab_scale")

        self.draw_common_scene(layout)

    def execute(self, context):
        overrides, err = load_user_elements(self.user_elements_path)
        if err:
            self.report({'ERROR'}, err)
            return {'CANCELLED'}
        elements = merge_elements(overrides)

        try:
            frames = parsing.parse_xyz_frames(
                self.filepath, known_symbols=elements,
                which=self.frames, step=self.frame_step,
            )
        except (OSError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        structure = frames[0]
        if not structure["symbols"]:
            self.report({'ERROR'}, "No atoms found in file.")
            return {'CANCELLED'}

        settings = self.base_settings()
        settings.update({
            "kind": 'XYZ',
            "source": self.filepath,
            "periodic_bonds": self.periodic_bonds,
            "show_unit_cell": self.show_unit_cell,
            "scale_slab": self.scale_slab,
            "slab_elements": [
                normalize_symbol(s)
                for s in self.slab_elements.split(",") if s.strip()
            ],
            "slab_scale": self.slab_scale,
            "frames": self.frames,
            "frame_step": self.frame_step,
            "unwrap": self.unwrap,
        })

        name = os.path.splitext(os.path.basename(self.filepath))[0]
        try:
            obj, info = builder.build_molecule(
                context, structure, elements, name, settings, frames=frames,
            )
        except Exception as exc:  # keep Blender responsive on failure
            self.report({'ERROR'}, f"Failed to build molecule: {exc}")
            return {'CANCELLED'}

        if overrides:
            self.report(
                {'INFO'},
                "Element overrides applied for: "
                + ", ".join(sorted(overrides)),
            )
        warnings = list(structure.get("warnings", ()))
        n_file = structure.get("n_frames_in_file", 1)
        if n_file > 1 and self.frames != 'ALL':
            warnings.append(
                f"File has {n_file} frames; imported only the "
                f"{self.frames.lower()} one"
            )
        self.finish_import(context, obj, info, warnings=warnings)
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# CSV import
# ---------------------------------------------------------------------------

def resolve_csv_labels(structure, elements):
    """Map CSV labels onto element entries; invent pseudo-elements (with
    distinct colors) for labels that are not chemical symbols."""
    label_map = {}
    unknown = []
    for label in sorted(set(structure["symbols"])):
        norm = normalize_symbol(label)
        if norm in elements:
            label_map[label] = norm
        else:
            unknown.append(label)
    for i, label in enumerate(unknown):
        elements[label] = {
            "number": 0,
            "radius": 60,        # 0.6 Angstrom instances by default
            "vdw_radius": 120,
            "color": pseudo_element_color(i),
            "object": "",
        }
        label_map[label] = label
    structure["symbols"] = [label_map[s] for s in structure["symbols"]]
    return unknown


class IMPORT_MESH_OT_csv_points(Operator, ImportSettingsMixin,
                                ImportHelper):
    """Import a label,x,y,z[,extra...] CSV as instanced points.

    Labels matching chemical elements get chemistry treatment; other
    labels get auto-assigned colors. Extra columns become attributes"""
    bl_idname = "import_mesh.xyz_csv_points"
    bl_label = "Import Points CSV"
    bl_options = {'REGISTER', 'UNDO'}

    filename_ext = ".csv"
    filter_glob: StringProperty(
        default="*.csv;*.tsv;*.txt", options={'HIDDEN'}
    )

    detect_bonds: BoolProperty(
        name="Detect Bonds",
        description="Connect nearby points with bonds/edges "
                    "(distance-based, like the molecular importer)",
        default=False,
    )

    def invoke(self, context, event):
        # CSVs usually carry generic data: sensible different defaults.
        self.world_scale = 1.0
        self.atom_scale = 1.0
        prefs = get_prefs(context)
        if prefs is not None and not self.user_elements_path:
            self.user_elements_path = prefs.user_elements_path
        return super().invoke(context, event)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        box = layout.box()
        box.label(text="Element Data", icon='PRESET')
        box.prop(self, "user_elements_path")
        self.draw_common_geometry(layout)
        layout.prop(self, "detect_bonds")
        self.draw_common_scene(layout)

    def execute(self, context):
        overrides, err = load_user_elements(self.user_elements_path)
        if err:
            self.report({'ERROR'}, err)
            return {'CANCELLED'}
        elements = merge_elements(overrides)

        try:
            structure = parsing.parse_csv(self.filepath)
        except (OSError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        if not structure["symbols"]:
            self.report({'ERROR'}, "No points found in file.")
            return {'CANCELLED'}

        unknown = resolve_csv_labels(structure, elements)

        settings = self.base_settings()
        settings.update({
            "kind": 'CSV',
            "source": self.filepath,
            "detect_bonds": self.detect_bonds,
        })

        name = os.path.splitext(os.path.basename(self.filepath))[0]
        try:
            obj, info = builder.build_molecule(
                context, structure, elements, name, settings
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Failed to build point cloud: {exc}")
            return {'CANCELLED'}

        if unknown:
            self.report(
                {'INFO'},
                "Non-element labels got auto colors: "
                + ", ".join(unknown),
            )
        self.finish_import(context, obj, info)
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Rebuilding from the source file (shared by Reload and Build Supercell)
# ---------------------------------------------------------------------------

def has_lattice(obj):
    """True when the molecule was imported from a file with a Lattice."""
    try:
        return bool(json.loads(obj.get(builder.PROP_LATTICE, '""')))
    except (TypeError, ValueError):
        return False


def stored_repeats(obj):
    try:
        return json.loads(obj.get(builder.PROP_REPEATS, "[1, 1, 1]"))
    except (TypeError, ValueError):
        return (1, 1, 1)


def normalize_repeats_for_ui(obj):
    """Current supercell repeats of ``obj`` as a plain int tuple."""
    return builder.normalize_repeats(stored_repeats(obj))


def primitive_atom_count(obj):
    """Atoms in one primitive cell (total real atoms / current repeats)."""
    counts = sum(json.loads(obj.get(builder.PROP_COUNTS, "{}")).values())
    nx, ny, nz = builder.normalize_repeats(stored_repeats(obj))
    cells = nx * ny * nz
    return counts // cells if cells and counts else counts


def rebuild_from_source(context, obj, setting_overrides=None):
    """Re-parse the object's source file and rebuild it in place.

    ``setting_overrides`` is merged into the stored import settings, so
    callers can change one thing (supercell repeats, say) while keeping
    everything else the user chose at import time.

    Returns ``(obj, info)``; raises OSError/ValueError on bad input.
    """
    source = obj.get(builder.PROP_SOURCE, "")
    if not source or not os.path.isfile(source):
        raise ValueError(f"Source file not found: {source}")

    try:
        settings = json.loads(obj.get(builder.PROP_SETTINGS, "{}"))
    except ValueError:
        settings = {}
    settings.setdefault("kind", 'XYZ')
    settings["source"] = source
    settings["supercell"] = list(
        builder.normalize_repeats(stored_repeats(obj))
    )
    if setting_overrides:
        settings.update(setting_overrides)

    overrides, err = load_user_elements(settings.get("user_elements_path",
                                                     ""))
    if err:
        raise ValueError(err)
    elements = merge_elements(overrides)

    if settings["kind"] == 'CSV':
        structure = parsing.parse_csv(source)
        resolve_csv_labels(structure, elements)
        frames = [structure]
    else:
        # Objects imported before multi-frame support read frame 1 only.
        frames = parsing.parse_xyz_frames(
            source, known_symbols=elements,
            which=settings.get("frames", 'FIRST'),
            step=settings.get("frame_step", 1),
        )
        structure = frames[0]
    if not structure["symbols"]:
        raise ValueError("No atoms found in file.")

    mol_name = obj.get(builder.PROP_NAME) or obj.name
    obj, info = builder.build_molecule(
        context, structure, elements, mol_name, settings, existing_obj=obj,
        frames=frames,
    )
    sync_display_to_modifier(obj)
    return obj, info


# ---------------------------------------------------------------------------
# Reload
# ---------------------------------------------------------------------------

class XYZMOL_OT_reload(Operator):
    """Re-import the molecule's source file with the same settings,
    replacing the geometry in place (styles, visibility, and instance
    swaps are kept)"""
    bl_idname = "xyz_molecule.reload"
    bl_label = "Reload From File"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = active_molecule(context)
        return (obj is not None and context.mode == 'OBJECT'
                and obj.get(builder.PROP_SOURCE))

    def execute(self, context):
        obj = active_molecule(context)
        source = obj.get(builder.PROP_SOURCE, "")
        try:
            obj, info = rebuild_from_source(context, obj)
        except (OSError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        except Exception as exc:
            self.report({'ERROR'}, f"Failed to rebuild molecule: {exc}")
            return {'CANCELLED'}

        self.report(
            {'INFO'},
            f"Reloaded {info['n_real']} atoms, {info['n_bonds']} bonds "
            f"from {os.path.basename(source)}"
        )
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Supercell
# ---------------------------------------------------------------------------

class XYZMOL_OT_build_supercell(Operator):
    """Repeat this periodic structure along its lattice vectors, replacing
    the current object. Bonds form throughout the interior and the cell is
    redrawn as a grid of the original cells"""
    bl_idname = "xyz_molecule.build_supercell"
    bl_label = "Build Supercell"
    bl_options = {'REGISTER', 'UNDO'}

    nx: IntProperty(
        name="A", description="Repeats along the first lattice vector",
        default=2, min=1, max=64,
    )
    ny: IntProperty(
        name="B", description="Repeats along the second lattice vector",
        default=2, min=1, max=64,
    )
    nz: IntProperty(
        name="C", description="Repeats along the third lattice vector",
        default=1, min=1, max=64,
    )

    @classmethod
    def poll(cls, context):
        obj = active_molecule(context)
        if obj is None or context.mode != 'OBJECT':
            return False
        return bool(obj.get(builder.PROP_SOURCE)) and has_lattice(obj)

    def invoke(self, context, event):
        obj = active_molecule(context)
        self.nx, self.ny, self.nz = builder.normalize_repeats(
            stored_repeats(obj)
        )
        return context.window_manager.invoke_props_dialog(self)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        col = layout.column(align=True)
        col.prop(self, "nx")
        col.prop(self, "ny")
        col.prop(self, "nz")

        obj = active_molecule(context)
        per_cell = primitive_atom_count(obj)
        if per_cell:
            total = per_cell * self.nx * self.ny * self.nz
            box = layout.box()
            box.label(
                text=f"{per_cell} atoms per cell "
                     f"-> {total:,} atoms",
                icon='INFO',
            )
            if total > 50000:
                box.label(text="Large scenes may be slow to edit",
                          icon='ERROR')

    def execute(self, context):
        obj = active_molecule(context)
        repeats = (self.nx, self.ny, self.nz)
        try:
            obj, info = rebuild_from_source(
                context, obj, {"supercell": list(repeats)}
            )
        except (OSError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        except Exception as exc:
            self.report({'ERROR'}, f"Failed to build supercell: {exc}")
            return {'CANCELLED'}

        nx, ny, nz = repeats
        self.report(
            {'INFO'},
            f"{nx}x{ny}x{nz} supercell: {info['n_real']} atoms, "
            f"{info['n_bonds']} bonds"
        )
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Display presets
# ---------------------------------------------------------------------------

class XYZMOL_OT_set_style(Operator):
    """Switch the molecule's display style"""
    bl_idname = "xyz_molecule.set_style"
    bl_label = "Set Display Style"
    bl_options = {'REGISTER', 'UNDO'}

    style: EnumProperty(name="Style", items=nodes.STYLE_ITEMS)

    @classmethod
    def poll(cls, context):
        return nodes.find_molecule_modifier(active_molecule(context)) \
            is not None

    def execute(self, context):
        obj = active_molecule(context)
        mod = nodes.find_molecule_modifier(obj)
        nodes.set_gn_input(mod, "Style", nodes.STYLE_TO_INT[self.style])
        display = getattr(obj, "xyz_mol_display", None)
        if self.style == 'BALL_AND_STICK':
            if display is not None:
                display.atom_scale = 0.7      # pushes to the modifier
            else:
                nodes.set_gn_input(mod, "Atom Scale", 0.7)
        elif self.style == 'SPACE_FILLING':
            if display is not None:
                display.atom_scale = 1.0
            else:
                nodes.set_gn_input(mod, "Atom Scale", 1.0)
        obj["xyz_style"] = self.style
        obj.update_tag()
        context.view_layer.update()
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Element tools
# ---------------------------------------------------------------------------

class XYZMOL_OT_toggle_element(Operator):
    """Show or hide all atoms of one element"""
    bl_idname = "xyz_molecule.toggle_element"
    bl_label = "Toggle Element Visibility"
    bl_options = {'REGISTER', 'UNDO'}

    element: StringProperty(name="Element")

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        hidden = builder.get_hidden_elements(obj)
        if self.element in hidden:
            hidden.discard(self.element)
        else:
            hidden.add(self.element)
        builder.set_hidden_elements(obj, hidden)
        builder.apply_hidden_elements(obj)
        obj.update_tag()
        context.view_layer.update()
        return {'FINISHED'}


class XYZMOL_OT_hide_hydrogens(Operator):
    """Hide (or show, when already hidden) all hydrogen atoms"""
    bl_idname = "xyz_molecule.hide_hydrogens"
    bl_label = "Hide Hydrogens"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        if "H" not in builder.get_molecule_elements(obj):
            self.report({'INFO'}, "This molecule has no hydrogens.")
            return {'CANCELLED'}
        return bpy.ops.xyz_molecule.toggle_element(element="H")


class XYZMOL_OT_show_all_elements(Operator):
    """Make every element visible again"""
    bl_idname = "xyz_molecule.show_all_elements"
    bl_label = "Show All"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        builder.set_hidden_elements(obj, set())
        builder.apply_hidden_elements(obj)
        obj.update_tag()
        context.view_layer.update()
        return {'FINISHED'}


class XYZMOL_OT_select_element(Operator):
    """Select all atoms of one element (works in Object and Edit Mode)"""
    bl_idname = "xyz_molecule.select_element"
    bl_label = "Select Element"
    bl_options = {'REGISTER', 'UNDO'}

    element: EnumProperty(name="Element", items=element_enum_items)
    extend: BoolProperty(
        name="Extend",
        description="Add to the current selection instead of replacing it",
        default=False,
    )

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        elements_in_file = builder.get_molecule_elements(obj)
        if self.element not in elements_in_file:
            self.report({'ERROR'}, f"No element {self.element} here.")
            return {'CANCELLED'}
        target = elements_in_file.index(self.element)

        if context.mode == 'EDIT_MESH' and context.edit_object == obj:
            bm = bmesh.from_edit_mesh(obj.data)
            layer = bm.verts.layers.int.get("IndexOfElement")
            if layer is None:
                self.report({'ERROR'}, "IndexOfElement attribute missing.")
                return {'CANCELLED'}
            count = 0
            for vert in bm.verts:
                match = vert[layer] == target
                if match:
                    count += 1
                if self.extend:
                    vert.select_set(vert.select or match)
                else:
                    vert.select_set(match)
            bm.select_flush_mode()
            bmesh.update_edit_mesh(obj.data)
        else:
            mesh = obj.data
            index_attr = mesh.attributes.get("IndexOfElement")
            if index_attr is None:
                self.report({'ERROR'}, "IndexOfElement attribute missing.")
                return {'CANCELLED'}
            n = len(index_attr.data)
            indices = [0] * n
            index_attr.data.foreach_get("value", indices)
            count = 0
            for i, vert in enumerate(mesh.vertices):
                match = indices[i] == target
                if match:
                    count += 1
                vert.select = (vert.select or match) if self.extend \
                    else match
            mesh.update()

        self.report(
            {'INFO'},
            f"Selected {count} {self.element} atoms"
            + ("" if context.mode == 'EDIT_MESH'
               else " (enter Edit Mode to see/edit the selection)"),
        )
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Instance swapping
# ---------------------------------------------------------------------------

class XYZMOL_OT_swap_instance(Operator):
    """Replace the instanced object for one element with any object in
    the scene (the original stays where it is; the molecule uses a copy)"""
    bl_idname = "xyz_molecule.swap_instance"
    bl_label = "Swap Instance"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        tools = context.scene.xyz_mol_tools
        return (active_molecule(context) is not None
                and tools.swap_object is not None)

    def execute(self, context):
        obj = active_molecule(context)
        tools = context.scene.xyz_mol_tools
        target = tools.swap_object
        element = tools.swap_element

        elements_in_file = builder.get_molecule_elements(obj)
        if element not in elements_in_file:
            self.report({'ERROR'}, f"No element {element} in molecule.")
            return {'CANCELLED'}
        index = elements_in_file.index(element)

        collection = builder.find_atom_collection(obj)
        if collection is None:
            self.report({'ERROR'}, "Template collection not found.")
            return {'CANCELLED'}
        if target.name in collection.objects:
            self.report({'ERROR'}, "That object is already a template.")
            return {'CANCELLED'}

        mol_name = obj.get(builder.PROP_NAME) or obj.name
        template_name = builder.template_object_name(
            mol_name, index, element
        )

        old = builder.find_template_object(collection, index)
        old_data = None
        if old is not None:
            old_data = old.data
            bpy.data.objects.remove(old)
        if (old_data is not None and old_data.users == 0
                and isinstance(old_data, bpy.types.Mesh)
                and old_data.name.startswith(f"{mol_name}_atom_")):
            bpy.data.meshes.remove(old_data)

        from mathutils import Matrix
        new_obj = target.copy()      # linked copy: shares data with source
        new_obj.parent = None
        new_obj.matrix_world = Matrix.Identity(4)
        new_obj.name = template_name
        new_obj.hide_render = True
        collection.objects.link(new_obj)

        obj.update_tag()
        context.view_layer.update()
        self.report(
            {'INFO'},
            f"{element} atoms now instance '{target.name}'."
        )
        return {'FINISHED'}


class XYZMOL_OT_reset_instance(Operator):
    """Restore the default sphere instance for one element"""
    bl_idname = "xyz_molecule.reset_instance"
    bl_label = "Reset to Sphere"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        tools = context.scene.xyz_mol_tools
        element = tools.swap_element

        elements_in_file = builder.get_molecule_elements(obj)
        if element not in elements_in_file:
            self.report({'ERROR'}, f"No element {element} in molecule.")
            return {'CANCELLED'}
        index = elements_in_file.index(element)

        collection = builder.find_atom_collection(obj)
        if collection is None:
            self.report({'ERROR'}, "Template collection not found.")
            return {'CANCELLED'}

        mol_name = obj.get(builder.PROP_NAME) or obj.name
        old = builder.find_template_object(collection, index)
        if old is not None:
            bpy.data.objects.remove(old)

        atom_mat = materials.create_atom_material(mol_name)
        sphere = builder.make_sphere_object(
            builder.template_object_name(mol_name, index, element),
            atom_mat, (index * 3.0, 0.0, 0.0),
        )
        collection.objects.link(sphere)
        sphere.hide_render = True

        obj.update_tag()
        context.view_layer.update()
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Color by attribute
# ---------------------------------------------------------------------------

class XYZMOL_OT_color_by_attribute(Operator):
    """Color atoms (and bonds) by a scalar attribute using a scientific
    colormap"""
    bl_idname = "xyz_molecule.color_by_attribute"
    bl_label = "Color by Attribute"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        tools = context.scene.xyz_mol_tools
        if tools.color_attribute in ("", 'NONE'):
            self.report({'ERROR'}, "No scalar attribute to color by.")
            return {'CANCELLED'}
        vmin = None if tools.color_auto_range else tools.color_min
        vmax = None if tools.color_auto_range else tools.color_max
        try:
            vmin, vmax = materials.apply_colormap_to_mesh(
                obj.data, tools.color_attribute,
                colormap=tools.colormap, vmin=vmin, vmax=vmax,
            )
        except ValueError as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        tools.color_min, tools.color_max = float(vmin), float(vmax)
        obj.update_tag()
        context.view_layer.update()
        self.report(
            {'INFO'},
            f"Colored by '{tools.color_attribute}' "
            f"[{vmin:.4g} .. {vmax:.4g}]",
        )
        return {'FINISHED'}


class XYZMOL_OT_color_by_element(Operator):
    """Restore the standard per-element (CPK) colors"""
    bl_idname = "xyz_molecule.color_by_element"
    bl_label = "Color by Element"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        elements_in_file = builder.get_molecule_elements(obj)
        try:
            table = json.loads(obj[builder.PROP_TABLE])
        except (KeyError, TypeError, ValueError):
            table = merge_elements({})

        mesh = obj.data
        index_attr = mesh.attributes.get("IndexOfElement")
        if index_attr is None:
            self.report({'ERROR'}, "IndexOfElement attribute missing.")
            return {'CANCELLED'}
        n = len(index_attr.data)
        indices = [0] * n
        index_attr.data.foreach_get("value", indices)
        symbols = [
            elements_in_file[i] if 0 <= i < len(elements_in_file) else "X"
            for i in indices
        ]
        fallback = merge_elements({})
        table = {**fallback, **table}
        materials.restore_element_colors(mesh, symbols, table)
        obj.update_tag()
        context.view_layer.update()
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Figure rendering
# ---------------------------------------------------------------------------

class XYZMOL_OT_render_figure(Operator):
    """Render a top-down orthographic figure of the molecule, framed on
    its visible atoms: ambient occlusion, transparent or white
    background, PNG saved next to the source file"""
    bl_idname = "xyz_molecule.render_figure"
    bl_label = "Render Figure"
    bl_options = {'REGISTER'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def execute(self, context):
        obj = active_molecule(context)
        tools = context.scene.xyz_mol_tools
        try:
            path = scene.render_figure(
                context, obj,
                resolution=tools.figure_resolution,
                background=tools.figure_background,
                high_quality=tools.figure_high_quality,
                frame_target=tools.figure_frame,
            )
        except Exception as exc:
            self.report({'ERROR'}, f"Render failed: {exc}")
            return {'CANCELLED'}
        self.report({'INFO'}, f"Figure saved: {path}")
        return {'FINISHED'}


# ---------------------------------------------------------------------------
# Edit-mode attribute stamp (kept from v1)
# ---------------------------------------------------------------------------

class MESH_OT_xyz_attribute_to_selection(Operator):
    """Add a named attribute and set its value on the selected vertices
    or edges (Edit Mode)"""
    bl_idname = "mesh.xyz_attribute_to_selection"
    bl_label = "Add Attribute to Selection"
    bl_options = {'REGISTER', 'UNDO'}

    attr_name: StringProperty(
        name="Name",
        description="Attribute name (created if it does not exist)",
        default="Tag",
    )
    domain: EnumProperty(
        name="Domain",
        items=[
            ('POINT', "Vertex", "Store the attribute on vertices"),
            ('EDGE', "Edge", "Store the attribute on edges"),
        ],
        default='POINT',
    )
    data_type: EnumProperty(
        name="Type",
        items=[
            ('FLOAT', "Float", "Floating-point attribute"),
            ('INT', "Integer", "Integer attribute"),
        ],
        default='FLOAT',
    )
    value_float: FloatProperty(name="Value", default=1.0)
    value_int: IntProperty(name="Value", default=1)
    only_selected: BoolProperty(
        name="Only Selected",
        description="Write the value only on selected elements "
                    "(otherwise on all elements)",
        default=True,
    )

    @classmethod
    def poll(cls, context):
        obj = context.edit_object
        return obj is not None and obj.type == 'MESH'

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.prop(self, "attr_name")
        layout.prop(self, "domain")
        layout.prop(self, "data_type")
        if self.data_type == 'FLOAT':
            layout.prop(self, "value_float")
        else:
            layout.prop(self, "value_int")
        layout.prop(self, "only_selected")

    def execute(self, context):
        name = self.attr_name.strip()
        if not name:
            self.report({'ERROR'}, "Attribute name must not be empty.")
            return {'CANCELLED'}

        obj = context.edit_object
        bm = bmesh.from_edit_mesh(obj.data)

        seq = bm.verts if self.domain == 'POINT' else bm.edges
        layer_pool = (
            seq.layers.float if self.data_type == 'FLOAT'
            else seq.layers.int
        )

        existing = obj.data.attributes.get(name)
        if existing is not None and (
            existing.domain != self.domain
            or existing.data_type != self.data_type
        ):
            self.report(
                {'ERROR'},
                f"Attribute '{name}' already exists as "
                f"{existing.data_type} on {existing.domain} - choose a "
                "different name or matching type/domain.",
            )
            return {'CANCELLED'}

        layer = layer_pool.get(name)
        if layer is None:
            layer = layer_pool.new(name)

        value = (self.value_float if self.data_type == 'FLOAT'
                 else self.value_int)
        count = 0
        for elem in seq:
            if not self.only_selected or elem.select:
                elem[layer] = value
                count += 1

        bmesh.update_edit_mesh(obj.data)
        target = "vertices" if self.domain == 'POINT' else "edges"
        self.report(
            {'INFO'}, f"Set '{name}' = {value} on {count} {target}."
        )
        return {'FINISHED'}


classes = (
    IMPORT_MESH_OT_xyz_molecule,
    IMPORT_MESH_OT_csv_points,
    XYZMOL_OT_reload,
    XYZMOL_OT_build_supercell,
    XYZMOL_OT_set_style,
    XYZMOL_OT_toggle_element,
    XYZMOL_OT_hide_hydrogens,
    XYZMOL_OT_show_all_elements,
    XYZMOL_OT_select_element,
    XYZMOL_OT_swap_instance,
    XYZMOL_OT_reset_instance,
    XYZMOL_OT_color_by_attribute,
    XYZMOL_OT_color_by_element,
    XYZMOL_OT_render_figure,
    MESH_OT_xyz_attribute_to_selection,
)
