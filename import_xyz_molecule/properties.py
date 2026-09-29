"""Scene-level tool settings shown in the sidebar (swap target, figure
options, color-by-attribute options) plus the dynamic enum callbacks that
list the active molecule's elements and attributes.
"""

import json

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)
from bpy.types import PropertyGroup

from . import builder, nodes
from .materials import COLORMAP_ITEMS

# Dynamic EnumProperty items must stay referenced from Python, or Blender
# may read freed strings; keep the last result alive here.
_element_items_cache = []
_attribute_items_cache = []

# Attributes that are pure UI/machinery and never worth coloring by.
# Everything else -- including the built-in physical ones below and any
# column imported from an extxyz/CSV file -- is offered to the user.
_HIDDEN_ATTRIBUTES = {
    "ElementColor",     # this is the colormap's output, not an input
    "ElementHidden",    # per-element visibility state
    "position",
    ".select_vert", ".select_edge", ".select_poly",
}

# Friendlier labels for the attributes the importer writes itself, in the
# order they should appear. Anything not listed keeps its raw name and is
# listed afterwards (those are the user's own imported columns).
_BUILTIN_ATTRIBUTE_LABELS = [
    ("AtomicNumber", "Atomic Number", "Element's atomic number Z"),
    ("AtomRadius", "Covalent Radius",
     "Per-atom covalent radius in Angstrom"),
    ("VdWRadius", "van der Waals Radius",
     "Per-atom van der Waals radius in Angstrom"),
    ("IndexOfElement", "Element Index",
     "Index into this molecule's element list (categorical coloring "
     "by element)"),
    ("RandomValue", "Random Value",
     "Per-atom random number, for artistic variation"),
    ("IsGhost", "Is Periodic Image",
     "1 for ghost atoms added across cell boundaries, 0 otherwise"),
]


def active_molecule(context):
    obj = context.active_object
    if builder.is_molecule(obj):
        return obj
    # Fall back to the molecule a selected child (e.g. unit cell) points to.
    if obj is not None and builder.is_molecule(obj.parent):
        return obj.parent
    return None


def element_enum_items(self, context):
    global _element_items_cache
    obj = active_molecule(context)
    items = []
    if obj is not None:
        for sym in builder.get_molecule_elements(obj):
            items.append((sym, sym, f"Element {sym}"))
    if not items:
        items = [('NONE', "(no molecule)", "")]
    _element_items_cache = items
    return items


def scalar_attribute_enum_items(self, context):
    """Every scalar point attribute worth coloring by.

    Built-in physical attributes come first under readable names, then
    any extra columns imported from the file (charges, forces, spins...).
    """
    global _attribute_items_cache
    obj = active_molecule(context)
    items = []
    if obj is not None and obj.data is not None:
        available = {
            attr.name for attr in obj.data.attributes
            if attr.domain == 'POINT'
            and attr.data_type in ('FLOAT', 'INT')
            and attr.name not in _HIDDEN_ATTRIBUTES
            and not attr.name.startswith(".")
        }
        for name, label, description in _BUILTIN_ATTRIBUTE_LABELS:
            if name in available:
                items.append((name, label, description))
                available.discard(name)
        for name in sorted(available):
            items.append((name, name,
                          f"Imported per-atom column '{name}'"))
    if not items:
        items = [('NONE', "(no molecule selected)",
                  "Select an imported molecule to list its attributes")]
    _attribute_items_cache = items
    return items


def element_counts(obj):
    try:
        return json.loads(obj["xyz_counts"])
    except (KeyError, TypeError, ValueError):
        return {}


# ---------------------------------------------------------------------------
# Per-object display settings (the sidebar's sliders and toggles)
# ---------------------------------------------------------------------------
#
# These are real add-on properties stored on the molecule object. Each one
# pushes its value into the Geometry Nodes modifier when changed. Drawing
# these -- instead of the modifier's own ID-properties -- means the sidebar
# never depends on when Blender chooses to create those ID-properties,
# which differs between 4.x and 5.x and left the controls blank.

# (property name, node group input socket name)
DISPLAY_SOCKETS = (
    ("atom_scale", "Atom Scale"),
    ("bond_radius", "Bond Radius"),
    ("show_atoms", "Show Atoms"),
    ("show_bonds", "Show Bonds"),
    ("show_ghosts", "Show Ghost Atoms"),
    ("random_rotation", "Random Rotation"),
    ("scale_jitter", "Scale Jitter"),
    ("seed", "Seed"),
    ("align_to_vector", "Align to Vector"),
    ("vector_attribute", "Vector Attribute"),
)

# While True, property updates do not push to the modifier (used when
# seeding many values at once, so the tree is not re-evaluated per value).
_seeding = False


def _push(self, context, prop_name, socket_name):
    """Update callback: write one display property into the modifier."""
    if _seeding:
        return
    obj = self.id_data
    mod = nodes.find_molecule_modifier(obj)
    if mod is None:
        return
    nodes.set_gn_input(mod, socket_name, getattr(self, prop_name))
    obj.update_tag()
    try:
        context.view_layer.update()
    except AttributeError:
        pass


def _updater(prop_name, socket_name):
    return lambda self, context: _push(self, context, prop_name, socket_name)


class XYZMoleculeDisplaySettings(PropertyGroup):
    """Display sliders/toggles for one molecule, mirrored into its
    Geometry Nodes modifier."""

    atom_scale: FloatProperty(
        name="Atom Scale",
        description="Sphere radius as a fraction of the element radius",
        default=0.7, min=0.0, soft_max=2.0,
        update=_updater("atom_scale", "Atom Scale"),
    )
    bond_radius: FloatProperty(
        name="Bond Radius",
        description="Radius of the bond cylinders (Angstrom)",
        default=0.1, min=0.0, soft_max=1.0, precision=3,
        update=_updater("bond_radius", "Bond Radius"),
    )
    show_atoms: BoolProperty(
        name="Atoms", description="Display atom spheres", default=True,
        update=_updater("show_atoms", "Show Atoms"),
    )
    show_bonds: BoolProperty(
        name="Bonds", description="Display bonds", default=True,
        update=_updater("show_bonds", "Show Bonds"),
    )
    show_ghosts: BoolProperty(
        name="Ghost Atoms (periodic)",
        description="Display periodic-image atoms at cell boundaries",
        default=True,
        update=_updater("show_ghosts", "Show Ghost Atoms"),
    )
    random_rotation: FloatProperty(
        name="Random Rotation",
        description="Random per-instance rotation "
                    "(0 = none, 1 = fully random)",
        default=0.0, min=0.0, max=1.0, subtype='FACTOR',
        update=_updater("random_rotation", "Random Rotation"),
    )
    scale_jitter: FloatProperty(
        name="Scale Jitter",
        description="Random per-instance scale variation",
        default=0.0, min=0.0, max=1.0, subtype='FACTOR',
        update=_updater("scale_jitter", "Scale Jitter"),
    )
    seed: IntProperty(
        name="Seed",
        description="Random seed for rotation/scale variation",
        default=0, min=0,
        update=_updater("seed", "Seed"),
    )
    align_to_vector: BoolProperty(
        name="Align to Vector",
        description="Align each instance's +Z axis to the vector "
                    "attribute below",
        default=False,
        update=_updater("align_to_vector", "Align to Vector"),
    )
    vector_attribute: StringProperty(
        name="Attribute",
        description="Per-atom vector attribute used for alignment "
                    "(BondDir = average bond direction; or an imported "
                    "column such as forces)",
        default="BondDir",
        update=_updater("vector_attribute", "Vector Attribute"),
    )


def sync_display_to_modifier(obj):
    """Push every display setting of ``obj`` into its modifier (one pass,
    no depsgraph update per value)."""
    display = getattr(obj, "xyz_mol_display", None)
    mod = nodes.find_molecule_modifier(obj)
    if display is None or mod is None:
        return False
    for prop_name, socket_name in DISPLAY_SOCKETS:
        nodes.set_gn_input(mod, socket_name, getattr(display, prop_name))
    obj.update_tag()
    return True


def seed_display_settings(obj, **values):
    """Set display settings without per-value pushes, then push once.

    Used right after import so the modifier and the sidebar agree from
    the first frame, whatever Blender did (or did not) seed itself.
    """
    global _seeding
    display = getattr(obj, "xyz_mol_display", None)
    if display is None:
        return False
    _seeding = True
    try:
        for prop_name, value in values.items():
            setattr(display, prop_name, value)
    finally:
        _seeding = False
    return sync_display_to_modifier(obj)


class XYZMoleculeToolSettings(PropertyGroup):
    """Sidebar tool settings (stored per scene)."""

    # ---- instance swapping ----------------------------------------------
    swap_element: EnumProperty(
        name="Element",
        description="Element whose instances will be replaced",
        items=element_enum_items,
    )
    swap_object: PointerProperty(
        name="Object",
        description="Any object in the scene to instance for this element "
                    "(spheres, cats, trees, ...)",
        type=bpy.types.Object,
    )

    # ---- color by attribute ---------------------------------------------
    color_attribute: EnumProperty(
        name="Attribute",
        description="Scalar per-atom attribute to color by "
                    "(e.g. charge from an extxyz file)",
        items=scalar_attribute_enum_items,
    )
    colormap: EnumProperty(
        name="Colormap",
        items=COLORMAP_ITEMS,
        default='VIRIDIS',
    )
    color_auto_range: BoolProperty(
        name="Auto Range",
        description="Use the attribute's min/max as the colormap range",
        default=True,
    )
    color_min: FloatProperty(name="Min", default=0.0)
    color_max: FloatProperty(name="Max", default=1.0)

    # ---- figure rendering ------------------------------------------------
    figure_resolution: IntProperty(
        name="Resolution",
        description="Pixels along the longer image side; the shorter "
                    "side follows the structure's footprint",
        default=1600, min=64, max=16384,
    )
    figure_frame: EnumProperty(
        name="Frame",
        description="What the top-down camera frames",
        items=[
            ('WHOLE', "Whole Object",
             "Every visible atom, slab included"),
            ('ADSORBATE', "Adsorbate Only",
             "Only atoms that are not slab elements (the Slab Elements "
             "list from the import); the whole object if none are left"),
        ],
        default='WHOLE',
    )
    figure_background: EnumProperty(
        name="Background",
        items=[
            ('TRANSPARENT', "Transparent",
             "PNG with alpha channel (drop onto any background)"),
            ('WHITE', "White", "Solid white background"),
        ],
        default='TRANSPARENT',
    )
    figure_high_quality: BoolProperty(
        name="High Quality (Cycles)",
        description="Path-traced render; slower but publication-grade",
        default=False,
    )


def register():
    bpy.utils.register_class(XYZMoleculeDisplaySettings)
    bpy.utils.register_class(XYZMoleculeToolSettings)
    bpy.types.Object.xyz_mol_display = PointerProperty(
        type=XYZMoleculeDisplaySettings
    )
    bpy.types.Scene.xyz_mol_tools = PointerProperty(
        type=XYZMoleculeToolSettings
    )


def unregister():
    del bpy.types.Scene.xyz_mol_tools
    del bpy.types.Object.xyz_mol_display
    bpy.utils.unregister_class(XYZMoleculeToolSettings)
    bpy.utils.unregister_class(XYZMoleculeDisplaySettings)
