"""Sidebar panels, add-on preferences, and File > Import menu entries."""

from bpy.props import StringProperty
from bpy.types import AddonPreferences, Panel

from . import builder, nodes, operators
from .properties import active_molecule, element_counts


class XYZMoleculeAddonPreferences(AddonPreferences):
    bl_idname = __package__

    user_elements_path: StringProperty(
        name="Element Overrides (JSON)",
        description=(
            "Default JSON file with element overrides, e.g. "
            '{"Au": {"radius": 150, "color": "FFD700"}}'
        ),
        subtype='FILE_PATH',
        default="",
    )
    default_slab_elements: StringProperty(
        name="Default Slab Elements",
        description="Comma-separated element symbols treated as slab atoms",
        default="Au,Cu,Pt,Ga,Si",
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "user_elements_path")
        layout.prop(self, "default_slab_elements")
        box = layout.box()
        box.label(text="JSON override format:", icon='INFO')
        box.label(text='{"Au": {"radius": 150, "vdw_radius": 214,')
        box.label(text='         "color": "FFD700", "object": "MyMesh"},')
        box.label(text=' "Cu": {"color": "AA4400"}}')
        box.label(text="Radius in picometers, color as sRGB hex; "
                       "\"object\" names a scene object to instance.")


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

class XYZMoleculePanelMixin:
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "XYZ Molecule"


def _display_settings(obj):
    """The molecule's sidebar display settings (None if unregistered)."""
    return getattr(obj, "xyz_mol_display", None)


class VIEW3D_PT_xyz_import(XYZMoleculePanelMixin, Panel):
    bl_label = "Import"
    bl_idname = "VIEW3D_PT_xyz_import"

    def draw(self, context):
        col = self.layout.column(align=True)
        col.operator("import_mesh.xyz_molecule",
                     text="Import XYZ File", icon='FILEBROWSER')
        col.operator("import_mesh.xyz_csv_points",
                     text="Import Points CSV", icon='SPREADSHEET')


class VIEW3D_PT_xyz_display(XYZMoleculePanelMixin, Panel):
    bl_label = "Display"
    bl_idname = "VIEW3D_PT_xyz_display"

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def draw(self, context):
        layout = self.layout
        obj = active_molecule(context)
        mod = nodes.find_molecule_modifier(obj)
        current_style = obj.get("xyz_style", 'BALL_AND_STICK')

        col = layout.column(align=True)
        for identifier, label, _desc in nodes.STYLE_ITEMS:
            op = col.operator("xyz_molecule.set_style", text=label,
                              depress=(identifier == current_style))
            op.style = identifier

        display = _display_settings(obj)
        if mod is not None and display is not None:
            col = layout.column(align=True)
            col.prop(display, "atom_scale")
            col.prop(display, "bond_radius")
            row = layout.row(align=True)
            row.prop(display, "show_atoms")
            row.prop(display, "show_bonds")
            layout.prop(display, "show_ghosts")

        layout.separator()
        col = layout.column(align=True)
        col.operator("xyz_molecule.reload", icon='FILE_REFRESH')

        repeats = operators.normalize_repeats_for_ui(obj)
        if operators.has_lattice(obj):
            nx, ny, nz = repeats
            col.operator(
                "xyz_molecule.build_supercell",
                text=(f"Supercell ({nx}x{ny}x{nz})"
                      if repeats != (1, 1, 1) else "Build Supercell"),
                icon='MOD_ARRAY',
            )


class VIEW3D_PT_xyz_elements(XYZMoleculePanelMixin, Panel):
    bl_label = "Elements"
    bl_idname = "VIEW3D_PT_xyz_elements"

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def draw(self, context):
        layout = self.layout
        obj = active_molecule(context)
        elements_in_file = builder.get_molecule_elements(obj)
        hidden = builder.get_hidden_elements(obj)
        counts = element_counts(obj)

        row = layout.row(align=True)
        row.operator("xyz_molecule.hide_hydrogens",
                     text="Show Hydrogens" if "H" in hidden
                     else "Hide Hydrogens",
                     icon='HIDE_ON' if "H" in hidden else 'HIDE_OFF')
        row.operator("xyz_molecule.show_all_elements", text="Show All")

        col = layout.column(align=True)
        for sym in elements_in_file:
            row = col.row(align=True)
            count = counts.get(sym)
            row.label(text=f"{sym}  ({count})" if count else sym)
            op = row.operator(
                "xyz_molecule.toggle_element", text="",
                icon='HIDE_ON' if sym in hidden else 'HIDE_OFF',
            )
            op.element = sym
            op = row.operator("xyz_molecule.select_element", text="",
                              icon='RESTRICT_SELECT_OFF')
            op.element = sym
            op.extend = False


class VIEW3D_PT_xyz_figure(XYZMoleculePanelMixin, Panel):
    bl_label = "Figure"
    bl_idname = "VIEW3D_PT_xyz_figure"

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def draw(self, context):
        layout = self.layout
        tools = context.scene.xyz_mol_tools
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(tools, "figure_frame")
        layout.prop(tools, "figure_resolution")
        layout.prop(tools, "figure_background")
        layout.prop(tools, "figure_high_quality")
        layout.operator("xyz_molecule.render_figure",
                        icon='RENDER_STILL')


class VIEW3D_PT_xyz_color(XYZMoleculePanelMixin, Panel):
    bl_label = "Color"
    bl_idname = "VIEW3D_PT_xyz_color"

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def draw(self, context):
        layout = self.layout
        tools = context.scene.xyz_mol_tools
        layout.use_property_split = True
        layout.use_property_decorate = False
        layout.prop(tools, "color_attribute")
        layout.prop(tools, "colormap")
        layout.prop(tools, "color_auto_range")
        sub = layout.column(align=True)
        sub.enabled = not tools.color_auto_range
        sub.prop(tools, "color_min")
        sub.prop(tools, "color_max")
        row = layout.row(align=True)
        row.operator("xyz_molecule.color_by_attribute",
                     text="Apply", icon='COLOR')
        row.operator("xyz_molecule.color_by_element",
                     text="By Element")


class VIEW3D_PT_xyz_instances(XYZMoleculePanelMixin, Panel):
    bl_label = "Instances"
    bl_idname = "VIEW3D_PT_xyz_instances"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return active_molecule(context) is not None

    def draw(self, context):
        layout = self.layout
        obj = active_molecule(context)
        tools = context.scene.xyz_mol_tools
        mod = nodes.find_molecule_modifier(obj)

        col = layout.column(align=True)
        col.label(text="Swap Instance", icon='OUTLINER_OB_GROUP_INSTANCE')
        col.prop(tools, "swap_element")
        col.prop(tools, "swap_object")
        row = col.row(align=True)
        row.operator("xyz_molecule.swap_instance", text="Swap")
        row.operator("xyz_molecule.reset_instance", text="Reset")

        display = _display_settings(obj)
        if mod is not None and display is not None:
            layout.separator()
            col = layout.column(align=True)
            col.label(text="Variation", icon='FORCE_TURBULENCE')
            col.prop(display, "random_rotation")
            col.prop(display, "scale_jitter")
            col.prop(display, "seed")

            layout.separator()
            col = layout.column(align=True)
            col.label(text="Orientation", icon='ORIENTATION_GIMBAL')
            col.prop(display, "align_to_vector")
            sub = col.column(align=True)
            sub.enabled = display.align_to_vector
            sub.prop_search(display, "vector_attribute",
                            obj.data, "attributes", text="Attribute")


class VIEW3D_PT_xyz_attributes(XYZMoleculePanelMixin, Panel):
    bl_label = "Attributes (Edit Mode)"
    bl_idname = "VIEW3D_PT_xyz_attributes"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        layout.operator("mesh.xyz_attribute_to_selection", icon='PLUS')
        if context.mode != 'EDIT_MESH':
            layout.label(text="Enter Edit Mode to tag a selection",
                         icon='INFO')


# ---------------------------------------------------------------------------
# Menu entries
# ---------------------------------------------------------------------------

def menu_func_import(self, context):
    self.layout.operator("import_mesh.xyz_molecule",
                         text="XYZ Molecule (.xyz)")
    self.layout.operator("import_mesh.xyz_csv_points",
                         text="Points CSV (.csv)")


classes = (
    XYZMoleculeAddonPreferences,
    VIEW3D_PT_xyz_import,
    VIEW3D_PT_xyz_display,
    VIEW3D_PT_xyz_elements,
    VIEW3D_PT_xyz_figure,
    VIEW3D_PT_xyz_color,
    VIEW3D_PT_xyz_instances,
    VIEW3D_PT_xyz_attributes,
)
