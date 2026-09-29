"""Materials and scientific colormaps.

One atom material and one bond material serve every element: both read the
``ElementColor`` attribute (per instance for atoms, per mesh point for
bonds, which gives bonds a smooth two-tone gradient).

Colormaps (viridis, plasma, coolwarm) are applied CPU-side by rewriting the
``ElementColor`` attribute, so they work with the same materials and are
carried through bonds automatically.
"""

import bpy

from .elements import hex_to_linear_rgb, srgb_to_linear

# ---------------------------------------------------------------------------
# Colormaps (sRGB control points, evenly spaced 0..1)
# ---------------------------------------------------------------------------

COLORMAPS = {
    # Matplotlib viridis -- perceptually uniform, colorblind-friendly.
    "VIRIDIS": ["440154", "482878", "3E4A89", "31688E", "26828E",
                "1F9E89", "35B779", "6ECE58", "B5DE2B", "FDE725"],
    # Matplotlib plasma.
    "PLASMA": ["0D0887", "46039F", "7201A8", "9C179E", "BD3786",
               "D8576B", "ED7953", "FB9F3A", "FDCA26", "F0F921"],
    # Moreland coolwarm -- diverging; ideal for signed data (charges).
    "COOLWARM": ["3B4CC0", "6688EE", "88BBFF", "B8D0F9", "DDDDDD",
                 "F5C4AD", "F49A7B", "DD604D", "B40426"],
}

COLORMAP_ITEMS = [
    ('VIRIDIS', "Viridis", "Perceptually uniform sequential colormap"),
    ('PLASMA', "Plasma", "Perceptually uniform sequential colormap"),
    ('COOLWARM', "Coolwarm",
     "Diverging blue-white-red colormap (good for signed data like "
     "charges)"),
]


def sample_colormap(name, t):
    """Sample colormap ``name`` at ``t`` in [0, 1] -> linear (r, g, b)."""
    stops = COLORMAPS[name]
    t = min(max(t, 0.0), 1.0)
    pos = t * (len(stops) - 1)
    i = min(int(pos), len(stops) - 2)
    frac = pos - i
    c0 = hex_to_linear_rgb(stops[i])
    c1 = hex_to_linear_rgb(stops[i + 1])
    return tuple(c0[k] + (c1[k] - c0[k]) * frac for k in range(3))


# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------

def _make_attribute_material(mat_name, attribute_type):
    """Principled material driven by the 'ElementColor' attribute."""
    mat = bpy.data.materials.get(mat_name)
    if mat is not None:
        return mat

    mat = bpy.data.materials.new(mat_name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()

    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (400, 0)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (100, 0)
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-220, 0)
    attr.attribute_type = attribute_type
    attr.attribute_name = "ElementColor"

    nt.links.new(attr.outputs["Color"], bsdf.inputs["Base Color"])
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    bsdf.inputs["Roughness"].default_value = 0.35
    if "Specular IOR Level" in bsdf.inputs:
        bsdf.inputs["Specular IOR Level"].default_value = 0.5
    return mat


def create_atom_material(name):
    """Instanced-geometry material: reads the instancer's ElementColor."""
    return _make_attribute_material(f"{name}_AtomMaterial", 'INSTANCER')


def create_bond_material(name):
    """Bond material: reads the mesh ElementColor, which propagates through
    Mesh to Curve -> Curve to Mesh, giving each bond a smooth gradient."""
    return _make_attribute_material(f"{name}_BondMaterial", 'GEOMETRY')


# ---------------------------------------------------------------------------
# Color-by-attribute
# ---------------------------------------------------------------------------

def apply_colormap_to_mesh(mesh, attr_name, colormap='VIRIDIS',
                           vmin=None, vmax=None):
    """Rewrite 'ElementColor' from a scalar point attribute + colormap.

    Returns ``(vmin, vmax)`` actually used, or raises ValueError.
    """
    attr = mesh.attributes.get(attr_name)
    if attr is None or attr.domain != 'POINT':
        raise ValueError(f"No point attribute named '{attr_name}'.")
    if attr.data_type not in ('FLOAT', 'INT'):
        raise ValueError(
            f"Attribute '{attr_name}' is {attr.data_type}; color-by needs "
            "a scalar FLOAT or INT attribute."
        )

    n = len(attr.data)
    values = [0.0] * n
    attr.data.foreach_get("value", values)
    values = [float(v) for v in values]

    if vmin is None:
        vmin = min(values)
    if vmax is None:
        vmax = max(values)
    # A constant attribute (e.g. Atomic Number on a single-element slab)
    # has no range: paint the middle of the ramp rather than the dark end,
    # which would look like a failed render.
    constant = abs(vmax - vmin) < 1e-30
    span = 1.0 if constant else (vmax - vmin)

    color_attr = mesh.attributes.get("ElementColor")
    if color_attr is None:
        color_attr = mesh.attributes.new(
            "ElementColor", 'FLOAT_VECTOR', 'POINT'
        )

    colors = [0.0] * (3 * n)
    for i, v in enumerate(values):
        t = 0.5 if constant else (v - vmin) / span
        r, g, b = sample_colormap(colormap, t)
        colors[3 * i] = r
        colors[3 * i + 1] = g
        colors[3 * i + 2] = b
    color_attr.data.foreach_set("vector", colors)
    mesh.update()
    return vmin, vmax


def restore_element_colors(mesh, symbols, elements):
    """Reset 'ElementColor' to the per-element CPK colors."""
    n = len(symbols)
    colors = [0.0] * (3 * n)
    for i, sym in enumerate(symbols):
        entry = elements.get(sym) or elements["X"]
        rgb = hex_to_linear_rgb(entry["color"])
        colors[3 * i: 3 * i + 3] = rgb[0:3]
    color_attr = mesh.attributes.get("ElementColor")
    if color_attr is None:
        color_attr = mesh.attributes.new(
            "ElementColor", 'FLOAT_VECTOR', 'POINT'
        )
    color_attr.data.foreach_set("vector", colors)
    mesh.update()


# ---------------------------------------------------------------------------
# World (studio lighting)
# ---------------------------------------------------------------------------

WORLD_NAME = "XYZ_Molecule_World"


def make_world(context):
    """Create (once) a blackbody + rim-highlight world and assign it."""
    if WORLD_NAME in bpy.data.worlds:
        context.scene.world = bpy.data.worlds[WORLD_NAME]
        return

    world = bpy.data.worlds.new(WORLD_NAME)
    world.use_nodes = True
    context.scene.world = world

    nodes = world.node_tree.nodes
    links = world.node_tree.links
    nodes.clear()

    w, h = 200, 100

    out = nodes.new("ShaderNodeOutputWorld")
    out.location = (4 * w, 0)

    blackbody = nodes.new("ShaderNodeBlackbody")
    blackbody.location = (0, 0)
    blackbody.inputs[0].default_value = 4000  # Kelvin

    bg_main = nodes.new("ShaderNodeBackground")
    bg_main.location = (w, 0)
    bg_main.inputs["Strength"].default_value = 1.0
    links.new(blackbody.outputs[0], bg_main.inputs["Color"])

    bg_highlight = nodes.new("ShaderNodeBackground")
    bg_highlight.location = (w, -2 * h)
    bg_highlight.inputs["Color"].default_value = (1.0, 0.863, 0.680, 1.0)
    bg_highlight.inputs["Strength"].default_value = 15.0

    tex_coord = nodes.new("ShaderNodeTexCoord")
    tex_coord.location = (-3 * w, 3 * h)

    normalize = nodes.new("ShaderNodeVectorMath")
    normalize.operation = 'NORMALIZE'
    normalize.location = (-2 * w, 3 * h)

    normal = nodes.new("ShaderNodeNormal")
    normal.location = (-2 * w, h)
    normal.outputs[0].default_value = (0.55, -0.558, 0.621)

    dot = nodes.new("ShaderNodeVectorMath")
    dot.operation = 'DOT_PRODUCT'
    dot.location = (-w, 3 * h)

    map_range = nodes.new("ShaderNodeMapRange")
    map_range.location = (0, 3 * h)
    map_range.inputs["From Min"].default_value = 0.85

    mixer = nodes.new("ShaderNodeMixShader")
    mixer.location = (2 * w, 0)

    links.new(tex_coord.outputs["Generated"], normalize.inputs[0])
    links.new(normalize.outputs["Vector"], dot.inputs[0])
    links.new(normal.outputs["Normal"], dot.inputs[1])
    links.new(dot.outputs["Value"], map_range.inputs["Value"])
    links.new(map_range.outputs["Result"], mixer.inputs["Fac"])
    links.new(bg_main.outputs["Background"], mixer.inputs[1])
    links.new(bg_highlight.outputs["Background"], mixer.inputs[2])
    links.new(mixer.outputs["Shader"], out.inputs["Surface"])
