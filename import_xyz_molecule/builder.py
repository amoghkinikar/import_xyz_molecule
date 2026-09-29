"""Molecule construction: mesh + attributes, template collection, unit
cell, custom properties. Ties the parser, bond, node, and material modules
together for both first import and reload.
"""

import json
import random

import bmesh
import bpy

from . import bonds as bonds_mod
from . import materials, nodes, trajectory
from .elements import hex_to_linear_rgb

TEMPLATE_ROOT_NAME = "XYZ_Molecule_Templates"

# Object custom-property keys.
PROP_FLAG = "xyz_molecule"
PROP_SOURCE = "xyz_source"
PROP_KIND = "xyz_kind"              # 'XYZ' or 'CSV'
PROP_SETTINGS = "xyz_settings"      # JSON dict of import settings
PROP_ELEMENTS = "xyz_elements"      # JSON list of element symbols (sorted)
PROP_HIDDEN = "xyz_hidden"          # JSON list of hidden element symbols
PROP_OFFSET = "xyz_offset"          # centering offset (Angstrom, xyz tuple)
PROP_LATTICE = "xyz_lattice"        # JSON 3x3 lattice or ""
PROP_COUNTS = "xyz_counts"          # JSON {symbol: atom count} (real atoms)
PROP_NAME = "xyz_name"              # base name used for datablocks
PROP_TABLE = "xyz_element_table"    # JSON {symbol: element entry} used here
PROP_REPEATS = "xyz_supercell"      # JSON [nx, ny, nz] supercell repeats
PROP_FRAMES = "xyz_n_frames"        # number of animated trajectory frames


def is_molecule(obj):
    return obj is not None and obj.type == 'MESH' and obj.get(PROP_FLAG)


def get_molecule_elements(obj):
    try:
        return json.loads(obj[PROP_ELEMENTS])
    except (KeyError, TypeError, ValueError):
        return []


def get_hidden_elements(obj):
    try:
        return set(json.loads(obj[PROP_HIDDEN]))
    except (KeyError, TypeError, ValueError):
        return set()


def set_hidden_elements(obj, hidden):
    obj[PROP_HIDDEN] = json.dumps(sorted(hidden))


# ---------------------------------------------------------------------------
# Template collection (per-element instance objects)
# ---------------------------------------------------------------------------

def _find_layer_collection(layer_col, collection):
    if layer_col.collection == collection:
        return layer_col
    for child in layer_col.children:
        found = _find_layer_collection(child, collection)
        if found:
            return found
    return None


def get_template_root_collection(context):
    """Hidden top-level collection that holds all per-molecule templates."""
    col = bpy.data.collections.get(TEMPLATE_ROOT_NAME)
    if col is None:
        col = bpy.data.collections.new(TEMPLATE_ROOT_NAME)
    if col.name not in context.scene.collection.children:
        context.scene.collection.children.link(col)
    layer_col = _find_layer_collection(
        context.view_layer.layer_collection, col
    )
    if layer_col is not None:
        layer_col.exclude = True
    return col


def template_object_name(mol_name, index, symbol):
    """Zero-padded index prefix keeps the collection's (name-sorted) order
    identical to the ``IndexOfElement`` attribute, even after swaps."""
    return f"{mol_name}_atom_{index:03d}_{symbol}"


def find_template_object(collection, index):
    prefix_end = f"_atom_{index:03d}_"
    for obj in collection.objects:
        if prefix_end in obj.name:
            return obj
    return None


def make_sphere_object(name, atom_material, location):
    mesh = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=1.0)
    bm.to_mesh(mesh)
    bm.free()
    for poly in mesh.polygons:
        poly.use_smooth = True
    mesh.materials.append(atom_material)
    obj = bpy.data.objects.new(name, mesh)
    obj.location = location
    return obj


def make_atom_collection(context, mol_name, elements_in_file, elements,
                         atom_material):
    """One template object per element (radius-1 sphere by default, or the
    object named by the element's JSON "object" key) in a hidden
    collection."""
    root = get_template_root_collection(context)

    col = bpy.data.collections.new(f"{mol_name}_Atoms")
    root.children.link(col)

    for i, symbol in enumerate(elements_in_file):
        name = template_object_name(mol_name, i, symbol)
        entry = elements.get(symbol, {})
        source_name = entry.get("object", "")
        source = bpy.data.objects.get(source_name) if source_name else None
        if source is not None:
            obj = source.copy()      # linked copy: shares mesh data
            obj.parent = None
            obj.name = name
        else:
            obj = make_sphere_object(
                name, atom_material, (i * 3.0, 0.0, 0.0)
            )
        col.objects.link(obj)
        obj.hide_render = True
    return col


def find_atom_collection(obj):
    mod = nodes.find_molecule_modifier(obj)
    if mod is None or mod.node_group is None:
        return None
    for node in mod.node_group.nodes:
        if node.bl_idname == "GeometryNodeCollectionInfo":
            return node.inputs["Collection"].default_value
    return None


# ---------------------------------------------------------------------------
# Supercell
# ---------------------------------------------------------------------------

MAX_SUPERCELL_ATOMS = 200000


def normalize_repeats(value):
    """Coerce a stored/passed repeat triple to ``(nx, ny, nz)`` ints >= 1."""
    try:
        nx, ny, nz = (max(1, int(v)) for v in value)
    except (TypeError, ValueError):
        return (1, 1, 1)
    return (nx, ny, nz)


def expand_supercell(structure, repeats):
    """Replicate a periodic structure along its lattice vectors.

    Returns a new structure dict whose ``lattice`` is the enlarged
    supercell and whose per-atom ``extra`` arrays are replicated to
    match. The input structure is not modified.

    Raises ValueError when the structure carries no lattice, or when the
    requested supercell would exceed ``MAX_SUPERCELL_ATOMS``.
    """
    nx, ny, nz = normalize_repeats(repeats)
    if (nx, ny, nz) == (1, 1, 1):
        return structure

    lattice = structure.get("lattice")
    if lattice is None:
        raise ValueError(
            "This structure has no unit cell, so it cannot be repeated. "
            "Supercells need an extended-xyz file with a Lattice=\"...\" "
            "entry in its comment line."
        )

    positions = structure["positions"]
    symbols = structure["symbols"]
    total = len(positions) * nx * ny * nz
    if total > MAX_SUPERCELL_ATOMS:
        raise ValueError(
            f"A {nx}x{ny}x{nz} supercell of this structure would have "
            f"{total:,} atoms (limit {MAX_SUPERCELL_ATOMS:,}). "
            "Use smaller repeats."
        )

    a, b, c = lattice
    new_positions = []
    new_symbols = []
    source_index = []          # original atom each copy came from
    for ia in range(nx):
        for ib in range(ny):
            for ic in range(nz):
                sx = ia * a[0] + ib * b[0] + ic * c[0]
                sy = ia * a[1] + ib * b[1] + ic * c[1]
                sz = ia * a[2] + ib * b[2] + ic * c[2]
                for i, (x, y, z) in enumerate(positions):
                    new_positions.append((x + sx, y + sy, z + sz))
                    new_symbols.append(symbols[i])
                    source_index.append(i)

    new_structure = dict(structure)
    new_structure["symbols"] = new_symbols
    new_structure["positions"] = new_positions
    new_structure["lattice"] = (
        tuple(v * nx for v in a),
        tuple(v * ny for v in b),
        tuple(v * nz for v in c),
    )
    new_structure["extra"] = {
        name: {"type": spec["type"],
               "values": [spec["values"][i] for i in source_index]}
        for name, spec in structure.get("extra", {}).items()
    }
    return new_structure


# ---------------------------------------------------------------------------
# Unit cell
# ---------------------------------------------------------------------------

def make_unit_cell_object(context, mol_name, lattice, offset,
                          parent, collection, world_scale,
                          repeats=(1, 1, 1)):
    """Wireframe cell box, drawn as a grid of ``repeats`` primitive cells.

    ``lattice`` is the full (super)cell; dividing it by ``repeats`` gives
    the primitive vectors, so a 3x3x1 supercell shows all nine original
    cells rather than one big box. With repeats (1, 1, 1) this is the
    familiar 8-vertex, 12-edge box.
    """
    verts, edges = unit_cell_geometry(lattice, offset, repeats)
    mesh = bpy.data.meshes.new(f"{mol_name}_cell")
    mesh.from_pydata(verts, edges, [])
    mesh.update()

    obj = bpy.data.objects.new(f"{mol_name}_cell", mesh)
    obj.display_type = 'WIRE'
    collection.objects.link(obj)
    obj.parent = parent
    # Parent inverse identity: the cell simply inherits the molecule's
    # transform (including world scale), matching the atom coordinates.
    return obj


def unit_cell_geometry(lattice, offset, repeats=(1, 1, 1)):
    """Vertices and edges of the cell grid (see make_unit_cell_object)."""
    nx, ny, nz = normalize_repeats(repeats)
    a = tuple(v / nx for v in lattice[0])
    b = tuple(v / ny for v in lattice[1])
    c = tuple(v / nz for v in lattice[2])
    ox, oy, oz = offset

    # Shared corner lattice: (nx+1) x (ny+1) x (nz+1) points, so interior
    # cell walls are drawn once instead of twice.
    def index(i, j, k):
        return (i * (ny + 1) + j) * (nz + 1) + k

    verts = []
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                verts.append((
                    i * a[0] + j * b[0] + k * c[0] - ox,
                    i * a[1] + j * b[1] + k * c[1] - oy,
                    i * a[2] + j * b[2] + k * c[2] - oz,
                ))

    edges = []
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                here = index(i, j, k)
                if i < nx:
                    edges.append((here, index(i + 1, j, k)))
                if j < ny:
                    edges.append((here, index(i, j + 1, k)))
                if k < nz:
                    edges.append((here, index(i, j, k + 1)))
    return verts, edges


# ---------------------------------------------------------------------------
# Attribute writing
# ---------------------------------------------------------------------------

def _write_point_attribute(mesh, name, attr_type, values):
    attr = mesh.attributes.get(name)
    if attr is not None and (attr.data_type != attr_type
                             or attr.domain != 'POINT'):
        mesh.attributes.remove(attr)
        attr = None
    if attr is None:
        attr = mesh.attributes.new(name, attr_type, 'POINT')
    if attr_type == 'FLOAT_VECTOR':
        flat = [c for v in values for c in v]
        attr.data.foreach_set("vector", flat)
    elif attr_type == 'INT':
        attr.data.foreach_set("value", [int(v) for v in values])
    elif attr_type == 'BOOLEAN':
        attr.data.foreach_set("value", [bool(v) for v in values])
    else:
        attr.data.foreach_set("value", [float(v) for v in values])
    return attr


def _write_edge_float(mesh, name, values):
    attr = mesh.attributes.get(name)
    if attr is not None and (attr.data_type != 'FLOAT'
                             or attr.domain != 'EDGE'):
        mesh.attributes.remove(attr)
        attr = None
    if attr is None:
        attr = mesh.attributes.new(name, 'FLOAT', 'EDGE')
    attr.data.foreach_set("value", [float(v) for v in values])
    return attr


def build_molecule_mesh(structure, elements, mol_name, settings,
                        frames=None):
    """Create the mesh datablock with all point attributes. Returns
    ``(mesh, info)`` where ``info`` carries derived data for the caller.

    ``frames``: optional list of per-frame structures (same atoms) for an
    animated trajectory; ``structure`` is then the first of them. Bonds
    and ghosts are the union over all frames (see ``trajectory``), and
    ``info["frame_verts"]`` holds the per-frame vertex positions for the
    shape keys.
    """
    frames = frames or [structure]
    symbols = structure["symbols"]
    lattices = trajectory.frame_lattices(frames)
    frames_positions = [f["positions"] for f in frames]
    pbc = structure.get("pbc")
    if (len(frames) > 1 and lattices[0] is not None
            and settings.get("unwrap", True)):
        frames_positions = trajectory.unwrap_positions(
            frames_positions, lattices, pbc or (True, True, True)
        )
    positions = frames_positions[0]
    n_real = len(symbols)

    # ---- centering -------------------------------------------------------
    if settings.get("center", True) and n_real:
        xs = [p[0] for p in positions]
        ys = [p[1] for p in positions]
        zs = [p[2] for p in positions]
        offset = (
            (min(xs) + max(xs)) / 2.0,
            (min(ys) + max(ys)) / 2.0,
            (min(zs) + max(zs)) / 2.0,
        )
    else:
        offset = (0.0, 0.0, 0.0)

    # ---- radii for bonding ----------------------------------------------
    slab_elements = set(settings.get("slab_elements") or [])
    scale_slab = settings.get("scale_slab", False)
    slab_scale = settings.get("slab_scale", 2.0)
    bond_scale = settings.get("bond_scale", 1.1)

    def entry_for(sym):
        return elements.get(sym) or elements["X"]

    # Two separate radius lists. Bond detection must always use the true
    # covalent radii: letting the slab scale inflate the bonding cutoff
    # makes every metal atom bond to its whole neighbourhood (a 2x scale
    # on a Cu slab turns ~4,000 bonds into ~20,000). The slab scale is a
    # presentation setting, so it only affects the drawn sphere size.
    bond_radii = []
    display_radii = []
    for sym in symbols:
        r = entry_for(sym)["radius"] * 0.01     # pm -> Angstrom
        bond_radii.append(r)
        if scale_slab and sym in slab_elements:
            display_radii.append(r * slab_scale)
        else:
            display_radii.append(r)

    # ---- bonds (plus periodic ghosts when a lattice is present) ---------
    # For a single frame this is exactly that frame's bonds and ghosts;
    # for a trajectory it is the union over all frames.
    assembled = trajectory.assemble(
        frames_positions, lattices, bond_radii, pbc=pbc, scale=bond_scale,
        detect_bonds=settings.get("detect_bonds", True),
        periodic_bonds=settings.get("periodic_bonds", True),
    )
    ghost_sources = assembled["ghost_sources"]
    edges = assembled["edges"]
    symbols = symbols + [symbols[s] for s in ghost_sources]
    display_radii = display_radii + [display_radii[s] for s in ghost_sources]
    positions = assembled["positions"][0]

    n_total = len(positions)

    # ---- mesh ------------------------------------------------------------
    def shifted(frame_positions):
        return [(x - offset[0], y - offset[1], z - offset[2])
                for x, y, z in frame_positions]

    verts = shifted(positions)
    mesh = bpy.data.meshes.new(f"{mol_name}_mesh")
    mesh.from_pydata(verts, edges, [])
    mesh.update()
    # Per-edge bond cutoff: lets the node tree drop bonds that have
    # stretched beyond it at the current frame (trajectories).
    _write_edge_float(mesh, "BondCutoff", assembled["cutoffs"])

    # ---- core attributes -------------------------------------------------
    numbers, colors, cov_radii, vdw_radii, indices = [], [], [], [], []
    elements_in_file = sorted(set(symbols))
    for i, sym in enumerate(symbols):
        entry = entry_for(sym)
        numbers.append(entry["number"])
        colors.append(hex_to_linear_rgb(entry["color"])[0:3])
        cov_radii.append(display_radii[i])
        vdw_radii.append(entry["vdw_radius"] * 0.01)
        indices.append(elements_in_file.index(sym))

    _write_point_attribute(mesh, "AtomicNumber", 'INT', numbers)
    _write_point_attribute(mesh, "ElementColor", 'FLOAT_VECTOR', colors)
    _write_point_attribute(mesh, "AtomRadius", 'FLOAT', cov_radii)
    _write_point_attribute(mesh, "VdWRadius", 'FLOAT', vdw_radii)
    _write_point_attribute(mesh, "IndexOfElement", 'INT', indices)
    _write_point_attribute(
        mesh, "IsGhost", 'INT',
        [0] * n_real + [1] * (n_total - n_real),
    )
    _write_point_attribute(mesh, "ElementHidden", 'INT', [0] * n_total)

    # ---- per-point random value (for material variation) ----------------
    rng = random.Random(settings.get("seed", 0))
    _write_point_attribute(
        mesh, "RandomValue", 'FLOAT',
        [rng.random() for _ in range(n_total)],
    )

    # ---- average bond direction (for oriented instances) ----------------
    _write_point_attribute(
        mesh, "BondDir", 'FLOAT_VECTOR',
        bonds_mod.compute_bond_directions(verts, edges),
    )

    # ---- extra columns from extxyz / csv --------------------------------
    for attr_name, spec in structure.get("extra", {}).items():
        values = list(spec["values"])
        if len(values) < n_total:  # extend to ghosts by copying sources
            for s in ghost_sources:
                values.append(spec["values"][s])
        _write_point_attribute(mesh, attr_name, spec["type"], values)

    mesh.update()

    info = {
        "offset": offset,
        "symbols": symbols,
        "elements_in_file": elements_in_file,
        "n_real": n_real,
        "n_ghosts": n_total - n_real,
        "n_bonds": len(edges),
        "extra_attributes": sorted(structure.get("extra", {})),
        "n_frames": len(frames),
        "lattices": lattices,
    }
    if len(frames) > 1:
        info["frame_verts"] = [shifted(p) for p in assembled["positions"]]
    return mesh, info


def apply_hidden_elements(obj):
    """Write the ElementHidden attribute from the object's hidden list."""
    hidden = get_hidden_elements(obj)
    elements_in_file = get_molecule_elements(obj)
    hidden_indices = {
        i for i, sym in enumerate(elements_in_file) if sym in hidden
    }
    mesh = obj.data
    index_attr = mesh.attributes.get("IndexOfElement")
    if index_attr is None:
        return
    n = len(index_attr.data)
    indices = [0] * n
    index_attr.data.foreach_get("value", indices)
    values = [1 if idx in hidden_indices else 0 for idx in indices]
    _write_point_attribute(mesh, "ElementHidden", 'INT', values)
    mesh.update()


# ---------------------------------------------------------------------------
# Trajectory animation (shape keys)
# ---------------------------------------------------------------------------

def _flat(verts):
    return [c for v in verts for c in v]


def _insert_linear_keyframes(context, id_block, data_path, frame_values,
                             index=-1):
    """Keyframe ``id_block.<data_path>`` at each ``(frame, value)`` with
    linear interpolation. Linear is set through the new-keyframe user
    preference (restored afterwards), which works the same with the old
    and the layered (4.4+ / 5.x) action API."""
    edit = getattr(getattr(context, "preferences", None), "edit", None)
    old = getattr(edit, "keyframe_new_interpolation_type", None)
    try:
        if old is not None:
            edit.keyframe_new_interpolation_type = 'LINEAR'
        for frame, value in frame_values:
            setattr(id_block, data_path, value)
            id_block.keyframe_insert(data_path, frame=frame, index=index)
    finally:
        if old is not None:
            edit.keyframe_new_interpolation_type = old


def add_frame_shape_keys(context, obj, frame_verts, frame_start):
    """One shape key per trajectory frame; frame k of the trajectory
    lands on scene frame ``frame_start + k``.

    Uses *absolute* shape keys driven by the key's Evaluation Time, which
    Blender interpolates between consecutive keys. If this Blender does
    not space the keys out in time as expected, falls back to relative
    keys with one-hot keyframed values (same result, more F-curves).
    """
    n = len(frame_verts)
    if n < 2:
        return None
    obj.shape_key_add(name="Frame 1", from_mix=False)   # basis = frame 1
    for k in range(1, n):
        block = obj.shape_key_add(name=f"Frame {k + 1}", from_mix=False)
        block.data.foreach_set("co", _flat(frame_verts[k]))
    key = obj.data.shape_keys
    for block in key.key_blocks:
        block.interpolation = 'KEY_LINEAR'

    key.use_relative = False
    times = [block.frame for block in key.key_blocks]
    increasing = all(t1 > t0 for t0, t1 in zip(times, times[1:]))
    if increasing:
        _insert_linear_keyframes(
            context, key, "eval_time",
            [(frame_start + k, times[k]) for k in range(n)],
        )
        return 'ABSOLUTE'

    key.use_relative = True
    blocks = list(key.key_blocks)
    for k in range(1, n):
        values = [(frame_start + k, 1.0)]
        values.append((frame_start + k - 1, 0.0))
        if k + 1 < n:
            values.append((frame_start + k + 1, 0.0))
        _insert_linear_keyframes(context, blocks[k], "value", values)
    return 'RELATIVE'


# ---------------------------------------------------------------------------
# Full build
# ---------------------------------------------------------------------------

def build_molecule(context, structure, elements, mol_name, settings,
                   existing_obj=None, frames=None):
    """Create (or rebuild, when ``existing_obj`` is given) the molecule.

    ``frames``: optional list of per-frame structures (the first equals
    ``structure``); two or more give an animated trajectory.

    Returns ``(obj, info)``.

    Raises ValueError if a supercell is requested for a structure that
    has no lattice, or if it would be too large.
    """
    frames = list(frames) if frames else [structure]
    repeats = normalize_repeats(settings.get("supercell", (1, 1, 1)))
    if repeats != (1, 1, 1):
        if len(frames) * len(frames[0]["positions"]) * \
                repeats[0] * repeats[1] * repeats[2] > \
                MAX_SUPERCELL_ATOMS * 20:
            raise ValueError(
                "This supercell of every trajectory frame would be too "
                "large; use smaller repeats or import a single frame."
            )
        frames = [expand_supercell(f, repeats) for f in frames]
    structure = frames[0]

    mesh, info = build_molecule_mesh(structure, elements, mol_name,
                                     settings, frames=frames)
    info["repeats"] = repeats
    world_scale = settings.get("world_scale", 0.01)

    if existing_obj is None:
        obj = bpy.data.objects.new(f"{mol_name}_obj", mesh)
        obj.scale = (world_scale, world_scale, world_scale)

        collection = bpy.data.collections.new(f"{mol_name}_collection")
        context.scene.collection.children.link(collection)
        collection.objects.link(obj)

        atom_mat = materials.create_atom_material(mol_name)
        bond_mat = materials.create_bond_material(mol_name)
        atom_col = make_atom_collection(
            context, mol_name, info["elements_in_file"], elements, atom_mat
        )
        nodes.create_molecule_geometry_nodes(
            obj, atom_mat, bond_mat, atom_col,
            group_name=f"{mol_name}_GN",
            atom_scale=settings.get("atom_scale", 0.7),
            bond_radius=settings.get("bond_radius", 0.1),
        )
    else:
        obj = existing_obj
        old_mesh = obj.data
        obj.data = mesh
        if old_mesh.users == 0:
            bpy.data.meshes.remove(old_mesh)

        # Rebuild templates only when the element set changed, so instance
        # swaps survive a reload.
        old_elements = get_molecule_elements(obj)
        atom_col = find_atom_collection(obj)
        if atom_col is None or old_elements != info["elements_in_file"]:
            atom_mat = materials.create_atom_material(mol_name)
            new_col = make_atom_collection(
                context, mol_name, info["elements_in_file"], elements,
                atom_mat,
            )
            mod = nodes.find_molecule_modifier(obj)
            if mod is not None:
                for node in mod.node_group.nodes:
                    if node.bl_idname == "GeometryNodeCollectionInfo":
                        node.inputs["Collection"].default_value = new_col
                        break   # don't keep iterating after writing
            if atom_col is not None:
                for template in list(atom_col.objects):
                    data = template.data
                    bpy.data.objects.remove(template)
                    if data is not None and data.users == 0:
                        if isinstance(data, bpy.types.Mesh):
                            bpy.data.meshes.remove(data)
                bpy.data.collections.remove(atom_col)

    # ---- object bookkeeping ---------------------------------------------
    obj[PROP_FLAG] = 1
    obj[PROP_SOURCE] = settings.get("source", "")
    obj[PROP_KIND] = settings.get("kind", "XYZ")
    obj[PROP_SETTINGS] = json.dumps(settings)
    obj[PROP_ELEMENTS] = json.dumps(info["elements_in_file"])
    obj[PROP_OFFSET] = list(info["offset"])
    obj[PROP_LATTICE] = json.dumps(structure.get("lattice") or "")
    obj[PROP_REPEATS] = json.dumps(list(repeats))
    obj[PROP_NAME] = mol_name
    counts = {}
    for sym in info["symbols"][:info["n_real"]]:
        counts[sym] = counts.get(sym, 0) + 1
    obj[PROP_COUNTS] = json.dumps(counts)
    obj[PROP_TABLE] = json.dumps({
        sym: elements.get(sym) or elements["X"]
        for sym in info["elements_in_file"]
    })
    if PROP_HIDDEN not in obj:
        obj[PROP_HIDDEN] = json.dumps([])

    apply_hidden_elements(obj)

    # ---- trajectory animation -------------------------------------------
    n_frames = info["n_frames"]
    obj[PROP_FRAMES] = n_frames
    mod = nodes.find_molecule_modifier(obj)
    if mod is not None:
        nodes.set_gn_input(mod, "Dynamic Bonds", n_frames > 1)
    frame_start = context.scene.frame_start
    if n_frames > 1:
        info["shape_key_mode"] = add_frame_shape_keys(
            context, obj, info["frame_verts"], frame_start
        )
        last = frame_start + n_frames - 1
        if context.scene.frame_end < last:
            context.scene.frame_end = last
            info["frame_end_extended"] = last

    # ---- unit cell -------------------------------------------------------
    old_cell = next(
        (child for child in obj.children
         if child.name.startswith(f"{mol_name}_cell")), None,
    )
    if old_cell is not None:
        cell_mesh = old_cell.data
        bpy.data.objects.remove(old_cell)
        if cell_mesh is not None and cell_mesh.users == 0:
            bpy.data.meshes.remove(cell_mesh)
    if structure.get("lattice") is not None and settings.get(
            "show_unit_cell", True):
        parent_col = (obj.users_collection[0]
                      if obj.users_collection else context.scene.collection)
        cell = make_unit_cell_object(
            context, mol_name, structure["lattice"], info["offset"],
            obj, parent_col, world_scale, repeats=repeats,
        )
        lattices = info["lattices"]
        if n_frames > 1 and trajectory.lattices_vary(lattices):
            # cell optimisation / NPT: animate the cell box too
            cell_verts = [
                unit_cell_geometry(lat, info["offset"], repeats)[0]
                for lat in lattices
            ]
            add_frame_shape_keys(context, cell, cell_verts, frame_start)

    context.view_layer.objects.active = obj
    obj.select_set(True)
    return obj, info
