"""Geometry Nodes tree for the molecule object.

One node group drives everything the sidebar exposes:

* Style switch: 0 Ball-and-Stick, 1 Space-Filling, 2 Licorice, 3 Wireframe.
* Per-element / per-atom hiding via the ``ElementHidden`` point attribute.
* Ghost-atom visibility (periodic images) via the ``IsGhost`` attribute.
* Per-instance variation: random rotation, scale jitter, seed.
* Optional alignment of instances to a per-point vector attribute
  (e.g. ``forces`` from extxyz, or the computed ``BondDir``).

Only nodes available since Blender 4.0 are used, so the tree builds on any
4.x/5.x Blender.
"""

import json

import bpy

# Keys under which the build stores a name -> identifier map and a
# name -> default map on the node group, so lookups never have to walk
# the live interface collection.
SOCKET_ID_KEY = "xyz_socket_ids"
SOCKET_DEFAULT_KEY = "xyz_socket_defaults"

# Style socket values (imported by the preset operator as well).
STYLE_BALL_AND_STICK = 0
STYLE_SPACE_FILLING = 1
STYLE_LICORICE = 2
STYLE_WIREFRAME = 3

STYLE_ITEMS = [
    ('BALL_AND_STICK', "Ball-and-Stick",
     "Atoms at covalent radius x atom scale, with bond cylinders"),
    ('SPACE_FILLING', "Space-Filling",
     "Atoms at their van der Waals radius, no bonds (CPK view)"),
    ('LICORICE', "Licorice / Sticks",
     "Uniform thin atoms and bonds of equal radius"),
    ('WIREFRAME', "Wireframe",
     "Thin bond lines only, no atom spheres"),
]

STYLE_TO_INT = {
    'BALL_AND_STICK': STYLE_BALL_AND_STICK,
    'SPACE_FILLING': STYLE_SPACE_FILLING,
    'LICORICE': STYLE_LICORICE,
    'WIREFRAME': STYLE_WIREFRAME,
}


def _socket(node, name, sock_types=None, output=False):
    """Find a (possibly type-multiplexed) node socket by name.

    Nodes like Compare, Switch, and Random Value expose one socket per data
    type under the same name; the enabled one depends on ``data_type``.
    """
    sockets = node.outputs if output else node.inputs
    matches = [s for s in sockets if s.name == name]
    if not matches:
        raise KeyError(f"{node.bl_idname} has no socket named '{name}'")
    if len(matches) == 1:
        return matches[0]
    if sock_types is not None:
        for s in matches:
            if s.type in sock_types and s.enabled:
                return s
        for s in matches:
            if s.type in sock_types:
                return s
    for s in matches:
        if s.enabled:
            return s
    return matches[0]


def _new_node(ng, *bl_idnames):
    """Create the first node type that exists in this Blender version."""
    last_error = None
    for idname in bl_idnames:
        try:
            return ng.nodes.new(idname)
        except RuntimeError as exc:
            last_error = exc
    raise last_error


def iter_input_sockets(ng):
    """Snapshot of a group's input sockets as plain tuples.

    Returning a fully materialised list matters for safety, not style:
    ``interface.items_tree`` is a live RNA collection, and touching the
    node group or its modifier while an iterator over it is still open
    can invalidate that iterator and crash Blender with a use-after-free.
    Callers therefore always work from this snapshot, never from the
    collection itself.
    """
    return [
        (item.identifier, item.name, item.socket_type,
         getattr(item, "default_value", None))
        for item in ng.interface.items_tree
        if item.item_type == 'SOCKET' and item.in_out == 'INPUT'
    ]


def _write_modifier_input(mod, identifier, value):
    """Write one value onto a Geometry Nodes modifier (4.x/5.x safe)."""
    try:
        mod[identifier] = value
        return True
    except (KeyError, TypeError):
        pass
    try:  # newer property-bag API
        socket = getattr(mod.properties.inputs, identifier)
        socket.value = value
        return True
    except AttributeError:
        return False


def set_gn_input(mod, socket_name, value):
    """Set a Geometry Nodes modifier input by socket name."""
    identifier = gn_input_identifier(mod, socket_name)
    if identifier is None:
        return False
    return _write_modifier_input(mod, identifier, value)


def push_interface_defaults(mod):
    """Write every non-geometry group input default onto the modifier.

    Blender seeds a modifier's input values from the node group's
    interface at the moment the group is assigned. Any socket added
    afterwards has no value on the modifier, and an unset boolean reads
    as False -- which would hide the atoms and bonds entirely. Pushing
    the defaults explicitly makes the result independent of ordering.
    """
    ng = mod.node_group
    if ng is None:
        return
    ids = stored_socket_map(ng, SOCKET_ID_KEY)
    defaults = stored_socket_map(ng, SOCKET_DEFAULT_KEY)
    if ids and defaults:
        for name, default in defaults.items():
            identifier = ids.get(name)
            if identifier is not None:
                _write_modifier_input(mod, identifier, default)
        return
    for identifier, _name, socket_type, default in iter_input_sockets(ng):
        if socket_type != 'NodeSocketGeometry' and default is not None:
            _write_modifier_input(mod, identifier, default)


def stored_socket_map(ng, key):
    """Name -> value map stored on the node group at build time, or {}."""
    try:
        raw = ng.get(key)
        return json.loads(raw) if raw else {}
    except (TypeError, ValueError):
        return {}


def gn_input_identifier(mod, socket_name):
    """The modifier property key ('Socket_2') for a named group input.

    Uses the map recorded when the tree was built; falls back to
    scanning the interface for node groups made by older versions.
    """
    ng = mod.node_group
    if ng is None:
        return None
    identifier = stored_socket_map(ng, SOCKET_ID_KEY).get(socket_name)
    if identifier:
        return identifier
    for identifier, name, _socket_type, _default in iter_input_sockets(ng):
        if name == socket_name:
            return identifier
    return None


def get_gn_input(mod, socket_name):
    identifier = gn_input_identifier(mod, socket_name)
    if identifier is None:
        return None
    try:
        return mod[identifier]
    except (KeyError, TypeError):
        return None


def find_molecule_modifier(obj):
    if obj is None:
        return None
    for mod in obj.modifiers:
        if mod.type == 'NODES' and mod.node_group is not None \
                and mod.node_group.get("xyz_molecule_tree"):
            return mod
    return None


def _new_float_input(iface, name, default, minimum=None, maximum=None,
                     description=""):
    sock = iface.new_socket(name, in_out='INPUT',
                            socket_type='NodeSocketFloat',
                            description=description)
    sock.default_value = default
    if minimum is not None:
        sock.min_value = minimum
    if maximum is not None:
        sock.max_value = maximum
    return sock


def create_molecule_geometry_nodes(obj, atom_material, bond_material,
                                   atom_collection, group_name,
                                   atom_scale=0.7, bond_radius=0.1):
    """Build the full atoms + bonds Geometry Nodes tree on ``obj``."""
    assert obj.type == 'MESH'

    ng = bpy.data.node_groups.new(group_name, 'GeometryNodeTree')
    ng["xyz_molecule_tree"] = True

    iface = ng.interface
    socket_ids = {}
    socket_defaults = {}

    def record(sock, default):
        socket_ids[sock.name] = sock.identifier
        socket_defaults[sock.name] = default
        return sock

    iface.new_socket("Geometry", in_out='INPUT',
                     socket_type='NodeSocketGeometry')

    sock = iface.new_socket("Style", in_out='INPUT',
                            socket_type='NodeSocketInt',
                            description=("0 Ball-and-Stick, 1 Space-"
                                         "Filling, 2 Licorice, 3 Wireframe"))
    sock.default_value = STYLE_BALL_AND_STICK
    sock.min_value = 0
    sock.max_value = 3
    record(sock, STYLE_BALL_AND_STICK)

    record(_new_float_input(iface, "Atom Scale", atom_scale, 0.0,
                            description="Sphere radius as a fraction of "
                                        "the element radius"),
           atom_scale)
    record(_new_float_input(iface, "Bond Radius", bond_radius, 0.0,
                            description="Radius of the bond cylinders "
                                        "(Angstrom)"),
           bond_radius)

    for name, default, desc in (
        ("Show Atoms", True, "Display atom spheres"),
        ("Show Bonds", True, "Display bonds"),
        ("Show Ghost Atoms", True,
         "Display periodic-image atoms at cell boundaries"),
    ):
        sock = iface.new_socket(name, in_out='INPUT',
                                socket_type='NodeSocketBool',
                                description=desc)
        sock.default_value = default
        record(sock, default)

    record(_new_float_input(iface, "Random Rotation", 0.0, 0.0, 1.0,
                            description="Random per-instance rotation "
                                        "(0 = none, 1 = fully random)"),
           0.0)
    record(_new_float_input(iface, "Scale Jitter", 0.0, 0.0, 1.0,
                            description="Random per-instance scale "
                                        "variation"),
           0.0)
    sock = iface.new_socket("Seed", in_out='INPUT',
                            socket_type='NodeSocketInt',
                            description="Random seed for rotation/scale "
                                        "variation")
    sock.default_value = 0
    record(sock, 0)

    sock = iface.new_socket("Align to Vector", in_out='INPUT',
                            socket_type='NodeSocketBool',
                            description="Align instance +Z to the vector "
                                        "attribute below")
    sock.default_value = False
    record(sock, False)
    sock = iface.new_socket("Vector Attribute", in_out='INPUT',
                            socket_type='NodeSocketString',
                            description="Per-point vector attribute used "
                                        "for alignment (e.g. BondDir, "
                                        "forces)")
    sock.default_value = "BondDir"
    record(sock, "BondDir")

    sock = iface.new_socket("Dynamic Bonds", in_out='INPUT',
                            socket_type='NodeSocketBool',
                            description="Hide bonds longer than their "
                                        "per-edge 'BondCutoff' at the "
                                        "current frame (animated "
                                        "trajectories)")
    sock.default_value = False
    record(sock, False)

    iface.new_socket("Geometry", in_out='OUTPUT',
                     socket_type='NodeSocketGeometry')

    ng[SOCKET_ID_KEY] = json.dumps(socket_ids)
    ng[SOCKET_DEFAULT_KEY] = json.dumps(socket_defaults)

    nodes = ng.nodes
    links = ng.links
    nodes.clear()

    group_in = nodes.new("NodeGroupInput")
    group_in.location = (-2200, 0)
    group_out = nodes.new("NodeGroupOutput")
    group_out.location = (1700, 0)

    def gi(name):
        return group_in.outputs[name]

    # ---- dynamic bonds (trajectories) -----------------------------------
    # A trajectory's mesh carries the union of the bonds found in any
    # frame. With "Dynamic Bonds" on, every edge whose current length
    # exceeds its own 'BondCutoff' (r_i + r_j) x scale is dropped, so
    # bonds form and break at the right frame. Vertices are kept.
    edge_verts = nodes.new("GeometryNodeInputMeshEdgeVertices")
    edge_verts.location = (-2000, 1100)

    edge_length = nodes.new("ShaderNodeVectorMath")
    edge_length.operation = 'DISTANCE'
    edge_length.location = (-1800, 1100)
    links.new(edge_verts.outputs["Position 1"], edge_length.inputs[0])
    links.new(edge_verts.outputs["Position 2"], edge_length.inputs[1])

    attr_cutoff = nodes.new("GeometryNodeInputNamedAttribute")
    attr_cutoff.location = (-1800, 940)
    attr_cutoff.data_type = 'FLOAT'
    attr_cutoff.inputs["Name"].default_value = "BondCutoff"

    too_long = nodes.new("FunctionNodeCompare")
    too_long.data_type = 'FLOAT'
    too_long.operation = 'GREATER_THAN'
    too_long.location = (-1600, 1040)
    links.new(_socket(edge_length, "Value", output=True),
              _socket(too_long, "A", {'VALUE'}))
    links.new(attr_cutoff.outputs["Attribute"],
              _socket(too_long, "B", {'VALUE'}))

    drop_edge = nodes.new("FunctionNodeBooleanMath")
    drop_edge.operation = 'AND'
    drop_edge.location = (-1420, 1040)
    links.new(gi("Dynamic Bonds"), drop_edge.inputs[0])
    links.new(too_long.outputs["Result"], drop_edge.inputs[1])

    del_long = nodes.new("GeometryNodeDeleteGeometry")
    del_long.location = (-1240, 1040)
    del_long.domain = 'EDGE'
    del_long.mode = 'EDGE_FACE'          # keep the vertices (atoms)
    links.new(gi("Geometry"), del_long.inputs["Geometry"])
    links.new(drop_edge.outputs["Boolean"], del_long.inputs["Selection"])
    source_geometry = del_long.outputs["Geometry"]

    # A ghost atom is only there to terminate a periodic bond: once that
    # bond is gone (dynamic mode), hide the ghost too.
    neighbors = nodes.new("GeometryNodeInputMeshVertexNeighbors")
    neighbors.location = (-1600, 860)

    no_bonds = nodes.new("FunctionNodeCompare")
    no_bonds.data_type = 'INT'
    no_bonds.operation = 'EQUAL'
    no_bonds.location = (-1420, 860)
    links.new(neighbors.outputs["Vertex Count"],
              _socket(no_bonds, "A", {'INT'}))
    _socket(no_bonds, "B", {'INT'}).default_value = 0

    # ---- shared "hidden" field ------------------------------------------
    # hidden = max(ElementHidden, IsGhost * (1 - ShowGhostAtoms))
    attr_hidden = nodes.new("GeometryNodeInputNamedAttribute")
    attr_hidden.location = (-2000, 700)
    attr_hidden.data_type = 'FLOAT'
    attr_hidden.inputs["Name"].default_value = "ElementHidden"

    attr_ghost = nodes.new("GeometryNodeInputNamedAttribute")
    attr_ghost.location = (-2000, 520)
    attr_ghost.data_type = 'FLOAT'
    attr_ghost.inputs["Name"].default_value = "IsGhost"

    ghosts_off = nodes.new("ShaderNodeMath")
    ghosts_off.operation = 'SUBTRACT'
    ghosts_off.location = (-1800, 560)
    ghosts_off.inputs[0].default_value = 1.0
    links.new(gi("Show Ghost Atoms"), ghosts_off.inputs[1])

    ghost_hidden = nodes.new("ShaderNodeMath")
    ghost_hidden.operation = 'MULTIPLY'
    ghost_hidden.location = (-1620, 560)
    links.new(attr_ghost.outputs["Attribute"], ghost_hidden.inputs[0])
    links.new(ghosts_off.outputs["Value"], ghost_hidden.inputs[1])

    hidden = nodes.new("ShaderNodeMath")
    hidden.operation = 'MAXIMUM'
    hidden.location = (-1440, 640)
    links.new(attr_hidden.outputs["Attribute"], hidden.inputs[0])
    links.new(ghost_hidden.outputs["Value"], hidden.inputs[1])

    hide_point = nodes.new("ShaderNodeMath")   # point domain: > 0.5
    hide_point.operation = 'GREATER_THAN'
    hide_point.location = (-1260, 700)
    hide_point.inputs[1].default_value = 0.5
    links.new(hidden.outputs["Value"], hide_point.inputs[0])

    hide_edge = nodes.new("ShaderNodeMath")    # edge domain: avg > 0.25
    hide_edge.operation = 'GREATER_THAN'
    hide_edge.location = (-1260, 520)
    hide_edge.inputs[1].default_value = 0.25
    links.new(hidden.outputs["Value"], hide_edge.inputs[0])

    # ---- style comparisons ----------------------------------------------
    def style_equals(value, location):
        cmp = nodes.new("FunctionNodeCompare")
        cmp.data_type = 'INT'
        cmp.operation = 'EQUAL'
        cmp.location = location
        links.new(gi("Style"), _socket(cmp, "A", {'INT'}))
        _socket(cmp, "B", {'INT'}).default_value = value
        return cmp

    is_sf = style_equals(STYLE_SPACE_FILLING, (-2000, 200))
    is_lic = style_equals(STYLE_LICORICE, (-2000, 40))
    is_wire = style_equals(STYLE_WIREFRAME, (-2000, -120))

    # =====================================================================
    # Atom branch
    # =====================================================================
    # orphan ghost = Dynamic Bonds AND no remaining bond AND IsGhost
    orphan_a = nodes.new("FunctionNodeBooleanMath")
    orphan_a.operation = 'AND'
    orphan_a.location = (-1240, 860)
    links.new(gi("Dynamic Bonds"), orphan_a.inputs[0])
    links.new(no_bonds.outputs["Result"], orphan_a.inputs[1])

    orphan_ghost = nodes.new("FunctionNodeBooleanMath")
    orphan_ghost.operation = 'AND'
    orphan_ghost.location = (-1240, 760)
    links.new(orphan_a.outputs["Boolean"], orphan_ghost.inputs[0])
    links.new(attr_ghost.outputs["Attribute"], orphan_ghost.inputs[1])

    point_hidden = nodes.new("FunctionNodeBooleanMath")
    point_hidden.operation = 'OR'
    point_hidden.location = (-1100, 700)
    links.new(hide_point.outputs["Value"], point_hidden.inputs[0])
    links.new(orphan_ghost.outputs["Boolean"], point_hidden.inputs[1])

    del_hidden_pts = nodes.new("GeometryNodeDeleteGeometry")
    del_hidden_pts.location = (-1050, 300)
    del_hidden_pts.domain = 'POINT'
    links.new(source_geometry, del_hidden_pts.inputs["Geometry"])
    links.new(point_hidden.outputs["Boolean"],
              del_hidden_pts.inputs["Selection"])

    mesh_to_points = nodes.new("GeometryNodeMeshToPoints")
    mesh_to_points.location = (-850, 300)
    links.new(del_hidden_pts.outputs["Geometry"],
              mesh_to_points.inputs["Mesh"])

    coll_info = nodes.new("GeometryNodeCollectionInfo")
    coll_info.location = (-650, 480)
    coll_info.inputs["Collection"].default_value = atom_collection
    coll_info.inputs["Separate Children"].default_value = True
    coll_info.inputs["Reset Children"].default_value = True

    elem_index = nodes.new("GeometryNodeInputNamedAttribute")
    elem_index.location = (-850, 120)
    elem_index.data_type = 'INT'
    elem_index.inputs["Name"].default_value = "IndexOfElement"

    # ---- per-point radius, per style ------------------------------------
    attr_radius = nodes.new("GeometryNodeInputNamedAttribute")
    attr_radius.location = (-1700, -50)
    attr_radius.data_type = 'FLOAT'
    attr_radius.inputs["Name"].default_value = "AtomRadius"

    attr_vdw = nodes.new("GeometryNodeInputNamedAttribute")
    attr_vdw.location = (-1700, -220)
    attr_vdw.data_type = 'FLOAT'
    attr_vdw.inputs["Name"].default_value = "VdWRadius"

    pick_radius = nodes.new("GeometryNodeSwitch")   # covalent vs vdW
    pick_radius.input_type = 'FLOAT'
    pick_radius.location = (-1500, -120)
    links.new(is_sf.outputs["Result"], pick_radius.inputs["Switch"])
    links.new(attr_radius.outputs["Attribute"],
              _socket(pick_radius, "False", {'VALUE'}))
    links.new(attr_vdw.outputs["Attribute"],
              _socket(pick_radius, "True", {'VALUE'}))

    scaled_radius = nodes.new("ShaderNodeMath")
    scaled_radius.operation = 'MULTIPLY'
    scaled_radius.location = (-1320, -120)
    links.new(_socket(pick_radius, "Output", output=True),
              scaled_radius.inputs[0])
    links.new(gi("Atom Scale"), scaled_radius.inputs[1])

    licorice_radius = nodes.new("GeometryNodeSwitch")  # licorice: bond size
    licorice_radius.input_type = 'FLOAT'
    licorice_radius.location = (-1140, -120)
    links.new(is_lic.outputs["Result"], licorice_radius.inputs["Switch"])
    links.new(scaled_radius.outputs["Value"],
              _socket(licorice_radius, "False", {'VALUE'}))
    links.new(gi("Bond Radius"),
              _socket(licorice_radius, "True", {'VALUE'}))

    # ---- alignment rotation ---------------------------------------------
    attr_vector = nodes.new("GeometryNodeInputNamedAttribute")
    attr_vector.location = (-850, -300)
    attr_vector.data_type = 'FLOAT_VECTOR'
    links.new(gi("Vector Attribute"), attr_vector.inputs["Name"])

    align = _new_node(ng,
                      "FunctionNodeAlignRotationToVector",   # 4.1+
                      "FunctionNodeAlignEulerToVector")      # 4.0 fallback
    align.location = (-650, -300)
    align.axis = 'Z'
    links.new(gi("Align to Vector"), _socket(align, "Factor"))
    links.new(attr_vector.outputs["Attribute"], _socket(align, "Vector"))

    instance = nodes.new("GeometryNodeInstanceOnPoints")
    instance.location = (-380, 300)
    instance.inputs["Pick Instance"].default_value = True
    links.new(mesh_to_points.outputs["Points"], instance.inputs["Points"])
    links.new(coll_info.outputs["Instances"], instance.inputs["Instance"])
    links.new(elem_index.outputs["Attribute"],
              instance.inputs["Instance Index"])
    links.new(_socket(licorice_radius, "Output", output=True),
              instance.inputs["Scale"])
    links.new(align.outputs["Rotation"], instance.inputs["Rotation"])

    # ---- per-instance variation -----------------------------------------
    rot_amount = nodes.new("ShaderNodeVectorMath")
    rot_amount.operation = 'SCALE'
    rot_amount.location = (-560, 0)
    rot_amount.inputs[0].default_value = (3.14159, 3.14159, 3.14159)
    links.new(gi("Random Rotation"), _socket(rot_amount, "Scale"))

    rot_amount_neg = nodes.new("ShaderNodeVectorMath")
    rot_amount_neg.operation = 'SCALE'
    rot_amount_neg.location = (-560, -160)
    links.new(rot_amount.outputs["Vector"], rot_amount_neg.inputs[0])
    _socket(rot_amount_neg, "Scale").default_value = -1.0

    rand_rot = nodes.new("FunctionNodeRandomValue")
    rand_rot.data_type = 'FLOAT_VECTOR'
    rand_rot.location = (-380, 0)
    links.new(rot_amount_neg.outputs["Vector"],
              _socket(rand_rot, "Min", {'VECTOR'}))
    links.new(rot_amount.outputs["Vector"],
              _socket(rand_rot, "Max", {'VECTOR'}))
    links.new(gi("Seed"), _socket(rand_rot, "Seed"))

    rotate_instances = nodes.new("GeometryNodeRotateInstances")
    rotate_instances.location = (-160, 300)
    links.new(instance.outputs["Instances"],
              rotate_instances.inputs["Instances"])
    links.new(_socket(rand_rot, "Value", {'VECTOR'}, output=True),
              rotate_instances.inputs["Rotation"])

    seed_plus_one = nodes.new("ShaderNodeMath")
    seed_plus_one.operation = 'ADD'
    seed_plus_one.location = (-560, -330)
    seed_plus_one.inputs[1].default_value = 1.0
    links.new(gi("Seed"), seed_plus_one.inputs[0])

    rand_scale = nodes.new("FunctionNodeRandomValue")
    rand_scale.data_type = 'FLOAT'
    rand_scale.location = (-380, -330)
    _socket(rand_scale, "Min", {'VALUE'}).default_value = -1.0
    _socket(rand_scale, "Max", {'VALUE'}).default_value = 1.0
    links.new(seed_plus_one.outputs["Value"], _socket(rand_scale, "Seed"))

    jitter = nodes.new("ShaderNodeMath")
    jitter.operation = 'MULTIPLY'
    jitter.location = (-200, -330)
    links.new(_socket(rand_scale, "Value", {'VALUE'}, output=True),
              jitter.inputs[0])
    links.new(gi("Scale Jitter"), jitter.inputs[1])

    jitter_scale = nodes.new("ShaderNodeMath")
    jitter_scale.operation = 'ADD'
    jitter_scale.location = (-20, -330)
    jitter_scale.inputs[1].default_value = 1.0
    links.new(jitter.outputs["Value"], jitter_scale.inputs[0])

    scale_instances = nodes.new("GeometryNodeScaleInstances")
    scale_instances.location = (60, 300)
    links.new(rotate_instances.outputs["Instances"],
              scale_instances.inputs["Instances"])
    links.new(jitter_scale.outputs["Value"],
              scale_instances.inputs["Scale"])

    # NOTE: no Set Material here on purpose. The sphere templates carry
    # the per-element atom material on their mesh, and swapped-in objects
    # (cats, trees, ...) keep their own materials this way.

    # atoms visible = Show Atoms AND NOT wireframe
    not_wire = nodes.new("FunctionNodeBooleanMath")
    not_wire.operation = 'NOT'
    not_wire.location = (280, 120)
    links.new(is_wire.outputs["Result"], not_wire.inputs[0])

    atoms_visible = nodes.new("FunctionNodeBooleanMath")
    atoms_visible.operation = 'AND'
    atoms_visible.location = (460, 120)
    links.new(gi("Show Atoms"), atoms_visible.inputs[0])
    links.new(not_wire.outputs["Boolean"], atoms_visible.inputs[1])

    atom_switch = nodes.new("GeometryNodeSwitch")
    atom_switch.input_type = 'GEOMETRY'
    atom_switch.location = (640, 300)
    links.new(atoms_visible.outputs["Boolean"],
              atom_switch.inputs["Switch"])
    links.new(scale_instances.outputs["Instances"],
              _socket(atom_switch, "True", {'GEOMETRY'}))

    # =====================================================================
    # Bond branch
    # =====================================================================
    del_hidden_edges = nodes.new("GeometryNodeDeleteGeometry")
    del_hidden_edges.location = (-1050, -520)
    del_hidden_edges.domain = 'EDGE'
    links.new(source_geometry, del_hidden_edges.inputs["Geometry"])
    links.new(hide_edge.outputs["Value"],
              del_hidden_edges.inputs["Selection"])

    mesh_to_curve = nodes.new("GeometryNodeMeshToCurve")
    mesh_to_curve.location = (-850, -520)
    links.new(del_hidden_edges.outputs["Geometry"],
              mesh_to_curve.inputs["Mesh"])

    # wireframe bonds are thinner
    wire_radius = nodes.new("ShaderNodeMath")
    wire_radius.operation = 'MULTIPLY'
    wire_radius.location = (-850, -740)
    wire_radius.inputs[1].default_value = 0.35
    links.new(gi("Bond Radius"), wire_radius.inputs[0])

    bond_radius_final = nodes.new("GeometryNodeSwitch")
    bond_radius_final.input_type = 'FLOAT'
    bond_radius_final.location = (-670, -740)
    links.new(is_wire.outputs["Result"],
              bond_radius_final.inputs["Switch"])
    links.new(gi("Bond Radius"),
              _socket(bond_radius_final, "False", {'VALUE'}))
    links.new(wire_radius.outputs["Value"],
              _socket(bond_radius_final, "True", {'VALUE'}))

    circle = nodes.new("GeometryNodeCurvePrimitiveCircle")
    circle.location = (-480, -740)
    links.new(_socket(bond_radius_final, "Output", output=True),
              circle.inputs["Radius"])

    curve_to_mesh = nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.location = (-380, -520)
    curve_to_mesh.inputs["Fill Caps"].default_value = True
    links.new(mesh_to_curve.outputs["Curve"],
              curve_to_mesh.inputs["Curve"])
    links.new(circle.outputs["Curve"],
              curve_to_mesh.inputs["Profile Curve"])

    smooth = nodes.new("GeometryNodeSetShadeSmooth")
    smooth.location = (-160, -520)
    links.new(curve_to_mesh.outputs["Mesh"], smooth.inputs["Geometry"])

    set_bond_material = nodes.new("GeometryNodeSetMaterial")
    set_bond_material.location = (60, -520)
    set_bond_material.inputs["Material"].default_value = bond_material
    links.new(smooth.outputs["Geometry"],
              set_bond_material.inputs["Geometry"])

    # bonds visible = Show Bonds AND NOT space-filling
    not_sf = nodes.new("FunctionNodeBooleanMath")
    not_sf.operation = 'NOT'
    not_sf.location = (280, -700)
    links.new(is_sf.outputs["Result"], not_sf.inputs[0])

    bonds_visible = nodes.new("FunctionNodeBooleanMath")
    bonds_visible.operation = 'AND'
    bonds_visible.location = (460, -700)
    links.new(gi("Show Bonds"), bonds_visible.inputs[0])
    links.new(not_sf.outputs["Boolean"], bonds_visible.inputs[1])

    bond_switch = nodes.new("GeometryNodeSwitch")
    bond_switch.input_type = 'GEOMETRY'
    bond_switch.location = (640, -520)
    links.new(bonds_visible.outputs["Boolean"],
              bond_switch.inputs["Switch"])
    links.new(set_bond_material.outputs["Geometry"],
              _socket(bond_switch, "True", {'GEOMETRY'}))

    # =====================================================================
    # Join
    # =====================================================================
    join = nodes.new("GeometryNodeJoinGeometry")
    join.location = (1100, 0)
    links.new(_socket(atom_switch, "Output", output=True),
              join.inputs["Geometry"])
    links.new(_socket(bond_switch, "Output", output=True),
              join.inputs["Geometry"])
    links.new(join.outputs["Geometry"], group_out.inputs["Geometry"])

    # Assign the finished group only now, so the modifier is seeded from
    # a complete interface, then push every default explicitly.
    mod = obj.modifiers.new("Molecule", 'NODES')
    mod.node_group = ng
    push_interface_defaults(mod)
    return mod
