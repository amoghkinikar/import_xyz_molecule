"""Multi-frame (trajectory) assembly, pure Python.

A Blender mesh has one fixed set of vertices and edges, but a trajectory's
bonds (and periodic ghost atoms) can change from frame to frame. The
approach here:

* vertices = the real atoms, plus every ghost atom needed in *any* frame
  (a ghost is identified by its source atom and lattice image, so it can
  follow that atom through the trajectory);
* edges = the union of the bonds found in every frame, each with its own
  cutoff ``(r_i + r_j) x scale`` stored as the ``BondCutoff`` edge
  attribute;
* per-frame vertex positions become shape keys.

The Geometry Nodes tree then hides every edge that is longer than its
cutoff at the current frame (and any ghost left without a bond), which
shows bonds forming and breaking at the right time.

A single frame goes through the same code and gives exactly the static
result: that frame's bonds and ghosts, in the same order as before.
"""

from . import bonds as bonds_mod


# ---------------------------------------------------------------------------
# Small 3x3 helpers (no numpy / mathutils needed)
# ---------------------------------------------------------------------------

def _inverse_transpose_rows(lattice):
    """Inverse of the matrix whose *columns* are the lattice vectors, so
    that ``frac = inv @ cartesian``. Returns None for a singular cell."""
    (a0, a1, a2), (b0, b1, b2), (c0, c1, c2) = lattice
    # M = [[a0, b0, c0], [a1, b1, c1], [a2, b2, c2]]
    m = ((a0, b0, c0), (a1, b1, c1), (a2, b2, c2))
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    if abs(det) < 1e-12:
        return None
    inv = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r = [x for x in range(3) if x != j]
            c = [y for y in range(3) if y != i]
            minor = (m[r[0]][c[0]] * m[r[1]][c[1]]
                     - m[r[0]][c[1]] * m[r[1]][c[0]])
            inv[i][j] = ((-1) ** (i + j)) * minor / det
    return inv


def _lattice_of(frame, fallback):
    return frame.get("lattice") or fallback


# ---------------------------------------------------------------------------
# Unwrapping
# ---------------------------------------------------------------------------

def unwrap_positions(frames_positions, lattices, pbc=(True, True, True)):
    """Undo periodic wrapping between consecutive frames.

    MD codes often wrap atoms back into the cell; between two stored
    frames an atom can then jump by a whole lattice vector, which a
    linearly interpolated animation would show as the atom flying across
    the cell. Each frame's atoms are shifted by whole lattice vectors so
    that every displacement from the previous frame is the minimum image.

    ``frames_positions``: list (per frame) of lists of (x, y, z).
    ``lattices``: per-frame lattice (or None). Returns a new list.
    """
    if not frames_positions:
        return []
    result = [list(frames_positions[0])]
    for k in range(1, len(frames_positions)):
        lattice = lattices[k]
        prev = result[k - 1]
        cur = frames_positions[k]
        inv = _inverse_transpose_rows(lattice) if lattice else None
        if inv is None:
            result.append(list(cur))
            continue
        a, b, c = lattice
        new = []
        for (px, py, pz), (x, y, z) in zip(prev, cur):
            d = (x - px, y - py, z - pz)
            n = [0, 0, 0]
            for axis in range(3):
                if pbc[axis]:
                    f = (inv[axis][0] * d[0] + inv[axis][1] * d[1]
                         + inv[axis][2] * d[2])
                    n[axis] = round(f)
            new.append((
                x - n[0] * a[0] - n[1] * b[0] - n[2] * c[0],
                y - n[0] * a[1] - n[1] * b[1] - n[2] * c[1],
                z - n[0] * a[2] - n[1] * b[2] - n[2] * c[2],
            ))
        result.append(new)
    return result


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble(frames_positions, lattices, bond_radii, pbc=None, scale=1.1,
             detect_bonds=True, periodic_bonds=True):
    """Union of bonds and ghosts over all frames.

    Returns a dict with

    * ``positions``: per frame, the ``n_real + n_ghost`` vertex positions,
    * ``ghost_sources``: real atom index for each ghost,
    * ``edges``: ``[(i, j), ...]`` into those vertices,
    * ``cutoffs``: per edge, ``(r_i + r_j) * scale`` in Angstrom,
    * ``n_real``.
    """
    n = len(frames_positions[0]) if frames_positions else 0
    edges = []
    edge_set = set()

    def add_edge(i, j):
        key = (i, j) if i < j else (j, i)
        if key not in edge_set:
            edge_set.add(key)
            edges.append((i, j))

    ghost_index = {}        # (source, image) -> ghost number (0-based)
    ghost_keys = []
    ghost_edge_keys = []    # (real j, ghost key) in discovery order
    ghost_edge_set = set()

    if detect_bonds and n:
        for k, positions in enumerate(frames_positions):
            for i, j in bonds_mod.compute_bonds(positions, bond_radii,
                                                scale=scale):
                add_edge(i, j)
            lattice = lattices[k]
            if lattice is None or not periodic_bonds:
                continue
            _gpos, gsrc, gbonds, gshift = bonds_mod.compute_periodic_ghosts(
                positions, bond_radii, lattice, pbc=pbc, scale=scale,
                return_shifts=True,
            )
            local = []                  # local ghost index -> key
            for src, image in zip(gsrc, gshift):
                key = (src, image)
                if key not in ghost_index:
                    ghost_index[key] = len(ghost_keys)
                    ghost_keys.append(key)
                local.append(key)
            for j, g in gbonds:
                pair = (j, local[g - n])
                if pair not in ghost_edge_set:
                    ghost_edge_set.add(pair)
                    ghost_edge_keys.append(pair)

    for j, key in ghost_edge_keys:
        edges.append((j, n + ghost_index[key]))

    ghost_sources = [src for src, _image in ghost_keys]

    positions_out = []
    for k, positions in enumerate(frames_positions):
        frame = list(positions)
        lattice = lattices[k]
        for src, (ia, ib, ic) in ghost_keys:
            a, b, c = lattice
            x, y, z = positions[src]
            frame.append((
                x + ia * a[0] + ib * b[0] + ic * c[0],
                y + ia * a[1] + ib * b[1] + ic * c[1],
                z + ia * a[2] + ib * b[2] + ic * c[2],
            ))
        positions_out.append(frame)

    radii_all = list(bond_radii) + [bond_radii[s] for s in ghost_sources]
    cutoffs = [(radii_all[i] + radii_all[j]) * scale for i, j in edges]

    return {
        "positions": positions_out,
        "ghost_sources": ghost_sources,
        "edges": edges,
        "cutoffs": cutoffs,
        "n_real": n,
    }


def frame_lattices(frames):
    """Per-frame lattice, falling back to the first frame's lattice for
    frames whose comment line carries none."""
    first = frames[0].get("lattice") if frames else None
    return [_lattice_of(f, first) for f in frames]


def lattices_vary(lattices, tol=1e-6):
    base = lattices[0]
    if base is None:
        return False
    for lat in lattices[1:]:
        if lat is None:
            continue
        for u, v in zip(base, lat):
            if any(abs(p - q) > tol for p, q in zip(u, v)):
                return True
    return False
