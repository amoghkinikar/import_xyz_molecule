"""Bond detection: KD-tree distance search plus periodic (ghost-image)
bonding for structures that carry a lattice (extended XYZ).

Periodic bonding strategy
-------------------------
For each of the 26 neighboring lattice images, atoms are translated by the
image shift and kept as *ghost candidates* when they land within the bond
cutoff of the real structure's bounding box. Ghosts that bond to at least
one real atom are appended to the structure (tagged with the ``IsGhost``
attribute) and the bonds are added as ordinary edges. Because the search
runs symmetrically, both sides of a periodic bond receive their stub --
exactly what viewers like VESTA display at cell boundaries.
"""

try:
    from mathutils.kdtree import KDTree
except ImportError:          # outside Blender (tests)
    KDTree = None


def _build_kdtree(positions):
    if KDTree is None:
        return None
    kd = KDTree(len(positions))
    for i, co in enumerate(positions):
        kd.insert(co, i)
    kd.balance()
    return kd


def compute_bonds(positions, radii, scale=1.1):
    """Return ``[(i, j), ...]`` index pairs (i < j) for bonded atoms.

    A bond is made when the interatomic distance is smaller than the sum of
    the two covalent radii times ``scale``. Uses ``mathutils.kdtree`` inside
    Blender and a brute-force loop as a testing fallback.
    """
    n = len(positions)
    if n < 2:
        return []
    max_cutoff = 2.0 * max(radii) * scale

    edges = []
    kd = _build_kdtree(positions)
    if kd is not None:
        for i, co in enumerate(positions):
            for _co, j, dist in kd.find_range(co, max_cutoff):
                if j > i and dist < (radii[i] + radii[j]) * scale:
                    edges.append((i, j))
    else:
        for i in range(n):
            xi, yi, zi = positions[i]
            for j in range(i + 1, n):
                xj, yj, zj = positions[j]
                dist = ((xi - xj) ** 2 + (yi - yj) ** 2
                        + (zi - zj) ** 2) ** 0.5
                if dist < (radii[i] + radii[j]) * scale:
                    edges.append((i, j))
    return edges


def compute_periodic_ghosts(positions, radii, lattice, pbc=None, scale=1.1,
                            return_shifts=False):
    """Find ghost atoms and periodic bonds across cell boundaries.

    Returns ``(ghost_positions, ghost_sources, ghost_bonds)``:

    * ``ghost_positions`` -- coordinates of ghost atoms (to be appended
      after the ``len(positions)`` real atoms),
    * ``ghost_sources``   -- for each ghost, the index of the real atom it
      is an image of (used to copy attributes),
    * ``ghost_bonds``     -- ``(real_index, ghost_index)`` pairs where
      ``ghost_index`` already counts from ``len(positions)``.

    With ``return_shifts=True`` a fourth list is returned: the integer
    lattice image ``(ia, ib, ic)`` of each ghost, which identifies it
    across the frames of a trajectory.
    """
    n = len(positions)
    empty = ([], [], [], []) if return_shifts else ([], [], [])
    if n == 0 or lattice is None:
        return empty
    if pbc is None:
        pbc = (True, True, True)

    a, b, c = lattice
    cutoff = 2.0 * max(radii) * scale

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    lo = (min(xs) - cutoff, min(ys) - cutoff, min(zs) - cutoff)
    hi = (max(xs) + cutoff, max(ys) + cutoff, max(zs) + cutoff)

    shift_range = [
        (-1, 0, 1) if periodic else (0,) for periodic in pbc
    ]

    candidates = []      # (pos, source_index, image)
    for ia in shift_range[0]:
        for ib in shift_range[1]:
            for ic in shift_range[2]:
                if ia == 0 and ib == 0 and ic == 0:
                    continue
                sx = ia * a[0] + ib * b[0] + ic * c[0]
                sy = ia * a[1] + ib * b[1] + ic * c[1]
                sz = ia * a[2] + ib * b[2] + ic * c[2]
                for i, (x, y, z) in enumerate(positions):
                    gx, gy, gz = x + sx, y + sy, z + sz
                    if (lo[0] <= gx <= hi[0]
                            and lo[1] <= gy <= hi[1]
                            and lo[2] <= gz <= hi[2]):
                        candidates.append(((gx, gy, gz), i,
                                           (ia, ib, ic)))
    if not candidates:
        return empty

    kd = _build_kdtree(positions)
    ghost_positions = []
    ghost_sources = []
    ghost_bonds = []
    ghost_shifts = []

    for (gpos, src, image) in candidates:
        found = []
        if kd is not None:
            for _co, j, dist in kd.find_range(gpos, cutoff):
                if dist > 1e-9 and dist < (radii[j] + radii[src]) * scale:
                    found.append(j)
        else:
            gx, gy, gz = gpos
            for j, (x, y, z) in enumerate(positions):
                dist = ((gx - x) ** 2 + (gy - y) ** 2 + (gz - z) ** 2) ** 0.5
                if 1e-9 < dist < (radii[j] + radii[src]) * scale:
                    found.append(j)
        if found:
            ghost_index = n + len(ghost_positions)
            ghost_positions.append(gpos)
            ghost_sources.append(src)
            ghost_shifts.append(image)
            ghost_bonds.extend((j, ghost_index) for j in found)

    if return_shifts:
        return ghost_positions, ghost_sources, ghost_bonds, ghost_shifts
    return ghost_positions, ghost_sources, ghost_bonds


def compute_bond_directions(positions, bonds):
    """Average bond direction per atom, normalized.

    Returns a list of unit vectors, one per position; atoms without bonds
    (or with perfectly canceling bonds) get ``(0, 0, 1)``. Written as the
    ``BondDir`` attribute so instances can be aligned along their local
    bond axis (arrows, struts, cables).
    """
    n = len(positions)
    sums = [[0.0, 0.0, 0.0] for _ in range(n)]
    for i, j in bonds:
        xi, yi, zi = positions[i]
        xj, yj, zj = positions[j]
        dx, dy, dz = xj - xi, yj - yi, zj - zi
        length = (dx * dx + dy * dy + dz * dz) ** 0.5
        if length < 1e-9:
            continue
        dx, dy, dz = dx / length, dy / length, dz / length
        si, sj = sums[i], sums[j]
        si[0] += dx
        si[1] += dy
        si[2] += dz
        sj[0] += dx
        sj[1] += dy
        sj[2] += dz

    result = []
    for s in sums:
        length = (s[0] ** 2 + s[1] ** 2 + s[2] ** 2) ** 0.5
        if length < 1e-6:
            result.append((0.0, 0.0, 1.0))
        else:
            result.append((s[0] / length, s[1] / length, s[2] / length))
    return result
