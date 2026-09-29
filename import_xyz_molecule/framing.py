"""Top-down figure framing, pure Python (no bpy), so it can be tested
outside Blender.

The figure camera looks straight down the world -Z axis (image +Y is
world +Y), orthographic. The structure's footprint in XY -- atom centers
plus the radius each atom is drawn with -- sets the orthographic scale
and the image aspect ratio, so a long nanoribbon gets a wide image rather
than a square one with empty bands.
"""


def topdown_frame(points, radii, long_side=1600, margin=0.08,
                  min_side=16):
    """Camera parameters for a top-down orthographic view.

    ``points``: world-space ``(x, y, z)`` atom centers.
    ``radii``: world-space drawn radius per point.
    ``long_side``: pixels of the longer image side.
    ``margin``: empty border on every side, as a fraction of the longer
    footprint side.

    Returns a dict: ``center`` (x, y), ``ortho_scale``, ``resolution``
    (x, y), ``camera_z``, ``clip_start``, ``clip_end``.
    """
    if not points:
        raise ValueError("Nothing visible to frame.")
    xmin = min(p[0] - r for p, r in zip(points, radii))
    xmax = max(p[0] + r for p, r in zip(points, radii))
    ymin = min(p[1] - r for p, r in zip(points, radii))
    ymax = max(p[1] + r for p, r in zip(points, radii))
    zmin = min(p[2] - r for p, r in zip(points, radii))
    zmax = max(p[2] + r for p, r in zip(points, radii))

    width = xmax - xmin
    height = ymax - ymin
    span = max(width, height, 1e-6)
    pad = margin * span
    width += 2.0 * pad
    height += 2.0 * pad
    ortho_scale = max(width, height)

    long_side = int(long_side)
    short = int(round(long_side * min(width, height) / ortho_scale))
    short = max(min_side, min(short, long_side))
    resolution = (long_side, short) if width >= height else \
        (short, long_side)

    depth = zmax - zmin
    standoff = max(depth, ortho_scale) * 0.5 + 1e-3
    camera_z = zmax + standoff
    clip_start = standoff * 0.25
    clip_end = camera_z - zmin + standoff

    return {
        "center": ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0),
        "ortho_scale": ortho_scale,
        "resolution": resolution,
        "camera_z": camera_z,
        "clip_start": clip_start,
        "clip_end": clip_end,
    }


def drawn_radius_factor(style, atom_scale, bond_radius, scale_jitter=0.0):
    """How each style draws an atom, as ``(attribute, multiplier)``:

    the drawn radius is ``attribute value * multiplier`` (attribute None
    means a constant ``multiplier``). Mirrors the node tree in nodes.py.
    """
    jitter = 1.0 + max(0.0, scale_jitter)
    if style == 'SPACE_FILLING':
        return "VdWRadius", atom_scale * jitter
    if style == 'LICORICE':
        return None, bond_radius * jitter
    if style == 'WIREFRAME':
        return None, bond_radius * 0.35
    return "AtomRadius", atom_scale * jitter


def interpolate_absolute(key_times, key_coords, eval_time):
    """Positions of absolute shape keys at ``eval_time`` (linear between
    the two bracketing keys; clamped at the ends)."""
    if eval_time <= key_times[0]:
        return key_coords[0]
    if eval_time >= key_times[-1]:
        return key_coords[-1]
    for k in range(len(key_times) - 1):
        t0, t1 = key_times[k], key_times[k + 1]
        if t0 <= eval_time <= t1:
            f = (eval_time - t0) / (t1 - t0) if t1 > t0 else 0.0
            a, b = key_coords[k], key_coords[k + 1]
            return [
                (p[0] + (q[0] - p[0]) * f, p[1] + (q[1] - p[1]) * f,
                 p[2] + (q[2] - p[2]) * f)
                for p, q in zip(a, b)
            ]
    return key_coords[-1]


def mix_relative(basis, key_coords, values):
    """Positions of relative shape keys: basis + sum(value * delta)."""
    out = [list(p) for p in basis]
    for coords, v in zip(key_coords, values):
        if not v:
            continue
        for i, (q, b) in enumerate(zip(coords, basis)):
            out[i][0] += v * (q[0] - b[0])
            out[i][1] += v * (q[1] - b[1])
            out[i][2] += v * (q[2] - b[2])
    return [tuple(p) for p in out]
