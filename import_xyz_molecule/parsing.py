"""File parsers: .xyz, extended XYZ (ASE extxyz), and generic CSV.

All parsers return a plain ``dict`` (``structure``) with the keys:

* ``symbols``   – list of normalized element symbols / labels,
* ``positions`` – list of ``(x, y, z)`` floats (Angstrom for xyz),
* ``comment``   – the raw comment line (xyz) or ``""``,
* ``lattice``   – ``((ax, ay, az), (bx, by, bz), (cx, cy, cz))`` or ``None``,
* ``pbc``       – ``(bool, bool, bool)`` or ``None``,
* ``extra``     – ``{name: {"type": 'FLOAT'|'INT'|'FLOAT_VECTOR',
  "values": [...]}}`` per-atom data (forces, charges, ...), later written
  as mesh point attributes.

Raises ``ValueError`` with a human-readable message on malformed files.
"""

import csv as _csv
import re

from .elements import NUMBER_TO_SYMBOL, normalize_symbol, resolve_symbol

# Attribute names that the importer writes itself; extra columns are renamed
# if they collide.
RESERVED_ATTRIBUTES = {
    "AtomicNumber", "ElementColor", "AtomRadius", "VdWRadius",
    "IndexOfElement", "RandomValue", "BondDir", "ElementHidden", "IsGhost",
    "position",
}

_KEY_VALUE_RE = re.compile(
    r"""(?P<key>[A-Za-z_][A-Za-z0-9_-]*)      # key
        =
        (?:"(?P<quoted>[^"]*)"                # "quoted value"
          |'(?P<squoted>[^']*)'               # 'quoted value'
          |(?P<plain>\S+))                    # bare value
    """,
    re.VERBOSE,
)

_TRUE_STRINGS = {"t", "true", "1", "yes"}
_FALSE_STRINGS = {"f", "false", "0", "no"}


def _new_structure():
    return {
        "symbols": [],
        "positions": [],
        "comment": "",
        "lattice": None,
        "pbc": None,
        "extra": {},
    }


def _safe_attr_name(name, taken):
    """Sanitize a column name so it is a usable, non-colliding attribute."""
    name = re.sub(r"[^A-Za-z0-9_]", "_", str(name).strip()) or "column"
    base = name
    i = 1
    while name in RESERVED_ATTRIBUTES or name in taken:
        name = f"{base}_{i}"
        i += 1
    return name


# ---------------------------------------------------------------------------
# Extended XYZ comment line
# ---------------------------------------------------------------------------

def parse_extxyz_comment(comment):
    """Parse an extxyz comment line.

    Returns ``(lattice, pbc, properties, warnings)`` where ``properties``
    is a list of ``(name, kind, ncols)`` triples (kind in ``'S' 'R' 'I'
    'L'``) or ``None`` when no usable ``Properties=`` entry exists.

    The comment line of a *plain* .xyz file is free text, so entries that
    merely look like extxyz keys but do not parse (e.g. ``lattice=fcc``)
    are ignored with a warning instead of aborting the import.
    """
    lattice = None
    pbc = None
    properties = None
    warnings = []

    for match in _KEY_VALUE_RE.finditer(comment):
        key = match.group("key")
        value = (match.group("quoted")
                 if match.group("quoted") is not None
                 else match.group("squoted")
                 if match.group("squoted") is not None
                 else match.group("plain"))
        low = key.lower()
        if low == "lattice":
            numbers = value.replace(",", " ").split()
            try:
                v = [float(x) for x in numbers]
            except ValueError:
                v = []
            if len(v) == 9:
                lattice = (tuple(v[0:3]), tuple(v[3:6]), tuple(v[6:9]))
            else:
                warnings.append(
                    f"Ignored '{key}={value}' in the comment line "
                    "(a Lattice needs 9 numbers)"
                )
        elif low == "pbc":
            flags = value.replace(",", " ").split()
            if len(flags) == 3:
                pbc = tuple(f.strip().lower() in _TRUE_STRINGS
                            for f in flags)
        elif low == "properties":
            fields = value.split(":")
            parsed = []
            ok = len(fields) % 3 == 0 and fields
            if ok:
                for i in range(0, len(fields), 3):
                    name, kind = fields[i], fields[i + 1]
                    try:
                        ncols = int(fields[i + 2])
                    except ValueError:
                        ok = False
                        break
                    parsed.append((name, kind.upper(), ncols))
            if ok:
                properties = parsed
            else:
                warnings.append(
                    f"Ignored malformed 'Properties={value}' in the "
                    "comment line"
                )
    return lattice, pbc, properties, warnings


# ---------------------------------------------------------------------------
# XYZ / extxyz
# ---------------------------------------------------------------------------

def parse_xyz(filepath, known_symbols=None):
    """Parse the first frame of an .xyz or extended-xyz file.

    ``known_symbols``: container of valid element symbols; labels that are
    not symbols are mapped onto one (``C1`` -> ``C``) or to ``"X"``. Pass
    ``None`` to keep labels as-is (CSV-style behavior).
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        return parse_xyz_stream(fh, known_symbols=known_symbols)


def parse_xyz_frames(filepath, known_symbols=None, which='ALL', step=1):
    """Parse a (possibly multi-frame) .xyz / extxyz file.

    ``which``: ``'ALL'`` (every ``step``-th frame), ``'FIRST'`` or
    ``'LAST'``. Returns a list of structure dicts (one per kept frame).
    Every frame must list the same atoms in the same order; otherwise a
    ValueError explains which frame differs. The first returned structure
    carries ``n_frames_in_file`` and the collected ``warnings``.
    """
    step = max(1, int(step))
    frames = []
    total = 0
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        first = parse_xyz_stream(fh, known_symbols=known_symbols)
        total = 1
        kept_last = first
        if which in ('ALL', 'FIRST'):
            frames.append(first)
        while which != 'FIRST':
            line = fh.readline()
            while line and not line.strip():
                line = fh.readline()
            if not line:
                break
            try:
                int(float(line))
            except ValueError:
                first.setdefault("warnings", []).append(
                    f"Stopped reading after frame {total}: unexpected "
                    f"line {line.strip()[:40]!r}"
                )
                break
            frame = parse_xyz_stream(fh, known_symbols=known_symbols,
                                     first_line=line)
            total += 1
            if frame["symbols"] != first["symbols"]:
                raise ValueError(
                    f"Frame {total} does not list the same atoms as frame "
                    "1 (different count or order). Import a single frame "
                    "instead (Frames: First or Last)."
                )
            kept_last = frame
            if which == 'ALL' and (total - 1) % step == 0:
                frames.append(frame)
    if which == 'LAST':
        frames = [kept_last]
    if frames[0] is not first:
        # keep file-level notes on the structure the caller looks at
        frames[0]["warnings"] = list(first.get("warnings", []))
        frames[0]["label_map"] = first.get("label_map", {})
    frames[0]["n_frames_in_file"] = total
    return frames


def parse_xyz_stream(fh, known_symbols=None, first_line=None):
    """Parse one frame from ``fh``. ``first_line``: the atom-count line,
    when the caller has already read it."""
    first = fh.readline() if first_line is None else first_line
    try:
        num_atoms = int(float(first))
    except ValueError:
        raise ValueError(
            "The first line of an .xyz file must be the atom count; "
            f"got: {first.strip()!r}"
        )

    structure = _new_structure()
    comment = fh.readline().rstrip("\n")
    structure["comment"] = comment

    lattice, pbc, properties, comment_warnings = \
        parse_extxyz_comment(comment)
    if comment_warnings:
        structure.setdefault("warnings", []).extend(comment_warnings)
    structure["lattice"] = lattice
    if pbc is not None:
        structure["pbc"] = pbc
    elif lattice is not None:
        structure["pbc"] = (True, True, True)

    # Column layout. Default plain xyz: species, x, y, z, then extras.
    if properties is None:
        properties = [("species", "S", 1), ("pos", "R", 3)]

    # Per-property column offsets.
    layout = []
    offset = 0
    for name, kind, ncols in properties:
        layout.append((name, kind, ncols, offset))
        offset += ncols
    min_columns = offset

    species_entry = next(
        (p for p in layout if p[0].lower() in ("species", "symbols", "element")
         and p[1] == "S"),
        None,
    )
    pos_entry = next(
        (p for p in layout
         if p[0].lower() in ("pos", "positions") and p[1] == "R"
         and p[2] == 3),
        None,
    )
    if species_entry is None or pos_entry is None:
        raise ValueError(
            "Properties entry must include 'species:S:1' and 'pos:R:3'; "
            f"got: {[p[0] for p in layout]}"
        )

    # Extra per-atom arrays (everything that is not species/pos).
    extra_defs = []      # (attr_name, kind, ncols, offset)
    taken_names = set()
    for name, kind, ncols, off in layout:
        if (name, kind, ncols, off) in (species_entry, pos_entry):
            continue
        if kind == "S":
            continue  # Blender has no useful per-point string attributes
        attr = _safe_attr_name(name, taken_names)
        taken_names.add(attr)
        extra_defs.append((attr, kind, ncols, off))

    extra_values = {attr: [] for attr, _k, _n, _o in extra_defs}
    unnamed_extra_count = 0  # plain xyz files with trailing numeric columns
    label_map = {}           # file label -> element, when they differ

    sp_off = species_entry[3]
    pos_off = pos_entry[3]

    for i in range(num_atoms):
        line = fh.readline()
        if not line:
            raise ValueError(
                f"File ended early: expected {num_atoms} atoms, got {i}."
            )
        parts = line.split()
        if len(parts) < min_columns:
            raise ValueError(
                f"Atom line {i + 1} has {len(parts)} columns, "
                f"expected at least {min_columns}: {line.strip()!r}"
            )

        token = parts[sp_off]
        if token.isnumeric():
            symbol = NUMBER_TO_SYMBOL.get(int(token), "X")
        elif known_symbols is not None:
            symbol, mapped = resolve_symbol(token, known_symbols)
            if mapped:
                label_map[token] = symbol
        else:
            symbol = normalize_symbol(token)
        try:
            x = float(parts[pos_off])
            y = float(parts[pos_off + 1])
            z = float(parts[pos_off + 2])
        except ValueError:
            raise ValueError(
                f"Atom line {i + 1} has non-numeric coordinates: "
                f"{line.strip()!r}"
            )
        structure["symbols"].append(symbol)
        structure["positions"].append((x, y, z))

        for attr, kind, ncols, off in extra_defs:
            try:
                if kind == "L":
                    vals = [1 if parts[off + c].strip().lower()
                            in _TRUE_STRINGS else 0 for c in range(ncols)]
                elif kind == "I":
                    vals = [int(float(parts[off + c])) for c in range(ncols)]
                else:  # 'R'
                    vals = [float(parts[off + c]) for c in range(ncols)]
            except (ValueError, IndexError):
                raise ValueError(
                    f"Atom line {i + 1}: could not read column(s) for "
                    f"'{attr}': {line.strip()!r}"
                )
            extra_values[attr].append(vals)

        # Plain xyz with anonymous trailing numeric columns.
        if not extra_defs and len(parts) > min_columns:
            unnamed_extra_count = max(
                unnamed_extra_count, len(parts) - min_columns
            )

    # Anonymous trailing columns (plain xyz only) cannot be named, so they
    # are skipped with a warning the operator reports.
    if unnamed_extra_count:
        structure.setdefault("warnings", []).append(
            f"{unnamed_extra_count} unnamed trailing column(s) ignored "
            "(add a Properties=... entry to import them)"
        )

    structure["label_map"] = label_map
    if label_map:
        pairs = ", ".join(f"{k}->{v}" for k, v in sorted(label_map.items()))
        structure.setdefault("warnings", []).append(
            f"Atom labels mapped to elements: {pairs}"
        )

    # Pack extra arrays into attribute specs.
    for attr, kind, ncols, _off in extra_defs:
        rows = extra_values[attr]
        if ncols == 1:
            structure["extra"][attr] = {
                "type": "INT" if kind in ("I", "L") else "FLOAT",
                "values": [row[0] for row in rows],
            }
        elif ncols == 3 and kind == "R":
            structure["extra"][attr] = {
                "type": "FLOAT_VECTOR",
                "values": [tuple(row) for row in rows],
            }
        else:  # split N-component data into scalar attributes
            for c in range(ncols):
                structure["extra"][f"{attr}_{c}"] = {
                    "type": "INT" if kind in ("I", "L") else "FLOAT",
                    "values": [row[c] for row in rows],
                }
    return structure


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def _looks_numeric(token):
    try:
        float(token)
        return True
    except (TypeError, ValueError):
        return False


def parse_csv(filepath):
    """Parse a ``label,x,y,z,extra...`` CSV/TSV file into a structure dict.

    * Delimiter is sniffed (comma, semicolon, tab, or whitespace).
    * A header row is auto-detected; header names become attribute names.
    * The label column is optional -- pure ``x,y,z`` files get label "X".
    * Extra numeric columns become per-point FLOAT attributes.
    """
    with open(filepath, "r", encoding="utf-8-sig", errors="replace",
              newline="") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        try:
            dialect = _csv.Sniffer().sniff(sample, delimiters=",;\t ")
            reader = _csv.reader(fh, dialect)
        except _csv.Error:
            # Fall back to whitespace splitting.
            reader = ([t for t in line.split()] for line in fh)

        rows = [
            [cell.strip() for cell in row]
            for row in reader
            if row and any(cell.strip() for cell in row)
        ]

    if not rows:
        raise ValueError("CSV file contains no data rows.")

    # Header detection: first row counts as a header when any of the cells
    # that should be coordinates are non-numeric AND the second row looks
    # like data.
    header = None
    first, rest = rows[0], rows[1:]
    numeric_in_first = sum(_looks_numeric(c) for c in first)
    if rest and numeric_in_first < max(3, len(first) - 1):
        header = first
        rows = rest
    if not rows:
        raise ValueError("CSV file has a header but no data rows.")

    ncols = max(len(r) for r in rows)
    if ncols < 3:
        raise ValueError(
            "CSV needs at least 3 columns (x, y, z), optionally preceded "
            "by a label column."
        )

    # Column roles. Prefer header names; else use position + numeric-ness.
    label_col = None
    xyz_cols = None
    if header:
        lower = [h.lower() for h in header]

        def find(*names):
            for n in names:
                if n in lower:
                    return lower.index(n)
            return None

        xi, yi, zi = find("x", "px", "pos_x"), find("y", "py", "pos_y"), \
            find("z", "pz", "pos_z")
        if xi is not None and yi is not None and zi is not None:
            xyz_cols = (xi, yi, zi)
        label_col = find("label", "element", "symbol", "name", "species",
                         "type", "id")

    if xyz_cols is None:
        # First row decides: a leading non-numeric column is the label
        # (unless the header already named one elsewhere).
        probe = rows[0]
        if label_col is None and not _looks_numeric(probe[0]):
            label_col = 0
        candidates = [c for c in range(ncols) if c != label_col]
        if len(candidates) >= 3:
            xyz_cols = tuple(candidates[:3])
    if xyz_cols is None or max(xyz_cols) >= ncols:
        raise ValueError("Could not locate x, y, z columns in the CSV.")

    used = set(xyz_cols) | ({label_col} if label_col is not None else set())
    extra_cols = [c for c in range(ncols) if c not in used]

    taken = set()
    extra_names = {}
    for c in extra_cols:
        raw = (header[c] if header and c < len(header) and header[c]
               else f"column_{c + 1}")
        name = _safe_attr_name(raw, taken)
        taken.add(name)
        extra_names[c] = name

    structure = _new_structure()
    extra_data = {name: [] for name in extra_names.values()}
    extra_is_int = {name: True for name in extra_names.values()}

    for i, row in enumerate(rows):
        if len(row) < max(xyz_cols) + 1:
            raise ValueError(
                f"CSV row {i + 1} has too few columns: {row!r}"
            )
        try:
            pos = tuple(float(row[c]) for c in xyz_cols)
        except ValueError:
            raise ValueError(
                f"CSV row {i + 1} has non-numeric coordinates: {row!r}"
            )
        label = row[label_col] if label_col is not None and \
            label_col < len(row) and row[label_col] else "X"
        structure["symbols"].append(label)
        structure["positions"].append(pos)

        for c, name in extra_names.items():
            token = row[c] if c < len(row) else ""
            if _looks_numeric(token):
                val = float(token)
                if val != int(val):
                    extra_is_int[name] = False
                extra_data[name].append(val)
            else:
                extra_data[name].append(0.0)
                if token:
                    extra_is_int[name] = False

    for name, values in extra_data.items():
        if extra_is_int[name]:
            structure["extra"][name] = {
                "type": "INT", "values": [int(v) for v in values],
            }
        else:
            structure["extra"][name] = {"type": "FLOAT", "values": values}
    return structure
