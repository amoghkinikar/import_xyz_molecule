"""Element data: symbols, radii, colors, and user JSON overrides.

The built-in table stores, per element:

* atomic number,
* covalent radius in picometers (used for bond detection and the
  Ball-and-Stick / Licorice styles),
* van der Waals radius in picometers (used for the Space-Filling style;
  values follow Bondi (1964) with the Alvarez (2013) extension for
  elements Bondi did not cover),
* an sRGB hex color following the common CPK/Jmol convention.

``"X"`` is the fallback entry for unknown symbols (magenta, tiny).

A user JSON file can override any of these per element and may also map an
element to a template object in the scene::

    {
        "Au": {"radius": 150, "vdw_radius": 214, "color": "FFD700"},
        "Cu": {"color": "AA4400", "object": "CatMesh"}
    }
"""

import json
import os

# (atomic number, covalent radius pm, van der Waals radius pm, sRGB hex)
_ELEMENT_TABLE = {
    "H":  (1, 32, 110, "FFFFFF"),   "He": (2, 120, 140, "D9FFFF"),
    "Li": (3, 145, 182, "CC80FF"),  "Be": (4, 105, 153, "C2FF00"),
    "B":  (5, 85, 192, "FFB5B5"),   "C":  (6, 75, 170, "211F20"),
    "N":  (7, 71, 155, "3050F8"),   "O":  (8, 63, 152, "FF0D0D"),
    "F":  (9, 64, 147, "90E050"),   "Ne": (10, 160, 154, "B3E3F5"),
    "Na": (11, 180, 227, "AB5CF2"), "Mg": (12, 150, 173, "8AFF00"),
    "Al": (13, 126, 184, "BFA6A6"), "Si": (14, 116, 210, "F0C8A0"),
    "P":  (15, 111, 180, "FF8000"), "S":  (16, 103, 180, "FFFF30"),
    "Cl": (17, 100, 175, "5FF025"), "Ar": (18, 96, 188, "80D1E3"),
    "K":  (19, 220, 275, "8F40D4"), "Ca": (20, 180, 231, "3DFF00"),
    "Sc": (21, 160, 215, "E6E6E6"), "Ti": (22, 140, 211, "BFC2C7"),
    "V":  (23, 135, 207, "A6A6AB"), "Cr": (24, 140, 206, "8A99C7"),
    "Mn": (25, 140, 205, "9C7AC7"), "Fe": (26, 140, 204, "E06633"),
    "Co": (27, 135, 200, "F090A0"), "Ni": (28, 135, 197, "50D050"),
    "Cu": (29, 135, 196, "C88033"), "Zn": (30, 135, 201, "7D80B0"),
    "Ga": (31, 130, 187, "C28F8F"), "Ge": (32, 125, 211, "668F8F"),
    "As": (33, 121, 185, "BD80E3"), "Se": (34, 116, 190, "FFA100"),
    "Br": (35, 115, 185, "A62929"), "Kr": (36, 117, 202, "5CB8D1"),
    "Rb": (37, 235, 303, "702EB0"), "Sr": (38, 200, 249, "00FF00"),
    "Y":  (39, 180, 232, "94FFFF"), "Zr": (40, 155, 223, "94E0E0"),
    "Nb": (41, 147, 218, "73C2C9"), "Mo": (42, 145, 217, "54B5B5"),
    "Tc": (43, 135, 216, "3B9E9E"), "Ru": (44, 130, 213, "248F8F"),
    "Rh": (45, 135, 210, "0A7D8C"), "Pd": (46, 140, 210, "006985"),
    "Ag": (47, 160, 211, "C0C0C0"), "Cd": (48, 155, 218, "FFD98F"),
    "In": (49, 155, 193, "A67573"), "Sn": (50, 145, 217, "668080"),
    "Sb": (51, 145, 206, "9E63B5"), "Te": (52, 140, 206, "D47A00"),
    "I":  (53, 140, 198, "940094"), "Xe": (54, 131, 216, "429EB0"),
    "Cs": (55, 265, 343, "57178F"), "Ba": (56, 215, 268, "00C900"),
    "La": (57, 195, 243, "70D4FF"), "Ce": (58, 185, 242, "FFFFC7"),
    "Pr": (59, 185, 240, "D9FFC7"), "Nd": (60, 185, 239, "C7FFC7"),
    "Pm": (61, 185, 238, "A3FFC7"), "Sm": (62, 185, 236, "8FFFC7"),
    "Eu": (63, 185, 235, "61FFC7"), "Gd": (64, 180, 234, "45FFC7"),
    "Tb": (65, 175, 233, "30FFC7"), "Dy": (66, 175, 231, "1FFFC7"),
    "Ho": (67, 175, 230, "00FF9C"), "Er": (68, 175, 229, "00E675"),
    "Tm": (69, 175, 227, "00D452"), "Yb": (70, 175, 226, "00BF38"),
    "Lu": (71, 175, 224, "00AB24"), "Hf": (72, 155, 223, "4DC2FF"),
    "Ta": (73, 146, 222, "4DA6FF"), "W":  (74, 137, 218, "2194D6"),
    "Re": (75, 135, 216, "267DAB"), "Os": (76, 130, 216, "266696"),
    "Ir": (77, 135, 213, "175487"), "Pt": (78, 135, 213, "D0D0E0"),
    "Au": (79, 135, 214, "FFC152"), "Hg": (80, 150, 223, "B8B8D0"),
    "Tl": (81, 190, 196, "A6544D"), "Pb": (82, 180, 202, "575961"),
    "Bi": (83, 160, 207, "9E4FB5"), "Po": (84, 190, 197, "AB5C00"),
    "At": (85, 147, 202, "754F45"), "Rn": (86, 142, 220, "428296"),
    "Fr": (87, 260, 348, "420066"), "Ra": (88, 215, 283, "007D00"),
    "Ac": (89, 195, 247, "70ABFA"), "Th": (90, 180, 245, "00BAFF"),
    "Pa": (91, 180, 243, "00A1FF"), "U":  (92, 175, 241, "008FFF"),
    "Np": (93, 175, 239, "0080FF"), "Pu": (94, 175, 243, "006BFF"),
    "Am": (95, 175, 244, "545CF2"), "Cm": (96, 176, 245, "785CE3"),
    "X":  (0, 10, 100, "FF00FF"),  # unknown-element fallback
}

# Reverse lookup: atomic number -> symbol (for .xyz files that use numbers).
NUMBER_TO_SYMBOL = {
    num: sym for sym, (num, _r, _v, _c) in _ELEMENT_TABLE.items()
}


def default_elements():
    """Return the built-in element dictionary in mutable dict-of-dicts form."""
    return {
        sym: {
            "number": num,
            "radius": rad,
            "vdw_radius": vdw,
            "color": col,
            "object": "",
        }
        for sym, (num, rad, vdw, col) in _ELEMENT_TABLE.items()
    }


def normalize_symbol(raw):
    """Normalize an element symbol: 'AU' -> 'Au', 'cl' -> 'Cl'."""
    raw = str(raw).strip()
    if not raw:
        return "X"
    return raw[0].upper() + raw[1:].lower()


def resolve_symbol(token, known_symbols):
    """Map an atom label from a structure file onto a known element.

    Returns ``(symbol, mapped)`` where ``mapped`` is True when the label
    was not itself a known symbol. Handles the labels DFT codes and
    hand-edited files use for inequivalent sites / custom kinds::

        'C'    -> ('C', False)       'AU' -> ('Au', False)
        'C1'   -> ('C', True)        'Au_s' -> ('Au', True)
        'Cu12' -> ('Cu', True)       'Hb'  -> ('H', True)

    Two-letter element matches are preferred over one-letter ones, so
    'Co1' is cobalt, not carbon. Labels matching nothing give ``'X'``.
    """
    symbol = normalize_symbol(token)
    if symbol in known_symbols:
        return symbol, False
    letters = ""
    for ch in str(token).strip():
        if ch.isalpha():
            letters += ch
        else:
            break
    for length in (2, 1):
        if len(letters) >= length:
            candidate = normalize_symbol(letters[:length])
            if candidate in known_symbols and candidate != "X":
                return candidate, True
    return "X", True


def load_user_elements(json_path):
    """Load a user JSON override dictionary. Returns ``(dict, error)``.

    Only the provided keys of each element are overridden; everything else
    keeps the built-in values. Unknown element symbols create new entries.
    """
    if not json_path:
        return {}, None
    try:
        import bpy
        path = bpy.path.abspath(json_path)
    except ImportError:  # outside Blender (tests)
        path = json_path
    if not os.path.isfile(path):
        return {}, f"User element file not found: {path}"
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        return {}, f"Could not read user element file: {exc}"
    if not isinstance(data, dict):
        return {}, ("User element file must contain a JSON object keyed by "
                    "element symbol.")

    cleaned = {}
    for sym, entry in data.items():
        if not isinstance(entry, dict):
            return {}, (f"Entry for '{sym}' must be an object, "
                        'e.g. {"radius": 150}.')
        item = {}
        for radius_key in ("radius", "vdw_radius"):
            if radius_key in entry:
                try:
                    item[radius_key] = float(entry[radius_key])
                except (TypeError, ValueError):
                    return {}, (f"'{radius_key}' for '{sym}' must be a "
                                "number (picometers).")
        if "color" in entry:
            col = str(entry["color"]).lstrip("#")
            try:
                int(col, 16)
            except ValueError:
                return {}, (f"'color' for '{sym}' must be a hex string "
                            'like "FFD700".')
            item["color"] = col.zfill(6)
        if "number" in entry:
            try:
                item["number"] = int(entry["number"])
            except (TypeError, ValueError):
                return {}, f"'number' for '{sym}' must be an integer."
        if "object" in entry:
            item["object"] = str(entry["object"])
        cleaned[normalize_symbol(sym)] = item
    return cleaned, None


def merge_elements(user_overrides):
    """Merge user overrides over the built-in table. Returns the final dict."""
    elements = default_elements()
    for sym, entry in (user_overrides or {}).items():
        base = elements.setdefault(
            sym,
            {"number": 0, "radius": 10, "vdw_radius": 100,
             "color": "FF00FF", "object": ""},
        )
        base.update(entry)
    return elements


# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------

def srgb_to_linear(c):
    if c < 0.0:
        return 0.0
    if c < 0.04045:
        return c / 12.92
    return ((c + 0.055) / 1.055) ** 2.4


def hex_to_linear_rgb(hex_string, alpha=1.0):
    """'FFC152' -> (r, g, b, a) in linear (scene) color space."""
    h = int(hex_string, 16)
    r = (h & 0xFF0000) >> 16
    g = (h & 0x00FF00) >> 8
    b = h & 0x0000FF
    return tuple(srgb_to_linear(c / 255.0) for c in (r, g, b)) + (alpha,)


def pseudo_element_color(index):
    """Deterministic, well-separated sRGB hex for the Nth unknown CSV label.

    Uses golden-ratio hue stepping so labels that are not chemical elements
    (e.g. CSV categories) each get a distinct, saturated color.
    """
    import colorsys
    hue = (index * 0.618033988749895) % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.95)
    return f"{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}"
