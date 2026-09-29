# XYZ Molecule Importer

Import `.xyz` / extended-xyz molecular structures — and generic
`label,x,y,z` point CSVs — into Blender as Geometry Nodes objects.
Built for scientists who do not use Blender: one import click gives a
centered, framed, correctly colored molecule; one more click renders a
publication-ready figure.

## Install

Requires Blender 5.1 or newer.

1. Open the [Releases page](https://github.com/amoghkinikar/import_xyz_molecule/releases)
   and, under **Assets** of the latest release, download
   `import_xyz_molecule-<version>.zip`. Do not unzip it.
2. In Blender: `Edit > Preferences > Get Extensions`, open the dropdown
   at the top right, choose **Install from Disk...** and select the zip.
3. The importer is now under `File > Import > XYZ Molecule (.xyz)` and
   in the **XYZ Molecule** tab of the 3D viewport sidebar (`N`).

> Use the zip from the Releases page, not the green **Code > Download
> ZIP** button: that one is the source code, and Blender will not
> install it as an add-on.

To update, install the newer release zip the same way; it replaces the
old version.

## Quick start

1. `File > Import > XYZ Molecule (.xyz)` (or use the **XYZ Molecule**
   tab in the 3D viewport sidebar, hotkey `N`).
2. The molecule is centered at the origin, the viewport zooms to it,
   and the default cube is removed (all optional).
3. In the sidebar: pick a display style (Ball-and-Stick, Space-Filling,
   Licorice, Wireframe), hide hydrogens, or hit **Render Figure** — a
   top-down view of the structure, saved as a PNG next to your `.xyz`
   file.

## Feature overview (sidebar panels)

- **Display** — style presets, atom/bond sizes, show/hide atoms, bonds
  and periodic ghost atoms; **Reload From File** re-imports the same
  file (same settings) after your calculation finishes;
  **Build Supercell** repeats a periodic structure (see below).
- **Elements** — per-element visibility, Hide Hydrogens, select all
  atoms of an element.
- **Figure** — a bird's-eye (top-down, orthographic) render framed on
  the visible atoms. **Frame** chooses the whole object or only the
  adsorbate (atoms that are not slab elements). **Resolution** is the
  longer image side; the shorter side follows the structure's
  footprint, so a long ribbon gives a wide image. Transparent or white
  background, optional Cycles quality. The camera and light exist only
  for the render, and your scene camera and render/colour settings are
  restored afterwards.
- **Color** — color atoms and bonds by any scalar per-atom attribute
  with viridis / plasma / coolwarm, or switch back to CPK element
  colors. The dropdown lists the built-in physical attributes (atomic
  number, covalent and van der Waals radius, element index, random
  value) followed by every column imported from the file (charges,
  spins, ...).
- **Instances** — swap the instanced object for any element (spheres →
  cats, trees, arrows...), per-instance random rotation / scale jitter
  with a seed, and alignment of instances to a vector attribute such as
  `forces` or the computed `BondDir`.

## Extended XYZ (ASE)

`Lattice="..."` in the comment line draws the unit cell as a wireframe
box and adds periodic bonds across cell boundaries (ghost atoms, tagged
with the `IsGhost` attribute and toggleable). Extra columns declared in
`Properties=...` (forces, charges, spins, ...) are imported as point
attributes automatically and can drive colors and orientations.

## Multi-frame files (trajectories)

An `.xyz` holding several frames (MD, a CP2K `GEO_OPT` `*-pos-1.xyz`,
an ASE trajectory) is imported as an animation by default: frame *k* of
the file is scene frame *k* (from the scene start frame), and the scene
end frame is extended if needed. The import dialog's **Frames** option
can instead take only the **First** or **Last** frame (for a geometry
optimisation, *Last* is the relaxed structure); **Frame Step** keeps
every Nth frame of long runs.

Bonds follow the atoms: every bond found in any frame is stored, and at
each frame the ones longer than their cutoff are hidden, so bonds form
and break at the right time. Periodic ghost atoms work the same way.
Atoms that an MD code wrapped back into the cell are unwrapped between
frames (**Unwrap Periodic Jumps**), so they do not fly across the cell.
If the cell itself changes (cell optimisation, NPT), the cell box is
animated too.

All frames must list the same atoms in the same order. Extra columns
(forces, charges) and `BondDir` are taken from the first frame.

## Atom labels

Labels that are not plain element symbols, as used for inequivalent
sites or custom kinds (`C1`, `Au_s`, `Cu12`, `Hb`), are mapped to their
element (two-letter symbols first, so `Co1` is cobalt) and the import
reports the mapping. Anything that matches no element becomes the
magenta `X` placeholder, with a warning.

## Supercells

For a structure with a `Lattice`, the **Build Supercell** button in the
Display panel repeats it along its lattice vectors. Enter repeats along
a, b and c (the dialog shows the resulting atom count before you commit)
and the object is rebuilt from its source file:

- bonds form normally throughout the interior, and periodic ghost stubs
  move out to the new outer boundary;
- per-atom attributes (charges, forces, ...) are replicated with their
  atoms, so color-by-attribute keeps working;
- the cell is drawn as a grid of the original cells, so the periodicity
  stays visible.

Repeats are stored on the object: **Reload From File** keeps the current
supercell (handy when a calculation finishes and you want the same 4x4
view of the new geometry), and running Build Supercell again *replaces*
the repeats rather than compounding them — so 1x1x1 returns to the
primitive cell. The limit is 200,000 atoms.

A plain `.xyz` has no cell, so the button stays hidden for those files;
give the file a `Lattice="..."` comment line (any ASE export does this)
to enable it.

## CSV import

`File > Import > Points CSV (.csv)`. A header row is auto-detected;
labels matching chemical elements get chemistry treatment, other labels
get distinct auto-assigned colors; extra numeric columns become point
attributes. Bond detection is optional and off by default.

## Element overrides (JSON)

Set a JSON file in the import dialog or the add-on preferences:

```json
{
    "Au": {"radius": 150, "vdw_radius": 214, "color": "FFD700"},
    "Cu": {"color": "AA4400", "object": "CatMesh"}
}
```

Radii in picometers, colors as sRGB hex. `"object"` names a scene
object to instance for that element right at import time.

## Attributes written per atom

`AtomicNumber`, `ElementColor`, `AtomRadius` (covalent, Å),
`VdWRadius` (Å), `IndexOfElement`, `IsGhost`, `ElementHidden`,
`RandomValue` (per-point random, for material variation), `BondDir`
(average bond direction), plus everything from extxyz/CSV columns.
Per bond (edge): `BondCutoff` (Å), used to hide stretched bonds in
animations.

## Development

Tests are in `tests/` (not part of the installed add-on). From the
repository root:

```
python -m unittest -v tests.test_core
blender -b --factory-startup --python-exit-code 1 -P tests/blender_smoke_test.py
```

The first needs only Python; the second runs inside Blender and checks
import, the node tree, trajectory animation and the figure render.

## License

GPL-3.0-or-later. Inspired by
[Molecular Nodes](https://github.com/BradyAJohnston/MolecularNodes).
