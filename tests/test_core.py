"""Tests for the parts of the add-on that do not need Blender.

Run from the repository root:

    python -m unittest -v tests.test_core

The Geometry Nodes tree, shape-key animation and rendering are covered by
tests/blender_smoke_test.py, which runs inside Blender.
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from import_xyz_molecule import (  # noqa: E402
    bonds, elements, framing, parsing, trajectory,
)

KNOWN = elements.merge_elements({})


def write_tmp(text):
    fd, path = tempfile.mkstemp(suffix=".xyz")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

class TestParsing(unittest.TestCase):

    def test_atom_labels_map_to_elements(self):
        """DFT site labels must not turn into magenta 'X' atoms, and real
        two-letter symbols must win over one-letter ones (Co1 is cobalt)."""
        cases = {"C": ("C", False), "AU": ("Au", False),
                 "C1": ("C", True), "Au_s": ("Au", True),
                 "Co1": ("Co", True), "Hb": ("H", True),
                 "Zz9": ("X", True)}
        for label, expected in cases.items():
            self.assertEqual(elements.resolve_symbol(label, KNOWN), expected,
                             label)

    def test_frames_all_first_last_and_step(self):
        """The Frames option picks the right geometry; for a geometry
        optimisation, Last must be the final (relaxed) frame."""
        text = "".join(f"2\nframe {k}\nC 0 0 0\nO {1.0 + 0.1 * k} 0 0\n"
                       for k in range(4)) + "\n\n"
        path = write_tmp(text)

        def x_of_o(frames):
            return [round(f["positions"][1][0], 6) for f in frames]

        self.assertEqual(x_of_o(parsing.parse_xyz_frames(path, KNOWN, 'ALL')),
                         [1.0, 1.1, 1.2, 1.3])
        self.assertEqual(x_of_o(parsing.parse_xyz_frames(path, KNOWN,
                                                         'FIRST')), [1.0])
        last = parsing.parse_xyz_frames(path, KNOWN, 'LAST')
        self.assertEqual(x_of_o(last), [1.3])
        self.assertEqual(last[0]["n_frames_in_file"], 4)
        self.assertEqual(x_of_o(parsing.parse_xyz_frames(path, KNOWN, 'ALL',
                                                         step=2)), [1.0, 1.2])

    def test_frames_with_different_atoms_are_rejected(self):
        """A mesh has one fixed vertex list, so frames must list the same
        atoms in the same order; anything else must fail loudly."""
        path = write_tmp("2\na\nC 0 0 0\nO 1 0 0\n2\nb\nO 0 0 0\nC 1 0 0\n")
        with self.assertRaises(ValueError):
            parsing.parse_xyz_frames(path, KNOWN, 'ALL')

    def test_lattice_is_read_per_frame(self):
        """Needed for periodic ghosts and the animated cell box in cell
        optimisations / NPT runs."""
        text = "".join(
            f'1\nLattice="{a} 0 0 0 {a} 0 0 0 {a}" '
            f'Properties=species:S:1:pos:R:3 pbc="T T T"\nC 0 0 0\n'
            for a in (10.0, 10.5))
        frames = parsing.parse_xyz_frames(write_tmp(text), KNOWN, 'ALL')
        lattices = trajectory.frame_lattices(frames)
        self.assertEqual([lat[0][0] for lat in lattices], [10.0, 10.5])
        self.assertTrue(trajectory.lattices_vary(lattices))


# ---------------------------------------------------------------------------
# Bonds across frames
# ---------------------------------------------------------------------------

def edge_visible(asm, frame, edge_index):
    """What the node tree does: hide an edge longer than its cutoff."""
    i, j = asm["edges"][edge_index]
    p = asm["positions"][frame]
    length = sum((a - b) ** 2 for a, b in zip(p[i], p[j])) ** 0.5
    return length <= asm["cutoffs"][edge_index]


class TestTrajectory(unittest.TestCase):

    def test_single_frame_gives_the_static_result(self):
        """Static imports go through the trajectory code too; they must get
        exactly the bonds and ghosts the static algorithm finds."""
        lattice = ((2.46, 0.0, 0.0), (-1.23, 2.1304, 0.0), (0.0, 0.0, 15.0))
        pos = [(0.0, 0.0, 0.0), (1.23, 0.7101, 0.0)]
        radii = [0.75, 0.75]
        asm = trajectory.assemble([pos], [lattice], radii, pbc=(1, 1, 0))
        gpos, gsrc, gbonds = bonds.compute_periodic_ghosts(
            pos, radii, lattice, pbc=(1, 1, 0))
        self.assertEqual(asm["edges"], bonds.compute_bonds(pos, radii) + gbonds)
        self.assertEqual(asm["ghost_sources"], gsrc)
        self.assertEqual(asm["positions"][0], pos + gpos)

    def test_bond_forms_then_breaks(self):
        """A bond present only in the middle frame is stored, and shown
        only in that frame."""
        frames = [[(0.0, 0.0, 0.0), (d, 0.0, 0.0)] for d in (3.0, 1.4, 3.0)]
        asm = trajectory.assemble(frames, [None] * 3, [0.75, 0.75])
        self.assertEqual(asm["edges"], [(0, 1)])
        self.assertEqual([edge_visible(asm, k, 0) for k in range(3)],
                         [False, True, False])

    def test_ghost_follows_its_source_atom(self):
        """A periodic ghost first needed in a later frame must exist in
        every frame and sit one lattice vector from its source atom."""
        lattice = ((5.0, 0.0, 0.0), (0.0, 50.0, 0.0), (0.0, 0.0, 50.0))
        f0 = [(1.0, 0.0, 0.0), (2.4, 0.0, 0.0)]     # no periodic bond
        f1 = [(0.2, 0.0, 0.0), (3.8, 0.0, 0.0)]     # bond across x boundary
        asm = trajectory.assemble([f0, f1], [lattice] * 2, [0.75, 0.75],
                                  pbc=(True, False, False))
        n = asm["n_real"]
        self.assertTrue(asm["ghost_sources"])
        for g, src in enumerate(asm["ghost_sources"]):
            for k in range(2):
                dx = asm["positions"][k][n + g][0] - asm["positions"][k][src][0]
                self.assertAlmostEqual(abs(dx), 5.0)

    def test_unwrap_triclinic_cell(self):
        """An atom wrapped back into a hexagonal cell between frames must be
        moved back, or the animation shows it flying across the cell."""
        lattice = ((2.46, 0.0, 0.0), (-1.23, 2.1304, 0.0), (0.0, 0.0, 15.0))
        start = (0.05, 0.05, 1.0)
        wrapped = (0.05 + 0.01 + 1.23, 0.05 + 0.02 - 2.1304, 1.0)  # moved, -b
        out = trajectory.unwrap_positions([[start], [wrapped]], [lattice] * 2)
        self.assertAlmostEqual(out[1][0][0], 0.06)
        self.assertAlmostEqual(out[1][0][1], 0.07)


# ---------------------------------------------------------------------------
# Figure framing
# ---------------------------------------------------------------------------

class TestFraming(unittest.TestCase):

    def test_topdown_frame_fits_a_ribbon(self):
        """Top-down camera: centred on the footprint, atoms' radii and the
        margin included, and a long ribbon gets a wide image."""
        points = [(-10.0, -2.0, 0.0), (10.0, 2.0, 1.0)]
        frame = framing.topdown_frame(points, [0.5, 0.5], long_side=1600,
                                      margin=0.05)
        # footprint 21 x 5 A, plus 5 % of 21 on every side
        self.assertEqual(frame["center"], (0.0, 0.0))
        self.assertAlmostEqual(frame["ortho_scale"], 23.1)
        self.assertEqual(frame["resolution"], (1600, 492))
        # camera above the highest atom, everything within the clip range
        self.assertGreater(frame["camera_z"] - frame["clip_start"], 1.5)
        self.assertLess(frame["camera_z"] - frame["clip_end"], -0.5)


if __name__ == "__main__":
    unittest.main()
