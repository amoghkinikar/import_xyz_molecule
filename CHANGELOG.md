# Changelog

## 2.2.0

### Added
- Multi-frame `.xyz` / extxyz files import as an animation (one shape key
  per frame). Bonds and periodic ghost atoms are the union over all
  frames; the node tree hides bonds longer than their per-edge
  `BondCutoff` at the current frame, and ghosts left without a bond.
  Import options: Frames (All / First / Last), Frame Step, Unwrap
  Periodic Jumps. A changing cell is animated too.
- Atom labels such as `C1`, `Au_s`, `Cu12` map to their element, with a
  report of the mapping.
- Figure panel: **Frame** option (whole object / adsorbate only).
- `blender_manifest.toml` (Blender 4.2+ extension), LICENSE (GPL-3.0).

### Changed
- Requires Blender 5.1 or newer (the version it is tested on).
- **Render Figure** now renders a top-down orthographic view framed on
  the visible atoms and their drawn radii. Resolution is the longer
  side; the image aspect follows the structure's footprint.
- Render Figure no longer leaves anything behind: its camera and sun are
  removed after the render, and the scene camera, engine, EEVEE AO /
  raytracing, Cycles samples, and view transform / look are restored.
  The render uses the Standard view transform so element colours are
  not washed out.

### Removed
- The angled perspective figure camera.
