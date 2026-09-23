# Geometry provenance — `turbinesFoam/geometry/`

Everything committed under `turbinesFoam/geometry/` is generated from committed
sources by `src/makeGeometry.py` through its per-component builder registry
(`BUILDERS`): the `nacelle` is meshed by the pinned gmsh from `src/nacelle.geo`,
the `phaseVI_blade` is built by the pure-Python structured loft in
`src/blade_phasevi.py`. This file records where each geometry comes from, what
generated it, and the hashes that link the committed artefacts to their inputs.

## Sources

| Component | Source | Builder |
|---|---|---|
| `nacelle` | Yang & Sotiropoulos, *A new class of actuator surface models for wind turbines*, arXiv:[1702.02108v4](https://arxiv.org/abs/1702.02108v4) (2018), Sec. 4.1 ("Flow over periodically placed nacelles") and Fig. 3 (schematic of the nacelle geometry with the triangular surface mesh) | gmsh, `src/nacelle.geo` |
| `phaseVI_blade` | NREL Phase VI rotor: Hand, Simms, Fingersh, Jager, Cotrell, Schreck & Larwood, *Unsteady Aerodynamics Experiment Phase VI: Wind Tunnel Test Configurations and Available Data Campaigns*, NREL/TP-500-29955 (2001), Table A-1 (the 26-station blade geometry, committed as `validation/phaseVI/data/geometry/phaseVI_blade.csv`); S809 profile: Somers, *Design and Experimental Results for the S809 Airfoil*, NREL/SR-440-6918 (1997), Table 2 (committed as `validation/phaseVI/data/s809/s809_somers_nlr.csv`, with its own `PROVENANCE.md`) | pure Python, `src/blade_phasevi.py` |

The MEXICO rotor blade is out of scope here: a future MEXICO blade is a separate
`mexico_blade` component owned by the `mexico-validation` change, not a
`blade0/1/2` reservation. The retired `blade0/1/2` names are documented in
`README.md` and rejected by the CLI; no per-blade STL is generated or committed.

The paper's nacelle is *"composed of a hemisphere facing upstream and a circular
cylinder in its downstream"* (Sec. 4.1) at `Re = 1000` based on the cylinder
radius `R` and the freestream velocity, in a `30R × 20R × 20R` domain with the
downstream end of the body at the streamwise origin.

### Cylinder length `L` (design §12 open item Q1)

Sec. 4.1 gives `R` but not the cylinder length. Fig. 3 dimensions the body as
`R` (hemisphere), `6R` (cylinder) and `2R` (diameter), so the model is

```
nose  <-- hemisphere, length R -->|<- cylinder, length L = 6R ->|  downstream end
x = -7R                            x = -6R                        x = 0
```

The figure labels are unambiguous, and a pixel cross-check of its dimension
arrows agrees: the `R` span in the exported figure is 260 px (`x = 29..289`) and
the `6R` span is 1568 px (`x = 289..1857`), a ratio of **6.03**. No assumed
aspect ratio was needed. The frame places the downstream end at `x = 0` and the
nose upstream, matching the paper's case setup.

Canonical units are metres with `R = 1 m` (the validation case's `Re = 1000` at
`U∞ = 1 m/s`, `ν = 1e-3 m²/s`).

## Generator

| Item | `nacelle` (gmsh route) | `phaseVI_blade` (Python route) |
|---|---|---|
| route | `gmsh` (pinned) | `python` (`src/blade_phasevi.py`) |
| gmsh | **4.15.2** (pip package `gmsh`, `libgmsh` 4.15) | not used |
| mesher options | `Mesh.Algorithm = 6` (Frontal-Delaunay), `Mesh.CharacteristicLengthMin = Mesh.CharacteristicLengthMax = 0.2125·R`, `Mesh.Binary = 1` | fixed 25-segment × 60-point structured loft, no mesher |
| CAD kernel | OpenCASCADE (`SetFactory("OpenCASCADE")`), hemisphere ∪ cylinder | n/a |
| geometry source | `src/nacelle.geo` (parameters `R = 1`, `cylinderRatio = 6`, `meshSizeRatio = 0.2125`) | committed station table + committed S809 coordinates |
| generator script | `src/makeGeometry.py` | `src/makeGeometry.py` → `src/blade_phasevi.py` |
| output | `stl/nacelle.stl` (binary STL), `metadata/nacelle.json` | `stl/phaseVI_blade.stl` (binary STL), `metadata/phaseVI_blade.json` |

Regeneration is byte-identical for both components (verified: two consecutive
runs produced the same sha256 for each STL and each metadata file).
`makeGeometry.py --check` re-runs the registered builder into a temporary
directory and byte-compares it against the committed artefacts; it never writes
to the tree. `--check --component phaseVI_blade` needs no gmsh at all.

gmsh's *element ordering* is not stable across environments: the same mesh
(same vertices, same triangles) is emitted in a different order depending on the
ambient environment, which would break byte-identical regeneration. The
generator therefore re-serializes the triangles before writing: every triangle
is rotated into its lexicographically smallest cyclic form (preserving the
winding, hence the outward normal) and the triangles are sorted by vertex
coordinates. The committed STLs also carry a fixed per-component header and
recomputed facet normals, so their bytes depend only on the mesh geometry —
verified identical in the plain and in the OpenFOAM
(`WM_PROJECT_DIR/etc/bashrc`) environments, and a duplicate triangle would be
rejected. The Python loft emits the same canonical form directly.

### Nacelle triangle count

The committed nacelle surface has **2616 triangles**. The paper reports 2652
triangles for its nacelle surface (Sec. 4.1) — a difference of **−1.4 %**,
within the `~2652` target. The count is set by the `meshSizeRatio` seed: with the
pinned gmsh the triangle count moves in plateaus as the seed varies
(`0.2095..0.2145` all give 2590–2616 triangles), so 2616 is the closest
reproducible value; the seed is deliberately not over-tuned to a single count.

Mesh statistics of the committed nacelle STL (all consistent with an *inscribed*
polyhedral approximation of the analytic surface):

| Quantity | Committed mesh | Analytic |
|---|---|---|
| triangle area sum | 46.9962 | `15πR² = 47.1239` (−0.27 %) |
| signed volume (outward) | +20.8102 | `(2/3)πR³ + πR²L = 20.9439` (−0.64 %) |
| bounds | `x ∈ [-6.996, 0]`, `y, z ∈ [-1, 1]` | `x ∈ [-7, 0]`, `y, z ∈ [-1, 1]` |

The extreme upstream point is a polygon ring just inside the analytic nose
(`-6.996` instead of `-7`): a triangulated sphere only reaches its pole if a
vertex sits exactly there, and the Frontal-Delaunay mesh puts none. The surface
is closed, manifold and consistently outward-oriented (checked by
`tests/test_nacelle_data.py`).

### Blade surface

The committed `phaseVI_blade` surface has **3000 triangles** (25 segments × 60
profile points × 2) and covers the full wetted span `r = 0.5083 → 5.029 m`,
including the cylinder stations at `r ≤ 0.8835 m` and the single transition
segment to the first S809 station at `1.0085 m`. It is the wetted surface only —
no root or tip caps; tip force zeroing, end effects and the inboard `chord_mount`
transition remain BEM-side.

| Quantity | Committed mesh |
|---|---|
| triangle area sum | 4.88423 m² |
| bounds | `x ∈ [-0.1061, 0.2163]`, `y ∈ [-0.5095, 0.1703]`, `z ∈ [0.5083, 5.029]` m |

The surface is an open, edge-manifold shell with exactly two boundary rings (one
60-edge ring at each end), a fixed outward winding and no duplicate or
degenerate triangles. The generation frame is documented in
`src/blade_phasevi.py` (design §4.4): `+z` is the outward span, `+y` the
reference chord direction (trailing edge → leading edge before twist/pitch) and
`+x = y × z` the section normal, with the chord-line quarter-chord point on the
radial reference line. The declared Phase VI pitch is `pitch_deg = 3.0`,
asserted against `validation/phaseVI/config/case.yaml` (`turbine.pitch_deg`) by
the blade data tests.

## Hashes

| Input / artefact | sha256 |
|---|---|
| nacelle input (`src/nacelle.geo` + parameter block) | `e0a69771736aa3b26a11b00880ddd0acf2cf45d8782f13d56fd78ccb89d20de5` |
| `stl/nacelle.stl` | `18c7beef9e323ba28da20f5c813232c0059693f221a594fd913c4e69c14cb74c` |
| `metadata/nacelle.json` | `f237fc40840e9d0970583e32f3b4f21427ab8c29d83755b8f6b55ebb37437b2d` |
| blade input (`phaseVI_blade.csv` + `s809_somers_nlr.csv` + parameter block) | `7a43047dc0728a04fb3ea5ca9510119242e9aca390bde6895c023140500fe456` |
| blade input file `validation/phaseVI/data/geometry/phaseVI_blade.csv` | `e678d6c52df9c8562361d61c65ad373535ea2351448874213b706bc6824e1576` |
| blade input file `validation/phaseVI/data/s809/s809_somers_nlr.csv` | `861f3fe5bbda4d5e5d8985e47c5805bada76fc8e1b2157b944ff78d500f4af8b` |
| `stl/phaseVI_blade.stl` | `77def499047fe0633e57564247a0fdd1330afbf4148b707db08ec79af0a2ec27` |
| `metadata/phaseVI_blade.json` | `01f75d17e2febb296a7b093f20ff011eb583de52e8f510815f876fcb8169e7dc` |

The nacelle `input_sha256` is computed over the raw `src/nacelle.geo` bytes, a
`NUL` byte, and the canonical parameter block
(`json.dumps(parameters, sort_keys=True, separators=(",", ":"))`), where
`parameters` is the `primary:` block declared in the `.geo`
(`{"R": 1.0, "cylinderRatio": 6.0, "meshSizeRatio": 0.2125}`). The blade
`input_sha256` is computed over the committed `phaseVI_blade.csv` bytes, a `NUL`
byte, the committed `s809_somers_nlr.csv` bytes, a `NUL` byte, and the canonical
parameter block (`{"airfoil_radius_m": 1.0085, "cylinder_radius_m": 0.8835,
"n_profile_points": 60, "n_stations": 26, "pitch_deg": 3.0, "sigma": -1.0}`).
Both recipes are implemented in `makeGeometry`/`blade_phasevi` and re-checked by
`--check` and the test suite.

## Environment requirement

The pip `gmsh` module links `libGLU.so.1`/`libGL.so.1` (and the MPI stack).
`import gmsh` fails — and therefore the nacelle route of `src/makeGeometry.py`
cannot regenerate — unless the system library environment is present. On
SDumont:

```sh
module load glu/9.0.2_gnu          # provides libGLU.so.1 for libgmsh
```

Without a runnable gmsh, the nacelle regeneration checks **skip** with an
explicit reason and still run the pure-Python checks on the committed artefacts.
The `phaseVI_blade` route has no gmsh or GL requirement.

## Reproducing

```sh
cd turbinesFoam/geometry
python3 src/makeGeometry.py --check --component phaseVI_blade   # blade, no gmsh
python3 src/makeGeometry.py --check --gmsh /path/to/gmsh        # nacelle, pinned gmsh
python3 src/makeGeometry.py --component phaseVI_blade           # regenerate the blade
python3 src/makeGeometry.py --gmsh /path/to/gmsh                # regenerate the nacelle
```

After regenerating, update the sha256 values printed by the generator in this
file — `--check` fails until they match.

## Known quirk

The OpenCASCADE boolean union of the hemisphere and the cylinder reports
`BOPAlgo_AlertUnableToOrientTheShape` from gmsh: the two primitives are tangent
along their shared circle, which is a degenerate case for the fuse. The union
still yields the correct single solid — the committed boundary mesh is closed,
manifold and outward-oriented.
