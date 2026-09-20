# Geometry provenance — `turbinesFoam/geometry/`

Everything committed under `turbinesFoam/geometry/` is generated from the gmsh
sources in `src/` by `src/makeGeometry.py`. This file records where the geometry
comes from, what generated it, and the hashes that link the committed artefacts
to their inputs.

## Sources

| Component | Source | Status |
|---|---|---|
| `nacelle` | Yang & Sotiropoulos, *A new class of actuator surface models for wind turbines*, arXiv:[1702.02108v4](https://arxiv.org/abs/1702.02108v4) (2018), Sec. 4.1 ("Flow over periodically placed nacelles") and Fig. 3 (schematic of the nacelle geometry with the triangular surface mesh) | generated in S1 |
| `blade0`, `blade1`, `blade2` | MEXICO rotor geometry tables (the same paper, Sec. 4.2) | **deferred to S2** — the layout and the metadata schema reserve them; no blade STL is generated in S1 |

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

| Item | Value |
|---|---|
| gmsh | **4.15.2** (pip package `gmsh`, `libgmsh` 4.15) |
| mesher options | `Mesh.Algorithm = 6` (Frontal-Delaunay), `Mesh.CharacteristicLengthMin = Mesh.CharacteristicLengthMax = 0.2125·R`, `Mesh.Binary = 1` |
| CAD kernel | OpenCASCADE (`SetFactory("OpenCASCADE")`), hemisphere ∪ cylinder |
| geometry source | `src/nacelle.geo` (parameters `R = 1`, `cylinderRatio = 6`, `meshSizeRatio = 0.2125`) |
| generator script | `src/makeGeometry.py` |
| output | `stl/nacelle.stl` (binary STL), `metadata/nacelle.json` |

Regeneration is byte-identical (verified: two consecutive runs produced the same
sha256 for the STL and the metadata). `makeGeometry.py --check` re-runs the
generation into a temporary directory and byte-compares it against the committed
artefacts; it never writes to the tree.

gmsh's *element ordering* is not stable across environments: the same mesh
(same vertices, same triangles) is emitted in a different order depending on the
ambient environment, which would break byte-identical regeneration. The
generator therefore re-serializes the mesh before writing: every triangle is
rotated into its lexicographically smallest cyclic form (preserving the winding,
hence the outward normal) and the triangles are sorted by vertex coordinates.
The committed STL also carries a fixed header and recomputed facet normals, so
its bytes depend only on the mesh geometry — verified identical in the plain and
in the OpenFOAM (`WM_PROJECT_DIR/etc/bashrc`) environments, and a duplicate
triangle would be rejected.

### Triangle count

The committed surface has **2616 triangles**. The paper reports 2652 triangles
for its nacelle surface (Sec. 4.1) — a difference of **−1.4 %**, within the
`~2652` target. The count is set by the `meshSizeRatio` seed: with the pinned
gmsh the triangle count moves in plateaus as the seed varies (`0.2095..0.2145`
all give 2590–2616 triangles), so 2616 is the closest reproducible value; the
seed is deliberately not over-tuned to a single count.

Mesh statistics of the committed STL (all consistent with an *inscribed*
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

## Hashes

| Input / artefact | sha256 |
|---|---|
| input (`src/nacelle.geo` + parameter block) | `e0a69771736aa3b26a11b00880ddd0acf2cf45d8782f13d56fd78ccb89d20de5` |
| `stl/nacelle.stl` | `18c7beef9e323ba28da20f5c813232c0059693f221a594fd913c4e69c14cb74c` |
| `metadata/nacelle.json` | `f237fc40840e9d0970583e32f3b4f21427ab8c29d83755b8f6b55ebb37437b2d` |

`input_sha256` is computed over the raw `src/nacelle.geo` bytes, a `NUL` byte,
and the canonical parameter block
(`json.dumps(parameters, sort_keys=True, separators=(",", ":"))`), where
`parameters` is the `primary:` block declared in the `.geo`
(`{"R": 1.0, "cylinderRatio": 6.0, "meshSizeRatio": 0.2125}`). The recipe is
implemented by `makeGeometry.input_sha256()` and re-checked by the test suite.

## Environment requirement

The pip `gmsh` module links `libGLU.so.1`/`libGL.so.1` (and the MPI stack).
`import gmsh` fails — and therefore `src/makeGeometry.py` cannot regenerate —
unless the system library environment is present. On SDumont:

```sh
module load glu/9.0.2_gnu          # provides libGLU.so.1 for libgmsh
```

Without a runnable gmsh, `tests/test_nacelle_data.py` **skips** the regeneration
checks with an explicit reason and still runs the pure-Python checks on the
committed artefacts.

## Reproducing

```sh
cd turbinesFoam/geometry
python3 src/makeGeometry.py --check --gmsh /path/to/gmsh   # verify, writes nothing
python3 src/makeGeometry.py --gmsh /path/to/gmsh           # regenerate in place
```

After regenerating, update the two sha256 values printed by the generator in
this file — `--check` fails until they match.

## Known quirk

The OpenCASCADE boolean union of the hemisphere and the cylinder reports
`BOPAlgo_AlertUnableToOrientTheShape` from gmsh: the two primitives are tangent
along their shared circle, which is a degenerate case for the fuse. The union
still yields the correct single solid — the committed boundary mesh is closed,
manifold and outward-oriented.
