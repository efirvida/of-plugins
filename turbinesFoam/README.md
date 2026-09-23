# turbinesFoam

[![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)
![OpenFOAM v2412](https://img.shields.io/badge/OpenFOAM-v2412-brightgreen.svg)
![OpenFOAM v2406](https://img.shields.io/badge/OpenFOAM-v2406-brightgreen.svg)
![OpenFOAM v2312](https://img.shields.io/badge/OpenFOAM-v2312-brightgreen.svg)
![OpenFOAM v2306](https://img.shields.io/badge/OpenFOAM-v2306-brightgreen.svg)
![OpenFOAM v2212](https://img.shields.io/badge/OpenFOAM-v2212-brightgreen.svg)

> **Vendored fork** — this directory is an independent copy of
> [turbinesFoam/turbinesFoam](https://github.com/turbinesFoam/turbinesFoam),
> included as a plugin base for research development.  The original git
> history was stripped and upstream changes are not tracked automatically.
> Licence: GPL-3.0 (see `LICENSE`).
>
> **Fork divergence:** this copy adds an optional **blade actuator surface
> model (ASM)** element (see "Actuator surface model" below) and an
> imported-surface **blade surface distributor** (see "Blade surface model"
> below).  The actuator line model (ALM) is unchanged and remains the default;
> both additions are opt-in extensions that do not exist upstream.

turbinesFoam is a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method.

[![](https://cloud.githubusercontent.com/assets/4604869/10141523/f2e3ad9a-65da-11e5-971c-b736abd30c3b.png)](https://www.youtube.com/watch?v=THZvV4R1vow)

Be sure to check out the
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).

## Installation

### Compile from source

```sh
cd $WM_PROJECT_USER_DIR
git clone https://github.com/efirvida/of-plugins.git
cd of-plugins/turbinesFoam
./Allwmake
```

or, from the repository root, `./Allwmake` builds this plugin together with
the other of-plugins libraries.

## Usage

See the tutorials located in the `tutorials` directory.

## Contributing

Pull requests are very welcome!
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details.

## Features

`fvOptions` classes for adding actuator lines and turbines constructed from
actuator lines to any compatible solver or turbulence model, e.g.,
`simpleFoam`, `pimpleFoam`, `interFoam`, etc.

## Actuator surface model (ASM)

An optional blade **actuator surface model** element
(`actuatorSurfaceElement`, Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 2.1)
can be selected per line/blade instead of the default actuator line element.
It computes the same blade element momentum loads as the line model (same
coefficient lookup, dynamic stall, added mass, end effects, and CSV output)
but differs in three ways:

1. **Chord-averaged inflow** — the inflow velocity is the arithmetic mean of
   the cell-interpolated velocity at the midpoint of each of `nChordwise`
   equal chord strips, instead of a single sample at the quarter chord.
2. **Uniform strip-wise force projection** — the element force is split
   evenly across the `nChordwise` strips (`forceVector_ / nChordwise` per
   strip), each projected with its own Gaussian kernel.
3. **Mesh-based projection width** — epsilon is `2 * cbrt(V) * meshFactor`
   (cell volume and `GaussianCoeffs.meshFactor` only), with no chord-length
   term, so the model couples to meshes finer than the chord.

### Configuration

In the line (`actuatorLineSourceCoeffs`) or blade subdictionary:

```
elementType actuatorSurfaceElement;
nChordwise 5;      // optional, default 5
```

Both keys are optional and reach every element of the line/blade. With
`elementType` absent (or set to `actuatorLineElement`) the case runs the
default, unchanged actuator line model. An unknown `elementType` fails with a
run-time selection error listing the registered types. See the
`tutorials/axialFlowTurbineASM` tutorial for a full HAWT example and
`tutorials/axialFlowTurbineASM/compareALMvsASM.py` for an ALM-vs-ASM
side-by-side comparison of two run directories. The comparison reads whatever
each run directory contains, so run both cases to the same `endTime` for a
meaningful comparison: the ASM tutorial ships with `endTime 0.1` to bound its
runtime, while the ALM tutorial defaults to `0.5`.

### Notes

- **Kernel adaptation:** the paper's 3D smoothed-cosine delta kernel
  (5-cell) is not implemented; the existing Gaussian projection kernel is
  reused and documented as the adaptation.
- **Force preservation:** the Gaussian kernel is evaluated at cell centres
  and not discretely renormalized; total-force preservation is approximate
  and inherited from the line model.
- **Turbulence injection:** only the momentum/force path uses the
  strip-wise spreading; turbulence injection (`addTurbulence`) keeps the
  inherited single-point kernel.
- **HAWT-only:** the surface element is turbine-agnostic in principle, but
  only the axial-flow HAWT path is tested here. CFTAL/VAWT usage is
  untested.
- **Epsilon semantics:** on the existing coarse tutorial mesh the ALM already
  selects the mesh-based projection width (`epsilon (mesh-based): 0.183` in
  both ALM and ASM debug output), so both models use the same epsilon there.
  The load shift in the tutorial comparison (mean Cp 0.67 for the ASM vs 0.55
  for the ALM over the same 0.1 s window) therefore comes from the
  chord-averaged inflow and the strip-wise spreading, not from epsilon. The
  surface model's mesh-only epsilon is intended to be re-validated on a finer
  mesh where it gives 1–2 cells of overlap.

## Nacelle surface model

An optional nacelle/hub **actuator surface model** (`nacelleSurfaceSource`,
Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 2.2) reads a triangulated surface
from an STL file and applies the paper's direct-forcing normal force (Eq. 19)
plus a friction-based tangential force (Eq. 21, Schultz–Grunow `cf` by default,
Eq. 22), spread onto the background grid with the smoothed four-point cosine
kernel (Eqs. 7, 8, 18). It carries no blade-element data and is usable
standalone (a rotor-less `fvOptions` entry — the paper's periodic-nacelle case)
or as an owned child of an axial-flow turbine via a `nacelle {}` subdictionary.

### Configuration

Standalone entry:

```
nacelle
{
    type            nacelleSurfaceSource;
    active          true;

    nacelleSurfaceSourceCoeffs
    {
        fieldNames          (U);
        selectionMode       cellSet;
        cellSet             nacelleCells;

        geometry            "geometry/nacelle.stl";  // required, ASCII/binary STL
        referenceVelocity   1.0;                     // U for Eqs. 21, 22 [m/s]
        rho                 1.0;                     // reference density for SI forces/CSV
        nu                  -1.0;                    // >= 0 uses it; else transportProperties.nu
        cfModel             schultzGrunow;           // or constant
        cf                  -1.0;                    // constant override / calibration
        bodyOrigin          (0 0 0);                 // sampling-contract body frame
        bodyAxis            (1 0 0);
        referenceArea       1.0;                     // cd = |F|/(0.5*rho*U^2*referenceArea)
        writeForceField     true;                    // write force.<name>
        writePerf           true;                    // postProcessing/nacelle/<name>.csv
        writeNodePerf       false;                   // postProcessing/nacelle/<name>_nodes.csv
    }
}
```

Composed through `axialFlowTurbineALSource`, which builds the source in
`createNacelle()` from the `nacelle {}` subdictionary (inheriting
`fieldNames`/`selectionMode`/`cellSet` and using `mag(freeStreamVelocity)` as
`referenceVelocity`):

```
nacelle
{
    geometry            "geometry/nacelle.stl";
    cfModel             schultzGrunow;
    // ... any nacelleSurfaceSourceCoeffs key except the inherited ones
}
```

Output: `postProcessing/nacelle/<name>.csv` with
`time,fx,fy,fz,f_n_mag,f_tau_mag,cd` (total force on the body in newtons,
`cd` from `referenceArea`) and, with `writeNodePerf true`,
`postProcessing/nacelle/<name>_nodes.csv` with the per-node body-frame
positions/normals/forces (`time,node,x,y,z,nx,ny,nz,fx,fy,fz,area`). The
distributed field is registered as `force.<name>` and is force per unit volume
per unit density, matching the incompressible `fvOptions` convention; the
compressible overload weights it by the local density.

### Notes

- **Friction-model limits:** Schultz–Grunow assumes a zero-pressure-gradient
  turbulent boundary layer and is keyed by the streamwise distance behind the
  nose (`Rex = U·x/ν`). In the hemisphere nose region `Rex → 0`, the relation
  is invalid, so `cf` is set to zero there and the nose friction is
  under-predicted by construction. Use a constant `cf` override for
  calibration studies.
- **MPI:** every rank holds the full node list and applies each node force to
  its own cells within the kernel support (a 5³-cell stencil). A node whose
  containing cell sits on another rank is resolved with a `findCell` and
  reduce/minimum sentinel; a sample that is unreachable on every rank is a
  fatal error.
- **`h` and `ũ`:** `h = cbrt(V_cell)` of the containing cell, and `ũ` is the
  current `eqn.psi()` iterand (PIMPLE does not expose the paper's predicted
  velocity) — both are documented adaptations.
- **Tests:** `tests/test_nacelle.py` runs the standalone case in serial and on
  two ranks (`tests/nacelleSurface`, a two-triangle plate with an analytic
  force), checks the CSV and total-force preservation, and asserts the fatal
  error paths for missing, empty and unparseable STL files.
- **Fork divergence:** `nacelleSurfaceSource`/`nacelleSurfaceSampler` and the
  `createNacelle()` wiring are new in this fork and are not part of upstream
  turbinesFoam. The nacelle is static — `rotate`, `tilt` and `yaw` move the
  blades and hub only.

### Validation

The fork adds a paper-faithful validation package for the nacelle model at
`validation/nacelle-asn/` (Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 4.1):
a rotor-less hemisphere + cylinder nacelle at Re = 1000 in a streamwise-periodic
30R x 20R x 20R cell, represented only by `nacelleSurfaceSource`.

- `config/case.yaml` is the single source of truth; `tools/generate_case.py`
  renders the committed `case/` skeleton (non-destructive `--check`) for the
  paper's coarse (153x80x80) and medium (115x151x151) grids, and the case
  probes the metric-station lines every time step for the profile series.
- Digitized wall-resolved-LES reference profiles (Fig. 5/Fig. 6, panels
  1R..19R) live in `data/reference/` with `PROVENANCE.md`;
  `scripts/compareNacelle.py` forms `⟨u⟩(z)`, the resolved `k(z)` and
  `CD = |F_drag| / (0.5 ρ U∞² πR²)` and enforces the per-station/per-grid
  acceptance. The paper's permeable-disk `CD = 0.48` is a datum, not a target.
- **Geometry:** `geometry/` is the shared deterministic gmsh + Python pipeline
  (`makeGeometry.py --check`) that produces the nacelle STL consumed by the
  case; blade STLs are reserved for a later change (the layout and metadata
  schema accommodate them).
- Staged Slurm: `scripts/slurm/stage0.slurm` (development partition: mesh +
  short stability run) is the only authorized execution; `production.slurm`
  (long queue: full wash-out + averaging) is **prepared only** and refuses to
  run without an explicit authorization.
- Closure adaptation: the paper's dynamic SGS model is not in standard
  OpenFOAM, so the headline is LES **WALE** with a documented URANS k-ω SST
  fallback; the nacelle has no body-fitted mesh, so near-wall quantities are
  not claimed.

## Blade surface model (imported-surface distributor)

An optional per-blade **blade actuator surface** (`bladeSurfaceSource`) imports
a triangulated blade surface (STL) and distributes the blade element momentum
loads over its triangles instead of the element strip projection (Yang &
Sotiropoulos, arXiv:1702.02108v4, Sec. 2.1). The owning `actuatorLineSource`
constructs it when the blade subdictionary carries `surfaceGeometry`; it is a
**distribution-only** model: the BEM chain (coefficient lookup, dynamic stall,
added mass, end effects, chord-averaged inflow) is unchanged and the surface
never samples per-node inflow or recomputes BEM loads.

### Configuration

In a blade subdictionary (`blades { blade1 { ... } }`) or a standalone
`actuatorLineSource`:

```
surfaceGeometry       "constant/triSurface/phaseVI_blade.stl"; // required; presence activates
surfaceOrigin         (0 0 0);        // construction frame; injected by AFTAL
surfaceSpanDirection  (0 0 1);        // outward root -> tip; injected by AFTAL
surfaceChordDirection (0 -1 0);       // trailing -> leading; injected by AFTAL
kernel                cosine;         // cosine (default) | gaussian (ablation)
meshFactor            1.0;            // gaussian only; fallback below
rho                   1.0;            // reference density for the SI contract/CSV
referenceVelocity     1.0;            // accepted for shared-template parity; unused
nu                    -1.0;           // accepted for shared-template parity; unused
bodyOrigin            (0 0 0);        // optional body-frame override
bodyAxis              (0 0 1);        // optional body-frame override
writePerf             true;           // per-station CSV
writeNodePerf         false;          // opt-in per-node CSV
logDistribution       true;           // per-addSup instrumentation line/CSV
```

`surfaceGeometry` accepts an ASCII or binary STL; it is read as given, or
resolved against the case directory when that is not an existing file (the S1
resolution). `axialFlowTurbineALSource` injects `surfaceOrigin`,
`surfaceSpanDirection` (outward root -> tip) and `surfaceChordDirection`
(trailing -> leading, before element pitch) from the blade construction frame
after cone/azimuth; a standalone source supplies them itself, and a user
`bodyOrigin`/`bodyAxis` overrides the body-frame defaults.

**`projectElementForce` (suppression semantics):** when the surface is active,
the blade source injects `projectElementForce false` into every element dict
of that blade, so the element strip projection is suppressed and the node
distribution reaches the momentum equation exactly once. The element still
computes its force and writes its element CSV and public `force()`; only the
strip projection into the field is skipped. Without `surfaceGeometry` the key
is not injected (element default `true`) and the delivered ALM/no-mesh ASM
behaviour is byte-identical. A user-set `projectElementForce` in the blade
subdictionary is passed through when no surface is configured.

**Kernel modes:** `kernel cosine` (default) is the paper's smoothed four-point
cosine kernel with support `2.5*h` (Eq. 8, `h = cbrt(V)` of the node's
containing cell). `kernel gaussian` selects the ablation that matches the
no-mesh ASM width: `eps = 2*cbrt(V)*meshFactor` with the truncated support
`eps*sqrt(ln 1000)` (D7). `meshFactor` is read from the surface subdictionary,
then from `profileData GaussianCoeffs.meshFactor` in that dictionary or,
failing that, in the owning elements (where AFTAL copies the blade
`profileData`), then the element default `2.0`. The Gaussian mode is an
**ablation** that isolates the
distribution geometry from the kernel/width confound, not a separate model.

### Output

Written on the master rank under `postProcessing/bladeSurface/`:

- `<owner>.surface.csv` (`writePerf`, default `true`), one row per element
  (its radial patch):
  `time,station,root_dist,area,force_x,force_y,force_z,c_ref_n,c_ref_t,f_ref_n,f_ref_t`.
  `force_*` is the patch force on the blade in newtons (`rho`-scaled);
  `root_dist` uses the element convention and `c_ref_*`/`f_ref_*` the element
  definitions, so `comparePhaseVI.py` consumes the file without a new
  conversion.
- `<owner>.surface_nodes.csv` (`writeNodePerf`, default `false`):
  `time,node,x,y,z,nx,ny,nz,fx,fy,fz,area,station,chord_fraction`, per node in
  the blade body frame with the force on the blade in SI units.
- `<owner>.surface_distribution.csv` (`logDistribution`, default `true`):
  `time,nodes,candidates,mean_candidates,max_candidates,seconds`, one row per
  `distribute()` call. `candidates` is the number of candidate cell entries
  visited in that call (the bounded-query evidence: it is far below the naive
  `nodes * N_local cells` scan), `mean_candidates`/`max_candidates` are the
  per-node mean and maximum, and `seconds` is the wall time of the call. The
  same counts are emitted as one `Info` line on the master rank
  (`Blade surface distribution '<owner>.surface': ...`). The counters are
  observational only and never change the distributed result.

Both files are written with 12 significant digits: the per-node file is the
audit source of the moment convention below and its reconstruction oracle
(`rtol 1e-6`) cannot be represented with the default six digits. The element
and turbine CSVs keep the delivered precision.

**Moment/torque convention (D8):** while the surface is active,
`actuatorLineSource::moment()` returns the distributed node moment
`sum_nodes X_i x F_i` (global frame, per unit density) **instead of** the
element moment, and the axial-flow turbine projects it on the rotor axis for
the torque, `ct` and `cp`. The reported torque therefore reflects the
application points of the load actually applied; adding the element moment
would double count the same force. The per-node CSV is the audit source for
the reported torque.

### Notes

- **Rotating frame:** `rotate`, `translate` and both `pitch` overloads are
  forwarded by the blade source, so the nodes stay in lockstep with the
  elements; `positions()` is the blade body-frame contract (it reflects the
  azimuth) and the node ordering is fixed at construction. `setSpeed`,
  `scaleVelocity` and `setOmega` do not move the geometry.
- **MPI:** every rank holds the full canonical node list and builds candidate
  lists over its own cells; each local cell is written exactly once and the
  totals are `returnReduce`d. A node whose containing cell is on another rank
  is resolved with the S1 `findCell` + reduce/minimum sentinel; an unreachable
  sample is a fatal error.
- **Bounded query:** each node needs one `cellSize` lookup and one support
  query at construction; each `addSup` iterates only the cached candidate
  list, never every local cell.
- **Fork divergence:** `bladeSurfaceSource`, `bladeSurfaceSampler` and the
  shared `surfaceSamplerBase` are new in this fork and are not part of
  upstream turbinesFoam. The nacelle sampler is re-based on
  `surfaceSamplerBase` with unchanged (byte-identical) output.
- **Tests:** `tests/test_blade_surface.py` runs the standalone
  `tests/bladeSurface` fixture (partition of unity, default path, suppression,
  CSV schemas, two-rank totals, Gaussian conservation, STL failure paths) and
  the minimal rotating `tests/bladeSurfaceAFTAL` fixture (rotation lockstep
  and the D8 surface moment against `turbine.csv`).

## Publications

Bachant, P., Goude, A., and Wosnik, M. (2016) [_Actuator line modeling of vertical-axis turbines_](https://arxiv.org/abs/1605.01449). arXiv preprint 1605.01449.

## How to cite

The latest release of turbinesFoam can be cited via DOI thanks to Zenodo: [![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)

## Acknowledgements

This work was funded through a National Science Foundation CAREER award,
principal investigator Martin Wosnik ([NSF CBET
1150797](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1150797), Energy for
Sustainability, original program manager Geoffrey A. Prentice, current program
manager Gregory L. Rorrer).

OpenFOAM is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

Interpolation, Gaussian projection, and vector rotation functions adapted from
NREL's [SOWFA](https://github.com/NREL/SOWFA).
