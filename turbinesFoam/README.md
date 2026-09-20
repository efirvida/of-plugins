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
> model (ASM)** element (see "Actuator surface model" below).  The actuator
> line model (ALM) is unchanged and remains the default; the ASM is an
> additive, opt-in extension that does not exist upstream.

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
