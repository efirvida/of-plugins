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
