// SPDX-License-Identifier: GPL-3.0-or-later
//
// nacelle.geo -- nacelle/hub surface for the nacelle actuator surface model.
//
// Geometry (Yang & Sotiropoulos, "A new class of actuator surface models for
// wind turbines", arXiv:1702.02108v4, Sec. 4.1 and Fig. 3): a hemisphere of
// radius R facing upstream, continued downstream by a circular cylinder of
// radius R and length L = 6R.  Fig. 3 dimensions the body as R (hemisphere),
// 6R (cylinder) and 2R (diameter); a pixel cross-check of the figure's
// dimension arrows reproduces 6.03 for the 6R/R ratio.
//
// Frame: the downstream end of the body sits at the streamwise origin x = 0
// and the nose points upstream, so the body spans x in [-(L + R), 0],
// matching the paper's periodic-nacelle case setup.  Units are metres; the
// canonical S1 parameters are R = 1 m (the validation case's Re = 1000 is
// based on R and the freestream velocity).
//
// Meshing is driven by makeGeometry.py, which runs
//     gmsh -2 -format stl -o stl/nacelle.stl src/nacelle.geo
// with the pinned gmsh version recorded in PROVENANCE.md.  The values marked
// "primary:" below are parsed by makeGeometry.py into the metadata parameter
// block and the input sha256; the derived values are kept in sync with them.
SetFactory("OpenCASCADE");

R = 1.0;                 // primary: hemisphere/cylinder radius [m]
cylinderRatio = 6.0;     // primary: cylinder length / R, Fig. 3 [-]
meshSizeRatio = 0.2125;  // primary: characteristic mesh length / R [-]

L = cylinderRatio*R;     // cylinder length [m]
lc = meshSizeRatio*R;    // characteristic mesh length [m]

// Hemisphere + cylinder.  The two primitives are tangent along their shared
// circle (x = -L), which is a degenerate case for the OCC boolean fuse: gmsh
// reports a benign `BOPAlgo_AlertUnableToOrientTheShape` warning and still
// produces the single solid whose boundary is the closed nacelle surface.
Sphere(1) = {-L, 0, 0, R};
Cylinder(2) = {-L, 0, 0, L, 0, 0, R};
BooleanUnion(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.Algorithm = 6;   // Frontal-Delaunay: uniform triangles, deterministic
Mesh.Binary = 1;      // binary STL
