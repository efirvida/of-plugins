# Provenance — S809 airfoil coordinates

## `s809_somers_nlr.csv`

| Field | Value |
|---|---|
| Source | Somers, D. M., *Design and Experimental Results for the S809 Airfoil*, NREL/SR-440-6918 (January 1997) |
| DOI | https://doi.org/10.2172/437668 |
| Local verified copy | `phasevi_research/verified/s809_somers_nlr.pdf` (not committed) |
| Source PDF sha256 | `a7496aad5680d9d4001e36cf8b13d9603ea31956093cded60976dd8fcaffeeb3` |
| Table | Table 2, "S809 Airfoil Coordinates" (PDF page index 19, printed page 15) |
| Committed CSV sha256 | `861f3fe5bbda4d5e5d8985e47c5805bada76fc8e1b2157b944ff78d500f4af8b` |
| Row selection | All 31 upper-surface and 30 lower-surface coordinate pairs, upper LE→TE then lower LE→TE |
| Extraction date | 2026-09-23 |
| Columns | `surface`, `x_over_c`, `z_over_c` |
| Units | fraction of chord (dimensionless) |
| Cross-check sources | NREL/TP-500-29955 Table A-2; NREL/TP-442-7817 Table A1 |

The committed coordinates are the data the blade geometry generator consumes
(`geometry/src/blade_phasevi.py`); no second copy is committed. The local,
unverified text dump `phasevi_research/verified/s809_ramsay.txt` is **not** the
committed source and was not used as one — every committed and cross-check
number was extracted from the verified PDFs with PyMuPDF.

### Extraction recipe (repeatable)

- Tool: PyMuPDF 1.27.2.3 (MuPDF 1.27.2); `pdftotext` is not available on the
  extraction host.
- Source: `s809_somers_nlr.pdf` (sha256 above), page index 19, `get_text()`.
- Tokens: `-?(?:\d+\.\d+|\.\d+)` (the integral page number is not matched).
- Row grouping: the token stream is read as rows of
  `(upper x/c, upper z/c, lower x/c, lower z/c)`. The table prints 30 such
  rows plus one final trailing pair `(1.0000, 0.0000)` — the shared trailing
  edge — which the PDF word coordinates place in the upper-surface columns and
  which completes the upper surface: 31 upper points and 30 lower points, both
  ending at the shared trailing edge.
- Committed form: values at the published precision (5 decimals), the maximum
  printed precision of the source table.

### Closure convention

The data contain no explicit leading-edge coordinate: the first upper point is
`(0.00037, 0.00275)` and the first lower point is `(0.00140, -0.00498)`. The
closed profile of `blade_phasevi.py` therefore closes across this small blunt
leading-edge gap between those two points and carries the shared trailing edge
once. No coordinate is invented and no value is smoothed.

### Cross-check 1 — NREL/TP-500-29955 Table A-2 (independent transcription)

| Field | Value |
|---|---|
| Source | NREL/TP-500-29955, Table A-2 "Airfoil profile coordinates" (PDF page index 74, printed page 65) |
| Local verified copy | `phasevi_research/verified/29955_nlr.pdf` (not committed) |
| Source PDF sha256 | `822ee980e97780599ad4ec1dc8b4cc8ad41131b0cbbdffa19999d571f6ce4958` |
| Method | Same token regex and row grouping, PyMuPDF text layer |
| Result | **All published digits agree**: 31/31 upper and 30/30 lower coordinate pairs are numerically equal (the two reports print the same numbers with different leading/trailing zero styles — `.00575` vs `0.00575`, `1.0000` vs `1.00000`). |

Documented source erratum: TP-500-29955 Table A-2 prints one extra trailing row
in its lower column, `0.00000, 0.00000`, which is geometrically impossible (the
lower surface already reaches the shared trailing edge `(1.0000, 0.0000)` in
the preceding row; a point at the origin is not on the profile). The gate
compares the 30 published lower points and requires that extra row to be
exactly this documented zero row; any other disagreement aborts the extraction
with a non-zero exit.

### Cross-check 2 — NREL/TP-442-7817 Table A1 (independent measurement)

| Field | Value |
|---|---|
| Source | Ramsay, R. R., Hoffmann, M. J., Gregorek, G. M., *Effects of Grit Roughness and Pitch Oscillations on the S809 Airfoil*, NREL/TP-442-7817 (December 1995; no DOI recorded in the report), Table A1 "S809 Measured Model Coordinates, 18-inch Desired Cord" (PDF page indices 34–40, printed pages A-3…A-9) |
| Local verified copy | `phasevi_research/verified/s809_ramsay.pdf` (not committed) |
| Source PDF sha256 | `9b5d151f7361ce11aef1bfce25447d1f7e4f34c961e1266064310945e42cf825` |
| Rows | 223 upper and 221 lower measured points (inches), classified into the four table columns by PDF word x position |
| Normalization | inch ordinates divided by the 18-inch desired chord |
| Gate | 1 % chord tolerance (`0.01` in `z/c`) at every design station with `x/c ≤ 0.9`; the last 10 % of chord is excluded because Ramsay documents a trailing-edge fabrication thickening of 1.25 mm (0.05 in) added to the upper surface over that region, which is therefore **not** a design-coordinate disagreement |
| Measured maximum deviation, gated | upper `0.000894` `z/c` at `x/c = 0.00037` (25 design stations); lower `0.000906` `z/c` at `x/c = 0.04223` (24 design stations) — both well inside 1 % chord |
| Measured maximum deviation, last 10 % (reported, not gated) | upper `0.002725` `z/c` at `x/c = 0.98528` (consistent with the documented 0.05 in / 18 in = 0.0028 thickening addition); lower `0.001078` `z/c` at `x/c = 0.98446` |

The measured ordinate at each design station is obtained by linear
interpolation of the monotone measured series (asserted monotone in chord
station within `x ≤ 18 in`); the table's last rows lie past the measured
trailing edge (`18.0241…18.0274 in`) and are recorded as a source feature but
are outside the design table and the gated region. A deviation above 1 %
chord, a non-monotone measured series, or a missing expected row aborts the
extraction with a non-zero exit.

## Attribution

- Adapted from `mttbrbr/single-actuator-line` (main `8284be8c...`),
  GPL-3.0-or-later. All NREL-derived numbers are rebuilt from the public
  sources listed above.
- NREL citations: NREL/SR-440-6918 (DOI 10.2172/437668) for the S809 design
  coordinates, NREL/TP-500-29955 (DOI 10.2172/15000240) and NREL/TP-442-7817
  for the cross-checks.
- Pinned `of-plugins` commit for this validation package: `52fd258...`
  (actuator surface model patch included in the vendored `turbinesFoam` fork).

No report PDF is committed in this repository.
