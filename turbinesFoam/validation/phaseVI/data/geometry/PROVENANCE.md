# Provenance — blade geometry

## `phaseVI_blade.csv`

| Field | Value |
|---|---|
| Source | NREL/TP-500-29955, *Unsteady Aerodynamics Experiment Phase VI: Wind Tunnel Test Configurations and Available Data Campaigns* (Hand et al., 2001) |
| DOI | https://doi.org/10.2172/15000240 |
| Local verified copy | `phasevi_research/verified/29955_nlr.pdf` (not committed) |
| Source sha256 | `822ee980e97780599ad4ec1dc8b4cc8ad41131b0cbbdffa19999d571f6ce4958` |
| Table | A-1 (blade geometry, 10.058 m diameter rotor, standard tip) |
| Committed file sha256 | `e678d6c52df9c8562361d61c65ad373535ea2351448874213b706bc6824e1576` |
| Row selection | All 26 radial stations of the standard-tip blade, r = 0.5083 m → 5.029 m |
| Extraction date | 2026-09-19 |
| Columns | `radius_m`, `chord_m`, `twist_deg_report`, `chord_mount` |
| Units | m, m, deg (positive towards feather), fraction of chord |

The station list and values are a byte-faithful copy of the committed upstream
SAL scaffold file `data/turbines/phaseVI_blade.csv` (identical sha256), which
was itself verified against TP-500-29955 Table A-1. `chord_mount` encodes the
pitch/twist axis location: 0.50 chord for the root cylinder transition and 0.30
chord through the airfoil (the report's 30 %-chord pitch/twist axis).
`twist_deg_report` is the report sign convention; `tools/element_data.py`
converts it to the turbinesFoam mounting angle with `-(twist + pitch)`.

No report PDF or workbook is committed in this repository.

## Attribution

- Adapted from `mttbrbr/single-actuator-line` (main `8284be8c...`),
  GPL-3.0-or-later. All NREL-derived numbers are rebuilt from the public
  sources listed above.
- Pinned `of-plugins` commit for this validation package: `52fd258...`
  (actuator surface model patch included in the vendored `turbinesFoam` fork).
