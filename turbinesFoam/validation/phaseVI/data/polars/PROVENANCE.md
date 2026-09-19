# Provenance — S809 polars

The committed polar data are rebuilt from the public NREL sources below. The
SAL scaffold polar `S809_Re1M_extended` is **not** used (it contains four
flipped `Cm` signs and unexplained `Cd` values, and a hand-made ±30..±180°
extension).

## Sources

| Source | DOI / URL | Local verified copy (not committed) | sha256 |
|---|---|---|---|
| NREL/TP-500-29955, Tables A-3..A-8 | https://doi.org/10.2172/15000240 | `phasevi_research/verified/29955_nlr.pdf` | `822ee980e97780599ad4ec1dc8b4cc8ad41131b0cbbdffa19999d571f6ce4958` |
| NREL/TP-442-7817, Tables B1..B4 (clean) | NREL/TP-442-7817 (December 1995; no DOI recorded in the report) | `phasevi_research/verified/s809_ramsay.pdf` | `9b5d151f7361ce11aef1bfce25447d1f7e4f34c961e1266064310945e42cf825` |

`raw/table_A3..A8.csv` transcribe the TP-500-29955 tables; the token-count
rule was cross-checked against the TP-442-7817 tables (which carry explicit
`--` placeholders) for the two OSU tables that appear in both reports
(A-6 = Table B1, A-7 = Table B2; all shared values agree numerically).
`raw/table_TP442-7817.csv` stacks the four clean (ungritted) OSU tables
B1..B4 at nominal Re = 0.75/1/1.25/1.5 × 10⁶.

### Raw tables

| File | Source table | Row selection | Columns (units) |
|---|---|---|---|
| `raw/table_A3.csv` | TP-500-29955 A-3 | CSU Re 0.3e6, all rows | `alpha_deg`, `cl`, `cdp` (-) |
| `raw/table_A4.csv` | TP-500-29955 A-4 | CSU Re 0.5e6, all rows | `alpha_deg`, `cl`, `cdp` (-) |
| `raw/table_A5.csv` | TP-500-29955 A-5 | CSU Re 0.65e6, all rows | `alpha_deg`, `cl`, `cdp` (-) |
| `raw/table_A6.csv` | TP-500-29955 A-6 | OSU Re 0.75e6, all rows | `alpha_deg`, `cl`, `cdp`, `cdw`, `cm`, `re_millions` |
| `raw/table_A7.csv` | TP-500-29955 A-7 | OSU Re 1e6, all rows | `alpha_deg`, `cl`, `cdp`, `cdw`, `cm`, `re_millions` |
| `raw/table_A8.csv` | TP-500-29955 A-8 | DUT Re 1e6, all rows | `alpha_deg`, `cl`, `cdw`, `cm` |
| `raw/table_TP442-7817.csv` | TP-442-7817 B1..B4 | OSU clean series, all rows | `re_millions`, `alpha_deg`, `cl`, `cdp`, `cdw`, `cm`, `re_measured_millions` |

Blank `cdw` cells are the source's blank wake-traverse entries (the
transcription preserves the blank, it does not invent a value). `cm` is the
quarter-chord moment coefficient about the 25 % chord axis; all coefficients
are non-dimensional and `alpha_deg` is in degrees. Extraction date:
2026-09-19.

## Derived files

| File | Built from | Rule |
|---|---|---|
| `S809_OSU_Re1M_total.dat` | `raw/table_A7.csv` (+ `raw/table_A3.csv` extension) | `Cd = Cdw` (total drag) where reported, else `Cdp` (documented fallback); `|alpha|` beyond the measured OSU range carries a clearly commented static extension from CSU Table A-3 (mirrored for negative alpha) plus a flat-plate closure at ±120/150/180° so the lookup never leaves the table. Format `(alpha Cl Cd Cm)` for the single-Re `profileData` table. |
| `S809_multiRe.dat` | `raw/table_TP442-7817.csv` | Same `Cdw`→`Cdp` fallback; linearly interpolated onto the shared Re = 1e6 angle grid restricted to the range every Reynolds level measures (no extrapolation). Emitted as `tableType multiRe; ReList (7.5e5 1e6 1.25e6 1.5e6); clData/cdData/cmData`. |

The CSU Tables A-3..A-5 and the DUT Table A-8 are committed raw data for
traceability and are **not** blended into the derived polars (cross-source
blending introduces kinks and is excluded by the design). The multi-Re file is
a sensitivity variant; the case generator renders the single-Re baseline.

Rebuild and staleness check:

```sh
python turbinesFoam/validation/phaseVI/scripts/buildPolars.py
python turbinesFoam/validation/phaseVI/scripts/buildPolars.py --check
```

## Attribution

- Adapted from `mttbrbr/single-actuator-line` (main `8284be8c...`),
  GPL-3.0-or-later.
- NREL citations: TP-500-29955 (DOI 10.2172/15000240), TP-442-7817.
- Pinned `of-plugins` commit: `52fd258...` (ASM patch included).
