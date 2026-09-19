# Provenance — experimental data

## Source

| Field | Value |
|---|---|
| Dataset | NWTC Unsteady Aerodynamics Experiment Phase VI statistics (WDH) |
| DOI | https://doi.org/10.21947/WDH-DAP/1910052 |
| Local verified copy | `phasevi_research/verified/wt_loads_statistics.xls` (17 MB, **not committed**) |
| Source sha256 | `8e1cfad21574aca8963b1eb73fc8c2e0c620006e81e6ee6ad733980da861d5f9` |
| Sheet | `ldsmean` (means of each 10-minute record) |
| Extraction script | `scripts/extractExperimentalData.py` (requires `xlrd`) |
| Extraction date | 2026-09-19 |
| Columns (0-based) | `File Name` 0, `Windspeed` 2, `Yaw` 3, `Repetition` 4, `WTBARO` 30, `WTATEMP` 34, `YAW` 47, `RPM` 52, `EAEROTH` 54, `LSSTQCOR` 58, `ROTPOW` 61, `CN30..CM95` 63..81 |

## Row selection

One canonical row per wind speed, 0° yaw only:

| Sequence | U∞ (m/s) | File Name | Sheet row (1-based) | Repetition | YAW (deg) |
|---|---|---|---|---|---|
| H | 7 | `h0700000` | 1204 | 0 | +0.0349 |
| H | 10 | `h1000000` | 1239 | 0 | −0.0487 |
| H | 13 | `h1300000` | 1272 | 0 | −0.0490 |
| H | 15 | `h1500000` | 1284 | 0 | −0.1290 |
| H | 20 | `h20m0000` | 1328 | 0 | −0.0540 |
| H | 25 | `h2500000` | 1342 | 0 | −0.0387 |
| S | 7 | `s0700000` | 1843 | 0 | +0.0110 |

Rule: keep only rows with `|YAW| ≤ 0.5°` and the lowest `Repetition`; the
workbook also contains yaw sweeps (2–180°) in the same sheet, which are
excluded. `m000` in `h20m0000` is the workbook's 0°-yaw token variant.

**Recorded correction:** the exploration table's 20 m/s values
(72.064 rpm / 9.219 kW / 1221.49 Nm / 3094.78 N) are the Sequence G teetered
row `g2000002`, not a Sequence H row. The canonical Sequence H 0° row is
`h20m0000`: 72.0839 rpm, ROTPOW 9.3661 kW, LSSTQCOR 1240.7201 Nm, EAEROTH
3152.3708 N, which yields TSR = 1.898 as required.

## Files and units

| File | Content | Units |
|---|---|---|
| `sequence_H_performance.csv` | 6 canonical Sequence H rows | `wind_speed_m_s` m/s, `rpm` rpm, `tsr` -, `rotpow_kw` kW, `lsstqcor_nm` Nm, `eaeroth_n` N, `wtbaro_pa` Pa, `wtatemp_degc` °C, `rho_kg_m3` kg/m³ |
| `sequence_H_spanwise.csv` | 6 rows × 5 stations (30/47/63/80/95 % span) | `r_over_R` fraction, `cn`/`ct`/`cm` non-dimensional |
| `sequence_S_performance.csv` | 1 Sequence S row (`s0700000`) | as above |
| `sequence_S_spanwise.csv` | 5 spanwise stations | as above |

Derived quantities:

- `tsr = rpm · 2πR / (60 · U∞)` with `R = 5.029 m` (the workbook has no TSR
  column); the values reproduce the frozen per-speed measured TSR table within
  ±0.002.
- `rho_kg_m3 = WTBARO / (287.058 · (WTATEMP + 273.15))`, the ideal-gas density
  from the row's own barometric pressure and air temperature (design decision
  D10; the workbook's `RHOTUN` column is not used).

Verified anchors (spec requirement): `h0700000` = 5.9458 kW / 789.8572 Nm /
1132.0692 N; `h2500000` = 11.9507 kW / 1580.4082 Nm / 4028.6306 N.

The .xls workbook is not committed; the derived numeric CSVs above are the
committed artifacts and can be regenerated with:

```sh
python turbinesFoam/validation/phaseVI/scripts/extractExperimentalData.py \
    --workbook /path/to/wt_loads_statistics.xls
```

## Attribution

- Adapted from `mttbrbr/single-actuator-line` (main `8284be8c...`),
  GPL-3.0-or-later.
- NREL citations: TP-500-29955 (DOI 10.2172/15000240), TP-442-7817, and the
  WDH Phase VI statistics workbooks (DOI 10.21947/WDH-DAP/1910052).
- Pinned `of-plugins` commit: `52fd258...` (ASM patch included).
