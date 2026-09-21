# Reference data provenance — wall-resolved-LES profiles

Digitized reference profiles for the periodic-nacelle validation case. Every
committed point comes from a marker actually drawn in the published figure; no
point is interpolated, smoothed, or invented, and the detector/calibration
cross-checks below are the ones that were run.

## Source

| Item | Value |
|---|---|
| Paper | Yang & Sotiropoulos, *A new class of actuator surface models for wind turbines*, [arXiv:1702.02108v4](https://arxiv.org/abs/1702.02108v4) (2018) |
| Section | Sec. 4.1, "Flow over periodically placed nacelles" |
| Velocity figure | Fig. 5 (`u_nacelle`) — "Vertical profiles of time-averaged streamwise velocity ⟨u⟩ at different downstream locations for the nacelle case" |
| TKE figure | Fig. 6 (`k_nacelle`) — "Vertical profiles of turbulence kinetic energy k at different downstream locations for the nacelle case" |
| Reference series | **Symbols (open circles): wall-resolved LES** — this is the digitized series |
| Other series in the figures | solid = actuator surface, coarse grid; dash-dot = actuator surface, medium grid; dashed = permeable disk, coarse grid (not digitized) |
| Panels | x = 1R, 3R, 5R, 7R, 9R, 11R, 13R, 15R, 17R, 19R (odd multiples of R; the paper has **no 10R panel**) |
| Axes | ⟨u⟩/U and k/U² vs z/D, with D = 2R (the figure labels `z/D`) |

The reference is the paper's **502×348×348 wall-resolved LES**. Its dynamic
subgrid-scale model (Eq. 25) is **not available in standard OpenFOAM**; the
validation case therefore runs the WALE adaptation (see the package `README.md`
and `config/case.yaml`). The reference is not reproduced by this package — it is
read off the published figure only.

## Upstream files

The figures ship with the arXiv LaTeX source as raster EPS converted to PDF:

| File | sha256 |
|---|---|
| `u_nacelle-eps-converted-to.pdf` (Fig. 5) | `363cb1a4861ba33a67e69dffd65b2671ae00d003de587154d6989e3480901e82` |
| `k_nacelle-eps-converted-to.pdf` (Fig. 6) | `f3282e66afd188a324fa3ef843f0e1cbd991a14661d9820e0934c66d20fc85fd` |
| embedded image `u_nacelle.png` (2966×1784, 8-bit grey) | `cac1f1347d0467ebd64d90f10a9cdc9e314bd6f655b7e81ac972ddfd181de1cc` |
| embedded image `k_nacelle.png` (2996×1796, 8-bit grey) | `3e3e094dc85a1ea8f685c46b109794f897e22f8e76b99b2d0ea4065862e5a2fa` |

Obtained from the arXiv v4 source package (`https://arxiv.org/e-print/1702.02108v4`).
The compiled paper (`1702.02108v4.pdf`, page 8/9) embeds the **same** PNG bytes
(verified: identical sha256), so there is no independent second pixel source —
the cross-check below is between two independent *digitizations* of the same
published figure, not between two upstream images.

## Method

`tools/digitize_reference.py` (requires numpy, Pillow, scipy; PyMuPDF for
`--pdf`) is the digitizer. In the figure image:

1. **Panels** are located from the frame lines (2×5 grid).
2. **Calibration** per panel, twice and independently:
   - `tick`: the labelled axis values are matched to their major tick marks via
     the x-axis label centroids, and a least-squares line is fitted through the
     tick positions (y: a fit through the three labelled major ticks);
   - `label`: a least-squares line through the x-axis label centroids alone.
3. **Markers** (open circles, outer diameter ≈ 29 px, inner ≈ 19 px) are
   detected twice and independently:
   - `hole`: the marker interior is an enclosed white region, so enclosed holes
     of marker size are labelled and their centroids taken (164 detections in
     Fig. 5, 280 in Fig. 6);
   - `filter`: a ring matched filter — dark fraction on the annulus plus half
     the white fraction of the interior — with local maxima above 1.05
     (278 detections in Fig. 5, 299 in Fig. 6). The filter also finds markers
     crossed by a curve, where the interior is split and the hole detector is
     blind.
4. **Merge**: a clean-interior hole within 4 px of a filter detection confirms
   it (the committed centre is the mean of the two); unmatched detections are
   kept with their detector recorded in `digitization.json`.
5. **Isolated detections** (no neighbour within 50 px along the marker chain)
   are dropped: in these figures the italic station labels contribute such
   detections (the bowl of "R", the counter of "9"). 12 per figure are recorded
   with their pixel positions in `digitization.json`.
6. Pixel → data conversion uses the `tick` calibration;
   `z_over_D = f(ypix)`, `u_over_U`/`k_over_U2` = `f(xpix)`.

### Cross-check (two independent digitizations)

| Metric | Fig. 5 (⟨u⟩) | Fig. 6 (k) |
|---|---|---|
| Detections: hole / filter | 164 / 278 | 280 / 299 |
| Matched pairs (≤ 4 px) | 158 | 274 |
| Median centre deviation | 0.82 px = **0.0026 u/U** | 0.82 px = **7.5e-5 k/U²** |
| 95th percentile | 1.40 px = 0.0045 u/U | 1.24 px = 1.1e-4 k/U² |
| Maximum | 1.99 px = 0.0064 u/U | 1.65 px = 1.5e-4 k/U² |

**Stated tolerance:** at every point where both detectors fired, the two
independent digitizations agree within **2 px** (≈ 0.007 u/U and ≈ 1.5e-4
k/U²).

### Calibration cross-check (tick vs label fit)

Evaluated across the panel width (tick minus label calibration):

| Figure | panel centre | panel left edge | panel right edge |
|---|---|---|---|
| Fig. 5 (u/U) | −0.0014 | +0.0056 | −0.0083 |
| Fig. 6 (k/U²) | −7e-5 | −4.4e-4 | +2.9e-4 |

The label fit is the weaker calibration (a label centroid is not exactly the
tick position; the "1" label is up to 2.5 px off), which is why the tick fit is
used for the committed values and the label fit is the cross-check.

### Pixel resolution

| Quantity | Figure | Per pixel |
|---|---|---|
| ⟨u⟩/U | Fig. 5 | 0.00323 (310.0 px per unit) |
| k/U² | Fig. 6 | 9.09e-5 (11 000 px per unit) |
| z/D | both | 0.00407 (246.0 px per unit) |

Digitization uncertainty is therefore ≈ ±0.003–0.007 in u/U, ≈ ±1e-4 in k/U²
and ≈ ±0.004 in z/D — larger than the graphical line thickness but small
compared with the validation bands.

## Files

One CSV per station and quantity, plain two-column tables with a header:

```
data/reference/u_<station>.csv   z_over_D, u_over_U
data/reference/k_<station>.csv   z_over_D, k_over_U2
```

Stations: `1R`, `3R`, `5R`, `7R`, `9R`, `11R`, `13R`, `15R`, `17R`, `19R`
(20 files, 565 points in total). Points are sorted by `z_over_D`, which spans
at least ±1.32 at every station (each marker chain covers almost the full panel
height; exact per-station ranges are in `digitization.json`).

`data/reference/digitization.json` is the machine-readable audit trail:
per-figure input hashes, detector counts, the agreement statistics above, the
calibration fits and residuals, the dropped false positives, and every
committed point with its detector (`hole`, `filter`) and the sha256 of each
written CSV.

`data/reference/PROVENANCE.md` is this file.

## Reproducing

```sh
# 1. fetch the arXiv source (or take the two figure PDFs from it)
curl -L -o source.tar.gz https://arxiv.org/e-print/1702.02108v4
tar xzf source.tar.gz -C src/

# 2. digitize (writes the CSVs and digitization.json into data/reference/)
python3 tools/digitize_reference.py \
    --pdf u=src/u_nacelle-eps-converted-to.pdf \
    --pdf k=src/k_nacelle-eps-converted-to.pdf \
    --out data/reference
```

The digitizer is deterministic: re-running it on the same inputs reproduces
byte-identical CSVs and an identical report (verified). The reference CSVs are
committed, so the compare tool and its tests never need the figure or the
digitizer dependencies.

## Limitations

- The reference is a **digitized figure**, not raw data: values inherit the
  figure's line thickness, marker size and the pixel grid above.
- Only the wall-resolved-LES marker series is digitized; the actuator-surface
  and permeable-disk curves in the same panels are not, because the case is
  compared against the paper's reference simulation, not against the paper's
  own actuator-model results.
- The 502×348×348 reference LES is not reproduced here; a coarse/medium-grid
  actuator-surface run cannot be expected to reproduce it at the same fidelity
  (the paper reports exactly that distinction at 3R–7R).
- The paper's panels are at odd multiples of R; a run profile at 10R has no
  exact reference panel. `compareNacelle.py` records the station offset it uses
  when mapping the far-wake stretch onto the nearest digitized station.
