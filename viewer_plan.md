# PlateDriver viewer — plan

Draft prepared with Claude Code.

A local Geode-style globe viewer for PlateDriver runs. It lives in `viewer/` in this repo, not in
GeodeViewers. It has two uses:
- a **teaching** tool: students explore one finished run;
- an **inspection** tool: look at any run's outputs.

## Decisions

| Topic | Decision |
|---|---|
| Purpose | Teaching and run inspection, in one viewer. |
| Route | **Hand-built**: Vite + three.js on Geode's `core/`, copied as PlateTree does. The generated route (`geode` package, petrify Explorer) has no grid or volume layers, so it can't show the chosen layers. |
| First-version layers | Plates and boundaries, surface grids, mantle volume. Point sets (earthquakes, seamounts, tracers, sutures) come later. |
| Comparison | A run picker (one globe; switching run keeps the time step). Linked side-by-side globes can come later. |
| Time | One slider labelled **model time t (Myr)**, running 0 → end, left to right. Every layer snaps to the nearest saved step. No Ma ages, so students don't read the run as real Earth. |
| Data flow | An export script, `viewer/prep/export_run.py <run_dir> [--name NAME]`, converts one run folder into the viewer's format and adds it to the run list. It works on any run folder. `run_SDP.py` is not changed. |
| Distribution | Local only: `npm run dev` on `localhost:5173`. A shareable build is decided later. |

## Geode conventions, and where this viewer departs from them

- **Kept:**
  - Layout: `prep/`, `public/data/`, `src/core/` (copied, not imported), `src/platedriver/`.
  - Runtime: node and the browser only. The export script runs in `pygmt17`, like PlateDriver.
  - Formats: volumes and grids use the Archive frame format (`ARCHIVE_FORMAT.md`), one uint8
    file per (variable, frame), decoded with `encode_min`/`encode_max`.
  - Every export is checked against its source (see Verification).
- **Different:**
  - The plates are not a published Reconstruction Model. Each run carries its own boundaries,
    plate-ID raster and rotations, and runs are never mixed. Geode's rule that a model is shown
    only under its own reconstruction holds automatically: a run is always shown with its own
    plates.
  - Time is forward model time, not Ma. The Ma readout in core, for example
    `tomography/instance.ts`, is replaced in the wrapper.
  - Boundary types: PlateDriver's `ridge` and `subduction` map to Geode's `ridge` and
    `subduction`. PlateDriver's `conservative` means "neither, by the divergence thresholds",
    not a checked transform, so it maps to `other`.

## Export format (per run)

```
viewer/public/data/runs.json                    list of exported runs (name, label, path)
viewer/public/data/runs/<name>/
  run.json          name, force model and modes, key params, saved steps and times, git commit,
                    layer list with encode ranges and colour ramps
  boundaries/NNNN.geojson                       boundary segments with type and plate IDs
  plates/NNNN.json                              per-plate pole, rate, RMS speed, area, continental
                                                fraction (from plate_stats.csv; velocity arrows)
  grids/<var>/NNNN.bin                          uint8, 0.5° lon/lat (see below): seafloor_age,
                                                topo_bathy, divergence, plate_id (categorical)
  volume/dVs/NNNN.bin                           uint8, the mantle dVs/Vs volume
  verification/                                 source-against-export figures and numbers
```

- **Grids** are resampled from 0.25° to 0.5° by default (`--grid-deg`; each file about 0.26 MB,
  with 51 steps × 4 variables ≈ 50 MB per run). This keeps exports small; 0.25° stays
  available.
- **Plate IDs** are exported as a categorical grid, so the plate under the cursor can be looked
  up.
- **The mantle volume** is saved on a 5° cell-centred grid (72 × 36 × 29 depths) with no pole
  rows. Geode's volume format needs both poles, so each pole row is set to the area-weighted
  mean of the ring next to it.
- **Spherical geometry:**
  - Grids are lon/lat, draped on the sphere by the core shaders.
  - Boundaries are great-circle segments. They are densified to at most 1° before rendering, so
    long segments don't cut through the globe.
  - Velocity arrows come from each plate's rotation, ω × r, evaluated on the sphere.

## Stages (one at a time, each with a gate)

| Stage | What | Gate |
|---|---|---|
| V0 | Scaffold `viewer/`: Vite + TypeScript, core copied from Geode's `viewer/src/core`, an empty globe with graticule, time slider over a fixed range. `node_modules/` and `public/data/runs/` gitignored | dev server shows the globe |
| V1 | `export_run.py` for boundaries and plates; the viewer shows boundaries (ridge red, trench blue with triangles on the overriding side, other black), plate IDs and velocity arrows, snapping to saved steps | screenshots at 0, 250 and 500 Myr beside the `diagnose_run` map sequence panels of the same run |
| V2 | Surface grids: seafloor age (batlow 0–150 Myr, as the maps), topo_bathy, divergence (diverging ramp centred on 0), plate ID; a layer dropdown; hover readout of value, plate ID and t | verification: decoded export against the NetCDF at three steps (difference maps, max error ≤ one uint8 step) |
| V3 | Mantle volume: depth slice, cutaway, isosurface (core `volume.ts`, `depthSlice.ts`, tomography cutaway), thresholds taken from the run's own recorded range | verification as V2, per depth; slabs under the trenches visible |
| V4 | Run picker (`runs.json`), run info panel (preset, modes, σc, tie seed…), legend, teaching defaults (which layer shows first) | you use it on the R2/R3 gate runs |

**Later, separate stages:** point sets (earthquakes by depth, seamounts, tracers, sutures),
linked globes to compare runs, interpolating plates between saved steps with the run's `.rot`
file, a shareable static build.

## Open questions (to decide at the stage that needs them)

- V2: are 0.5° grids good enough for teaching, or should 0.25° be the default?
- V4: which layer and run should students see first?
- Should the viewer later become a GeodeViewers member? That would add Pages deployment and a
  site card. It is out of scope for now.
