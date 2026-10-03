# PlateDriver viewer

A local globe viewer for PlateDriver runs. Plan and stages: [`../viewer_plan.md`](../viewer_plan.md).

```
cd viewer
npm install
conda run -n pygmt17 python prep/export_run.py <run_dir> [--name NAME] [--label LABEL]
npm run dev          # http://localhost:5173
npm run typecheck
```

Panel: Projection (Globe, Robinson, Plate Carrée), layer toggles (continents, belts, rift
strength with its opacity, boundaries, velocity arrows, plate IDs), Diagnostics (the plot above
the time bar and which series it shows) and Run (every exported run).

Diagnostics plot: one strip per series at every model step (number of plates, and those under
0.5% of the sphere; RMS speed; trench and ridge length; largest plate; net rotation; continent
hidden by overlap), with rift (red) and weld (blue) times as ticks along the top. Its time axis
is the slider's: hover to read every series, click to move the slider there.

URL parameters: `?run=<name>` (default: the last run exported), `?t=<Myr>`,
`?view=<lon>,<lat>` (globe camera at load), `?proj=robinson|plateCarree` (starting projection),
`?layer=strength` (rift strength on at load).

## Layout

- `src/core/` — copied (not imported) from Geode `viewer/src/core` at Geode commit `bbaa7ab`.
- `vendor/petrify/js/` — copied from petrify commit `45768d0` (v0.15.0+1), the one Geode `bbaa7ab` pins.
- `src/platedriver/` — this viewer's own code: `main.ts`, `timeBar.ts`, `plateOverlay.ts`
  (velocity arrows and plate-ID labels), `continentLayer.ts` (continent class grid, every
  projection), `scalarLayer.ts` (a scalar grid through a colour map, and its colour bar),
  `diagnosticsPanel.ts` (the diagnostics plot), `runData.ts`.
- `prep/export_run.py` — converts one run folder; runs in `pygmt17`.
- `public/data/runs.json`, `public/data/runs/<name>/` — exported runs (the run folders are gitignored).
- `scripts/shoot.mjs <out.png> [url]` — screenshot the running dev server.
- `scripts/check_scrub.mjs [url]` — scrub every saved step, fail on page errors.
- `scripts/check_projection_switch.mjs` — switch projection with the real dropdown, fail on page errors.
- `scripts/check_run_switch.mjs [url]` — switch through every run with the Run dropdown, fail on page errors.

## Exported run (so far)

```
run.json                  force model and modes, forces, params, saved steps (step, time, counts)
boundaries/manifest.json  petrify BoundarySeries frames, one per saved step
boundaries/NNNN.geojson   boundary polylines: boundary_type (ridge | subduction | other), plates,
                          and for trenches downgoing_plate, overriding_plate, polarity (Left | Right)
plates/NNNN.json          per plate: pole, rate, RMS speed, area and continental fraction, label lon/lat
velocities.json           petrify VelocityField: fixed points, east/north km/Myr per saved step
diagnostics.json          every step's whole-run series (step_stats.csv, plate_stats.csv) and the
                          rift and weld times, for the diagnostics plot
grids/continents/NNNN.bin uint8, 720 x 361 at 0.5°: 0 ocean, 1 continent, 2 craton core, 3 belt
grids/strength/NNNN.bin   uint8, 720 x 361: rift strength (self_driving_plates.rift_strength, the
                          rift test's strength per unit length of cut), byte x 0.02, 255 no data
```

Grids use the Geode Archive byte layout: longitude fastest from -180 (seam not repeated), then
latitude ascending from -90, both poles included. `run.json` gives `grid` (nlon, nlat) and `layers`.

- Continents are `craton_NNNN.nc > 0` (in PlateDriver the same cells as land and as no seafloor
  age). With `continent_strength: cores` each continental node takes the core/belt flag of its
  nearest `craton_cores_NNNN.csv` point, the lookup the model itself uses.

- PlateDriver's 2-point boundary segments are chained into polylines (petrify spaces trench
  triangles per feature) and densified along great circles to ≤ 1°.
- `conservative` is exported as `other`: it means "neither ridge nor subduction by the divergence
  thresholds", not a checked transform.
- Trench polarity, velocity arrows (ω × r) and label positions use the same rules as
  `diagnose_run.py`'s maps, so the two can be compared directly.
- Velocities are one `velocities.json` for the run (VelocityField's format), not part of
  `plates/NNNN.json` as the plan sketched.
