# Evidence behind README statements

Text records of the tests behind statements in `README.md`. Each record is self-contained: it says what
was run, how it was measured, the numbers, what they show and what they do not, and whether the code has
changed since. Figures and run folders (`show-me/`, `runs/`) are local working folders and are not in the
repository; everything needed to judge or repeat a test is written here.

Conventions used in every record:
- Time is model time elapsed since the start of the run (Myr), not Ma.
- Areas are shares of the sphere, cos(latitude)-weighted on the lon/lat grid.
- "Tie seed N" is `tie_break_seed: N`. "Step" is one `delta_time`.
- "atlantis66" is the bundled starting map (`examples/atlantis/`), resolved at `reconstruction_time: 0`.
- A "sliver" is a plate smaller than 0.5% of the sphere.
- Commits are on branch `mantle-extract`. To repeat a test exactly, check out the commit given; numbers
  from older records are not expected to match today's defaults (see each record's Status).

| ID | Topic | Date | README section citing it |
|---|---|---|---|
| E1 | Boundaries snap to meridians and parallels (`plate_rebuild: lonlat`); ridge-push direction artefact | 2026-09-28 | Force modes: `plate_rebuild` |
| E2 | A plate reaching from pole to pole is turned inside out by the dateline wrap | 2026-10-01 | Force modes: `rebuild_tiles` |
| E3 | Continental rift strength σc sweep (`rifting: trial`, straight cuts) | 2026-09-28 | Rifting: "How σc behaves" |
| E4 | Least-cost rift path against straight cuts | 2026-09-29 | Force modes: `rift_path` |
| E5 | Sensitivity to the time step and the grids | 2026-11-17 | Known issues: "Results depend on the time step and the plate grid" |

---

## E1. Boundaries snap to meridians and parallels with `plate_rebuild: lonlat`

**Claim (README).** With `plate_rebuild: lonlat`, about 7× more boundary length lies within 1° of N–S or
E–W than for random orientations, because rebuilding the plates from the lon/lat raster every step snaps
near-aligned boundaries onto meridians and parallels.

**Code.** Commit `cbe91a8` (2026-09-28). At that commit every preset used the `lonlat` rebuild
(`plate_rebuild: midpoint` did not exist yet). Rebuild code: `update_plate_polygons` in
`self_driving_plates.py`.

**Mechanism.** Every step `compute_timestep` rebuilds the plates from the 0.25° plate-ID raster:
1. the moved plate polygons are rasterised;
2. `polygonize` returns staircase outlines along meridians and parallels;
3. `topojson` `toposimplify(topo_simplify_tolerance)`, tolerance 0.4 in degrees of lon/lat
   (Douglas–Peucker), simplifies them.

A boundary within a small angle of a grid line has staircase steps smaller than the tolerance on one axis,
so the simplification replaces it with an exact meridian or parallel. The rebuild repeats every step, so a
snapped boundary stays snapped unless the plates move it far enough.

**Measure: alignment excess.**
- Each boundary piece gets the azimuth at its midpoint on the sphere, atan2(Δlon·cos(lat), Δlat), folded to
  0–90° (0 = N–S, 90 = E–W), weighted by great-circle length.
- Excess = (length share within 1° of 0° or 90°) / (2/90). 1 means no preference; 7.5 means about 17% of
  boundary length within 1° of a grid direction against 2.2% expected.
- Excluded: pieces poleward of 80° latitude, and pieces within 0.5° of the dateline (the dateline seams are
  exactly N–S by construction).
- Check of the measure: random great-circle pieces give 1.02.

**Runs.** atlantis66, 1-Myr steps.

| Runs | Mean excess over steps ≥ 50 | Range over those steps |
|---|---|---|
| `physical`, 500 Myr, no tracer thinning | 7.3 | 4.3–10.1 |
| `physical`, 500 Myr, `tracer_thinning: cell` | 7.7 | 5.0–11.5 |
| `physical`, 200 Myr, 3 runs, no thinning | 7.3, 7.8, 9.4 | 5.7–12.6 |
| `simple`, 200 Myr, 3 runs | 4.4, 4.5, 4.8 | 2.5–7.6 |

**Results.**
- The input topology, resolved at 0 Ma before any rebuild, has an excess of 1.76. After one rebuild (step 0)
  it is 2.41. Runs reach their working level within about 30 Myr.
- **No-force test (decisive).** The step-0 plates were given 25 passes of a random rigid rotation of the
  whole sphere, each followed by the model's own rebuild, with no forces. Rigid rotation alone cannot create
  a preferred orientation. The excess rose pass by pass from 2.4 to 12.8 at 0.3° per pass (about 33 km),
  and from 1.5 to 5.5 at 1° per pass (about 110 km). Measured on the polygons just after each exact
  rotation (before the rebuild), the excess equals the previous pass's, so the measure itself adds nothing.
  So the rebuild creates the alignment, and slow motion lets it accumulate.
- `physical` is about 1.6× more aligned than `simple`, consistent with the no-force test: `physical` plates
  move about 3× slower, so snapped boundaries stay snapped longer.
- Orientation histograms (500-Myr `physical`, steps 50–500): ridges and trenches both peak in the first 1°
  bin at 0° and at 90° at 6–9× the isotropic density, with a deficit at intermediate angles.
- Snapped boundaries form straight runs of up to several thousand km, for example about 45° of longitude
  along 35°N at 400 Myr in the run with thinning.
- Tracer thinning does not change it: 7.7 against 7.3, inside the 7.3–9.4 of three runs without thinning.

**Second artefact: ridge-push direction.** The body ridge-push force on the 1° physics grid, recomputed
offline from the seafloor tracers at 100 and 200 Myr, points within 1° of N–S or E–W 11.4 and 11.1× more
often than an isotropic force would (10.1 and 10.4 with thinned tracers). Cause: the age grid takes each
cell's nearest tracer, so it is piecewise constant, and its gradient is taken by finite differences along
longitude and latitude; where only one neighbour differs, the gradient is exactly E–W or N–S. Its effect on
plate motion was not measured; integrated over a plate, much of it may cancel.

**Limits.** Neither artefact's effect on plate motions was measured. The snapping does change the
boundaries that forces act on and how maps look.

**Status.**
- Snapping: fixed for `physical` by `plate_rebuild: midpoint` (commit `fb0d078`, 2026-09-28), which moves
  every staircase vertex to the midpoint of its edge before simplifying. Three 500-Myr runs gave an excess of
  1.2 against 7.3 with `lonlat`; `midpoint` boundaries are about 15% shorter in total and mean plate counts
  slightly lower (8.7–9.5 against 10.1–10.9). One `midpoint` run rose to about 4 in its last 40 Myr; this was
  not followed up. `simple` still uses `lonlat`.
- Ridge-push direction: open. `ridge_push_smooth_deg` (default 0; 1.5 costs about 0.2 s per step) smooths
  age per plate before the gradient. Stage F2 of `docs/public_release_plan.md` (ridge push from a gradient
  on the sphere) is the planned fix.

---

## E2. A plate reaching from pole to pole is turned inside out by the dateline wrap

**Claim (README).** With `rebuild_tiles: off`, a plate that touches both edges of the lon/lat map and a
pole is replaced by its complement when rebuilt, and its cells go to its neighbours, with no error.
`rebuild_tiles: on` prevents it.

**Code.** Found 2026-10-01 in run `c1_weld` (R5 set-up with `collision_weld: on`, code of 2026-09-30).
`rebuild_tiles` added in `f0dfce9` (2026-10-03); `physical` default since `c15aeb6` (2026-11-16). Code:
`_tile_index` and `update_plate_polygons(..., tiles=True)` in `self_driving_plates.py`.

**What happened in the run** (log and `plate_stats.csv`, every step):
- Plate 777 had absorbed most plates by welds and covered 62% of the sphere at 238 Myr. At 243 Myr plate 444
  (23%) was welded into it, so 777 covered 81–87% at 244–246 Myr.
- From 244 Myr plates split into "disconnected pieces" every step (6, 8, 11, 3 and 16 pieces at 244, 245,
  248, 249 and 250 Myr), and plates enclosed by 777 were absorbed. At 246 Myr two plates were left (777 87%,
  555 13%); at 247 Myr 555 was 48% and 777 52%, interleaved in long tongues with thin strips of new plates
  between them.

**Cause, tested on the saved 240-Myr plate grid with no motion** (one rebuild: grid → polygons → dateline
wrap → grid):
- As saved (777 = 62%): 0.2% of cells change plate.
- With 444 welded into 777 by hand (85%): 33% of cells change; 555 goes from 2.7% to 17%, 112 from 3.4% to
  7%, 788 and 790 vanish.
- With only 788 welded into 777 (69%, 777 then covers both poles): 51% of cells change. 777's outline then
  holds a segment between antipodal points, which pygplates rejects as a polygon.

**Which stage fails.** Area of plate 777 after each rebuild stage, measured by burning that stage's
polygons back onto the grid (cos-latitude weights):

| 240-Myr grid | 777 in the grid | after lon/lat polygonising | after the dateline wrap |
|---|---|---|---|
| as saved | 62.1% | 62.1% (all on 777's cells) | 62.1% (all on 777's cells) |
| 444 welded into 777 | 85.4% | 85.4% (all on 777's cells) | 14.6%, none on 777's cells |
| 788 welded into 777 | 68.9% | 68.9% (all on 777's cells) | 31.1%, none on 777's cells |

The lon/lat polygons are correct. pygplates' dateline wrap returns exactly the complement (100 − 85.4,
100 − 68.9). The outline of a plate touching both map edges and a pole runs up the dateline, down it again
(the same meridian on the sphere) and along the pole row (one point on the sphere), so the wrap decides the
inside wrongly. All other plates came through unchanged. The code already caught one form of this and
stopped with an error (`63a237b`); these forms passed silently and gap filling handed the unowned cells to
neighbours, compounding over the following steps.

**Fix tested on the same grids.** Cut every plate into tiles of 90° longitude × 45° latitude before
polygonising (a temporary ID per plate and tile), then give the pieces back their plate ID. No piece can
hold a pole or touch both map edges.

| 240-Myr grid | cells changing plate in one no-motion rebuild, untiled | tiled |
|---|---|---|
| as saved | 0.22% | 0.14% |
| 444 welded into 777 | 32.9% | 0.08% |
| 788 welded into 777 | 51.0% | 0.11% |

**Status.** Fixed by `rebuild_tiles: on`, the `physical` default. The tiles' shared edges are not plate
boundaries and are left out of the force sum. `simple` keeps `off`.

---

## E3. Continental rift strength σc sweep (`rifting: trial`, straight cuts)

**Claim (README).** On atlantis66 over 500 Myr, σc 1.0 gives practically no rifting, σc 0.5 gives a few
dozen rifts with the large continent 222 rifting and continents dispersing, and σc 0.25 runs away.

**Code.** `86bd3b2` (2026-09-28) for the three σc runs; the no-rifting reference is `fb0d078` (identical
config apart from rifting). The `physical` preset of that date: `plate_rebuild: midpoint`,
`overlap: plate_id`, no `rebuild_tiles`, no welding.

**The rift test as run.**
- Each step every plate is tested against 36 great-circle cuts: 12 azimuths × 3 offsets (−15°, 0°, +15°)
  about its centroid.
- A cut passes if the tension across it (from the plate's solved torque balance) exceeds the mean strength
  along it, the two free halves would open along at least 50% of the cut, and both halves are at least
  `min_child_area_deg2`.
- Strength: continent σc (`rift_strength_continent`); ocean 3σc·max(0.2, √(min(age, 80)/80)).
- A passing cut must persist for 3 Myr (`rift_persist_myr`); both halves then get a 10-Myr cooldown
  (`rift_cooldown_myr`).

**Runs.** atlantis66, `physical`, 500 steps of 1 Myr, `random_seed: 42`, tie seed 0 (σc 0.5 also tie seeds
1 and 2), 0.25° raster, saved every 10 steps.

**Results.**

| Rifts per 100 Myr | 0–100 | 100–200 | 200–300 | 300–400 | 400–500 | total |
|---|---|---|---|---|---|---|
| σc 1.0 | 0 | 0 | 1 | 0 | 0 | 1 |
| σc 0.5 | 2 | 7 | 7 | 3 | 5 | 24 |
| σc 0.25 | 30 | 34 | 68 | 133 | 265 | 530 |

| σc | new plates | median lifetime | share living ≥ 50 Myr | alive at 500 Myr |
|---|---|---|---|---|
| 1.0 | 1 | 11 Myr | 0% | 0 |
| 0.5 | 24 | 48 Myr | 46% | 8 |
| 0.25 | 530 | 22 Myr | 30% | 154 |

| | no rifting | σc 1.0 | σc 0.5 | σc 0.25 |
|---|---|---|---|---|
| RMS speed (km/Myr) | 37.9 | 38.0 | 42.0 | 38.6 |
| pole swings > 10° (% of area per step) | 2.0 | 2.0 | 3.0 | 10.7 |
| mean plates | 8.7 | 8.8 | 13.0 | 72.5 |
| r(speed, trench fraction) | +0.83 | +0.86 | +0.84 | +0.60 |
| r(speed, continental fraction) | −0.87 | −0.87 | −0.67 | −0.68 |
| seconds per step | 1.28 | 1.55 | 1.74 | 3.19 |
| rift test, seconds per step | none | 0.21 | 0.22 | 0.30 |

- σc 1.0: one rift (an oceanic sliver of plate 112 at 299 Myr that lasted 11 Myr); otherwise the same world
  as no rifting, so the trial mode leaves an unsplit world unchanged.
- σc 0.5: plate 222 rifts at 275 Myr (child 4.6% of the sphere) and 406 Myr (3.6%); 19 plates at 500 Myr
  against 8 without rifting. Plates of at least 1% of the sphere stay at 7–10, within 1–2 of no rifting; the
  extra plates are small. The largest continental block falls from 52% to 28% of continental area between
  380 and 430 Myr (7 blocks at 500 Myr; no-rifting runs end at 44–61% with 4 blocks).
- σc 0.25: cuts start within 10 Myr, children split again once their cooldown ends, 182 plates at 500 Myr.
  Rejected.
- **Most rifts at σc 0.5 are oceanic slivers.** 15 of 24 children have no continent; children are 0.1–1% of
  the sphere and most last 25–50 Myr. Cause, from the strength law: ocean younger than about 9 Myr
  (3√(age/80) < 1) is weaker than continent, reaching the floor of 0.6σc below about 3 Myr, so a cut along
  a ridge flank is cheaper than a cut through a continent. The continental rifts are the long-lived ones.

**σc 0.5 over three tie seeds.**

| tie seed | rifts | 222 rifts at (Myr) | new plates living ≥ 50 Myr | oceanic children | largest continent at 500 Myr (% of continental area) | blocks | mean plates |
|---|---|---|---|---|---|---|---|
| 0 | 24 | 275, 406 | 11 | 16 | 28 | 7 | 13.0 |
| 1 | 28 | 273, 379 | 13 | 21 | 36 | 6 | 13.4 |
| 2 | 21 | none | 9 | 20 | 48 | 3 | 11.7 |

Early events repeat in all three runs (333 rifts at 10 Myr; 115 at 6–81 Myr; 666 three times at 102–150
Myr; 555 at 179–236 Myr). Plate 222's rifting is not robust: absent in tie seed 2, which has no rifts
between 237 and 426 Myr. In every run 67–95% of new plates have no continent at birth.

**Conclusion.** σc 0.5 was chosen (`rift_strength_continent: 0.5`, still the default).

**Limits.** One run per σc for the sweep; runs diverge after about 100 Myr, so a single late rift is not a
reproducible event.

**Status.** The numbers predate `overlap: polarity`, `rebuild_tiles`, merging and the later rift options, so
they describe the behaviour of σc, not today's default world. Large ocean plates still never rift and small
plates still rift more easily (README known issues).

---

## E4. Least-cost rift path against straight cuts

**Claim (README).** On atlantis66, `rift_path: least_cost` gives 1.7× as many rifts as straight cuts, mostly
along young ocean, and continents break up less.

**Code.** Straight cuts: the E3 runs (`86bd3b2`). `least_cost`: `f1be8d8` (2026-09-29).

**What `least_cost` does.** Each candidate great circle whose straight cut is under tension gets the
cheapest path across the plate between the straight cut's two ends, within `rift_path_corridor_deg` (5°) of
the straight cut, on the 1° physics grid; a step's cost is its great-circle length × the mean strength of
its two cells. The trial-split test is then rerun on the path (halves by flood fill, tension from the path's
own normals, strength along the path). A split planned for the next step keeps both halves' cells and moves
them with the plate.

**Runs.** atlantis66, `physical`, 500 steps of 1 Myr, tie seeds 0, 1 and 2, σc 0.5, no welding.

**Results.**

| | straight | `least_cost` |
|---|---|---|
| rifts per run (tie 0, 1, 2) | 24, 28, 21 | 31, 49, 45 |
| continental rifts per run | 7, 7, 1 | 3, 12, 13 |
| new plates (3 runs) | 73 | 125 |
| median lifetime | 43 Myr | 46 Myr |
| share living ≥ 50 Myr | 45% | 45% |
| child ≥ 10% continent | 21% | 22% |
| median area at birth | 0.15% of sphere | 0.16% of sphere |
| plate 222 rifts at (tie 0, 1, 2) | 275 and 406; 273 and 379; never | 479; 314; never |
| largest continent at 500 Myr (% of continental area) | 28, 36, 48 | 64, 39, 66 |
| mean plates | 11.7–13.5 | 12.7–16.6 |
| largest plate ever (% of sphere) | 31–33 | 33–34 |
| RMS speed (km/Myr) | 40–42 | 39–41 |
| r(speed, continental fraction), plates living ≥ 50 Myr | −0.56 to −0.67 | −0.37 to −0.53 |
| rift test, seconds per step | 0.21 | 0.67–0.71 |
| seconds per step | 1.55–1.74 | 2.22–2.41 |

- `least_cost` makes 1.7× as many rifts, of the same kind as straight cuts: small, mostly oceanic, with the
  same share surviving.
- Most new `least_cost` boundaries are long ridge-parallel lines in the ocean: the cheapest path runs through
  the youngest seafloor, which is weaker than continent (E3). Continental rifts are mostly short cuts at
  continental margins.
- Continents stay more together with `least_cost` (largest continent 39–66% against 28–48%), because the
  cheapest cuts go through ocean instead.
- The path search adds about 0.5 s per step, about 35% more run time.

**Limits.** Welding was off, so following sutures (the reason for the path) was not tested. The path uses no
craton or suture weakening beyond what is in the strength field.

**Status.** `physical` keeps `rift_path: straight`; `configs/current_model.yaml` uses `least_cost` with
`continent_strength: cores_random`. Since `a815951` (2026-11-17) the path finds the side of each segment 0.7
grid cells off the path instead of a fixed 0.7°; no change on the 1° grid.

---

## E5. Sensitivity to the time step and the grids

**Claim (README).** Results depend on the time step and the plate grid but not on the force grid at its 1°
default. A finer time step leaves about twice as much of the starting ocean unsubducted at 200 Myr; the
effect is gone by 500 Myr. The trenches consume the right area at every step; what changes is which ocean
they consume.

**Code.** 2026-11-17, commits up to `1c51946` (the 2° force-grid row with the fixed offset uses `a815951`).
The reference is eleven seeds of the current model (`runs/overlap_polarity` and `runs/seeds_current`); the
code of 17 November reproduces the 16 November runs bit for bit.

**Runs.** `configs/current_model.yaml`, tie seeds 0–2 unless stated, one numerical setting changed at a time.
Parameters counted in steps were scaled with the time step so they cover the same time.

| Variant | Overrides of `configs/current_model.yaml` | Length |
|---|---|---|
| time step 0.25 Myr | `delta_time: 0.25`, `merge_min_steps: 20`, `tracer_thin_interval: 40`, `csv_save_interval: 40` | 200 Myr |
| time step 0.5 Myr | `delta_time: 0.5`, `merge_min_steps: 10`, `tracer_thin_interval: 20`, `csv_save_interval: 20` | 300 Myr (8 seeds to 200 Myr; 4 seeds to 500 Myr) |
| time step 2 Myr | `delta_time: 2.0`, `merge_min_steps: 3`, `tracer_thin_interval: 5`, `csv_save_interval: 5` | 300 Myr |
| plate grid 0.125° | `raster_spacing: 0.125` | 200 Myr (6 seeds) |
| plate grid 0.5° | `raster_spacing: 0.5` | 300 Myr |
| force grid 0.5° | `physics_grid_deg: 0.5` | 300 Myr |
| force grid 2° | `physics_grid_deg: 2.0` | 300 Myr (before and after `a815951`) |

Repeat with `gate.py`, for example:
`conda run -n pygmt17 python gate.py <variant config> --name sens_dt_half --seeds 0 1 2 --until 300`.
`gate.py` maps times to steps correctly for any `delta_time` only since `741e72b`.

**Measures** (from `gate.py`'s table): ocean older than 200 Myr (share of the ocean, from
`seafloor_age_NNNN.nc`; at 200 Myr this is exactly the starting ocean that has not subducted); median ocean
age; mean plate count; slivers; rifts.

**At 200 Myr** (mean over seeds 0–2, range in brackets; defaults over eleven seeds):

| | ocean > 200 Myr, % | ocean age median, Myr | plates (100–200 Myr) | slivers | rifts (100–200 Myr) |
|---|---|---|---|---|---|
| defaults (1 Myr, 0.25°, 1°) | 4.9 [0.5–9.6] | 68 [53–83] | 16.1 [11.7–18.1] | 2.5 [1.0–3.9] | 12.2 [8–16] |
| time step 0.25 Myr | 10.6 [9.8–11.4] | 61 [56–64] | 15.8 [14.4–17.1] | 3.1 [2.1–4.3] | 7.7 [7–9] |
| time step 0.5 Myr (8 seeds) | 11.0 [9.6–12.3] | 85 [74–91] | 17.6 [15.3–19.6] | 2.4 [1.4–3.3] | 13.8 [12–16] |
| time step 2 Myr | 3.0 [1.5–5.4] | 59 [57–60] | 11.3 [10.3–11.9] | 0.7 [0.1–1.1] | 3.7 [1–7] |
| plate grid 0.125° (6 seeds) | 11.5 [7.5–13.6] | 91 [86–97] | 16.4 [14.9–17.3] | 2.6 [2.2–3.4] | 13.5 [12–15] |
| plate grid 0.5° | 9.9 [5.5–12.3] | 68 [52–81] | 14.8 [13.5–16.8] | 2.3 [1.7–3.6] | 13.0 [12–14] |
| force grid 0.5° | 4.3 [0.8–9.2] | 69 [60–84] | 17.0 [15.6–19.4] | 3.0 [2.5–3.8] | 9.7 [8–12] |
| force grid 2° | 7.9 [3.7–12.0] | 69 [58–80] | 20.6 [19.3–21.8] | 4.6 [3.6–5.1] | 21.3 [17–27] |
| force grid 2°, offset fixed (`a815951`) | 9.0 [6.3–10.5] | 70 [51–81] | 18.0 [13.2–21.1] | 3.6 [1.0–5.9] | 13.7 [9–17] |

**Results.**
- **Force grid converged at 1°.** 0.5° changes no measure beyond the seed spread, at 2.4× the cost per step.
  2° adds plates, slivers and rifts; part of that came from the least-cost side offset being fixed at 0.7°
  (set for the 1° grid). With the offset at 0.7 grid cells (`a815951`), 2° rifts are within the spread at
  200 Myr but not at 300 Myr (37 against 29 over 0–300 Myr).
- **Time step not converged.** Starting ocean left at 200 Myr: 3.0% at 2 Myr, 4.9% at 1 Myr, 11.0% at
  0.5 Myr, 10.6% at 0.25 Myr. Monotonic, and the last two agree, which suggests convergence below 0.5 Myr at
  about twice the default's value. The median ocean at 0.5 Myr is about 20 Myr older.
- **Plate grid not converged, and not monotonic.** 0.125° (11.5%) and 0.5° (9.9%) both leave more starting
  ocean than the default 0.25° (4.9%), so the default sits at the low end.
- **Not an overlap-rule error.** At every step size and grid, the area settled as overlap per Myr equals the
  trenches' convergence × length to 98–99% (from `step_stats.csv` and `trench_forces_NNNN.csv`). What
  changes is which ocean subducts: from 100 to 200 Myr the ocean arriving at trenches averages 72–73 Myr old
  at 1 Myr, 63–65 Myr at 0.5 Myr and 66 Myr on the 0.125° grid. At 0.25 Myr it is 74 Myr old and 11% still
  survives, so no single mechanism explains every row.
- **The effect is transient.** Four seeds at 0.5 Myr to 500 Myr: 5.8% of the starting ocean left at
  250 Myr against 1.9% over eleven default seeds; by 500 Myr 0.2% at both. Ocean older than 200 Myr at
  500 Myr: 10.1% (5.4–13.5) against 15.3% (4.0–29.1), within the seed spread. The largest plate peaks
  somewhat higher at 0.5 Myr (32.6% against 27.5% over 250–500 Myr). So a finer step changes how fast the
  starting ocean is consumed, not the old ocean that builds up later; force model v3's case (old ocean
  accumulates without subduction initiation) does not rest on the time step.
- **Plate counts.** At 300 Myr a 0.5-Myr step or a 0.5° plate grid gives more plates (24–27 against 13–23)
  and about twice as many slivers. At 0.5 Myr the extra plates come from fewer merges (6–13 in 300 Myr
  against 14–15), not more rifts: a merge needs the pair's boundary to stay nearly conservative for a run of
  steps, and with more, shorter steps one noisy step breaks the run more often. Rifts fall at both 0.25 and
  2 Myr.

**Open.** Which per-step process steers trenches towards younger ocean at a finer step is unknown. Candidates
are rules that decide once per step: polarity for new plate pairs, the rift test's persistence, merging, and
rebuilding plates from the grid. Polarity churn and trench consumption were checked and ruled out.

**Guidance.** The defaults (1 Myr, 0.25°, 1°) are the values the layout rules were tuned with. Treat early
(< about 250 Myr) old-ocean shares and late plate counts from 1-Myr runs as unconverged; check such a claim at
`delta_time: 0.5` before relying on it.
