# Force model v1 — plan

Draft prepared with Claude Code. Decisions come from a design interview (Q1–Q18); the
question number is given in brackets.

## Goal

Plausible synthetic worlds first; speed matters more than strict realism (Q1). Forces
remain dimensionless, tunable weights, but take the physically correct *form*.
Every change is a mode flag. A `simple` preset reproduces today's model exactly: it
is the fastest option and useful for teaching.

## Decisions

| # | Topic | `physical` (default) | Also available | `simple` = today |
|---|---|---|---|---|
| Q2, Q19 | Plate mobility | `area_edge`: ω = k · D⁻¹ τ, D = basal area tensor + boundary (edge) tensor | `area` (basal only) | `none`: ω = τ |
| Q3 | Continental keels | w = `continent_drag_factor` (3) on continental (craton) cells, 1 on ocean | — | — |
| Q4 | Seafloor age | Existing age tracers, queried every step | — | age used for output only |
| Q5 | Slab pull magnitude | √min(age, 80)/√80 × min(1, L_slab/660 km) | age only | uniform |
| Q6 | Ridge push | Body force over ocean < 80 Myr | `boundary` | `boundary` |
| Q7 | Slab suction | Deep slab (> 660 km) → suction on *both* plates | `ratio` | overriding plate only, same weight as pull |
| Q8 | Pull on/off | Pull ∝ slab present (L_slab); convergence only seeds slab | `ramp`, `hard` | `hard` |
| Q9 | Polarity | Older ocean subducts; sticky per pair; continent rule first | `id` | lower plate ID |
| Q10 | Net rotation | Left free; logged per step | — | not logged |
| Q11 | Speed scale | `calibrate_once` at step 0 to `target_rms_speed`, then frozen; cap at 3× | `per_step`, `none`; `speed_scale` override | `none` |
| Q12 | Collisions | `breakoff`: slab detaches after `breakoff_delay_myr`, then coupling | `viscous` (replaces `resistance`), `coupling` | `resistance` |
| — | Transform resistance | Not in v1 (Forsyth & Uyeda: small) | — | — |
| Q13 | Physics grid | Separate `physics_grid_deg: 1.0`, built lazily | — | not built |

Rule (Q13): a field (physics grid, age raster, slab-length query) is built only if a
selected mode reads it. The `simple` preset must cost exactly what a step costs today.

## Config schema (Q14)

```yaml
force_model: physical            # preset: simple | physical
force_modes:                     # optional; any key here overrides the preset
  drag: area_edge                # none | area | area_edge
  ridge_push: body               # boundary | body
  slab_pull: slab_length         # hard | ramp | slab_length
  slab_suction: depth_split      # ratio | depth_split
  polarity: age                  # id | age
  collision: breakoff            # resistance | viscous | coupling | breakoff
  speed_normalisation: calibrate_once   # none | calibrate_once | per_step
params:
  continent_drag_factor: 3.0
  edge_drag_width_km: 3000       # boundary resistance (S1 sweep: 1000 / 3000 / 10000 km) = basal drag over a strip this wide
  target_rms_speed: 50           # km/Myr (= mm/yr)
  max_rms_speed_factor: 3.0      # cap = factor × target; each capped step is logged
  speed_scale: null              # a number fixes k and skips calibration (compare runs)
  slab_ref_depth_km: 660
  slab_age_saturation_myr: 80
  ridge_push_age_cutoff_myr: 80
  pull_ramp_km_per_myr: 10       # ramp mode only
  breakoff_delay_myr: 10
  physics_grid_deg: 1.0
```

When the run starts, the log prints the resolved modes and the k it used.

## Formulation notes

**Spherical geometry.** Everything here is on a sphere, and nothing uses a flat
approximation:

- **Tensors:** for each plate p, the drag tensor is I_p = Σ w (𝟙 − r rᵀ) dA over its
  cells; J_p is the same sum with w = 1. r is the unit position vector and dA uses
  cos(lat) weights on the physics grid.
- **Drag:** a rigid plate rotating at ω against a viscous mantle feels torque −c I ω,
  so the balance is ω = k I⁻¹ τ.
- **Boundary resistance (Q19):** with basal drag alone, resistance scales with area while
  driving forces scale with boundary length, so speed would go as 1/size. Forsyth & Uyeda
  (1975) found no clear relation between speed and plate area, and put the main
  resistance at plate edges. `area_edge` therefore adds an edge tensor
  ℓ Σ L (𝟙 − m mᵀ) over the plate's boundary segments, where ℓ = `edge_drag_width_km`.
  Setting `drag: area` drops it.
- **RMS speed:** found exactly from the tensors, with no velocity grid:
  RMS² = R² Σ_p ω_pᵀ J_p ω_p / A_total.
- **Net rotation:** ω_nr = (Σ_p J_p)⁻¹ Σ_p J_p ω_p. It is logged in `step_stats.csv`
  every step, as the rate and the pole, alongside the RMS speed.
- **Frame:** drag acts relative to the mantle, so the frame is fixed to the mantle. That
  is the same frame as the hotspots, and net rotation is a real output (Q10).

**Ridge push as a body force (Q6).** The force per unit area is f = +C ∇age for ocean
cells younger than the age cutoff, and 0 elsewhere. It points toward older seafloor, down
the subsidence slope. (An earlier version of this note had the sign wrong.) The plate's
torque is Σ r × f dA, and C = `forces.ridge_push` / cutoff, so one full flank gives the
boundary mode's force per unit ridge length.

- The gradient is taken within each plate only (never across a boundary). Per-plate
  spherical smoothing of the age before the gradient was tried in S4 and left off by
  default: it cost 0.2 s per step for little change (`ridge_push_smooth_deg`).
- The gradient is wrap-aware, with ∂/∂x = (1/(R cos φ)) ∂/∂λ.

*Correction to what I said in Q6:* the magnitude is ∝ |∇age|, not constant per unit
area. The force F(t) grows linearly with age, so the force per unit area is dF/dx. Its
integral from the ridge out to 80 Myr is therefore C·80 per unit length of ridge,
whatever the spreading rate. That matches Parsons & Richter (1980). A constant force
per area would give fast-spreading plates too much push.

**Age raster (Q4, Q13).** Assign each physics-grid cell the age of the nearest active
seafloor tracer, using one BallTree per step. Craton cells, and cells with no tracer
within 2 cells, get NaN.

**Slab length (Q5, Q7, Q8).**

- For each trench segment, take the deepest slab tracer within
  `boundary_search_radius_deg` of the segment, on its overriding side. Call its depth z.
  One BallTree on the slab tracers serves all segments.
- Pull = w_sp · f_age · min(1, z/660), where f_age is the age term from Q5.
- Suction = w_ss · min(1, max(0, z − 660)/660). It acts on both plates, toward the
  trench.
- The final definition will be fixed in S3, with a figure.

**Trench age (Q5).** This is the age of the nearest seafloor tracer on the downgoing
side, probed at an offset along the segment normal.

**Polarity (Q9).** Keep the existing continent rule in `resolve_downgoing_plate` as the
first test. For ocean–ocean pairs, at first convergence the plate whose ocean is older
along the contact goes down. The result is stored in `polarity_dict`, and its format
stays the same.

**Collisions (Q12).**

- `viscous`: push-back ∝ w_col · max(0, −div) / v_ref, which can stop plates but never
  reverse them.
- `breakoff`: once `collision_state` age exceeds `breakoff_delay_myr`, deactivate the
  slab tracers attached to the pair so the pull goes to zero, then apply coupling.

## Speed work (Q15)

1. Store tracers as numpy arrays instead of GeoDataFrames, in `advect_tracers`,
   `deactivate_subducted_tracers` and `seed_ridge_tracers`. They currently take about
   0.3–0.4 s per step, and the target is under 0.05 s (an estimate). The CSV output
   stays the same.
2. Share lookups within each step: one tree on the seafloor tracers and one on the slab
   tracers, reused by the age raster, trench ages and slab length.
3. After S5, spend a fixed amount of time profiling `resolve_overlaps`, `wrap_polygons`,
   `plate_id_raster` and `fragment`, at 0.1–0.3 s each.
4. No GPU work in v1: the work is many small geometry operations, and PyTorch in pygmt17
   is broken.

## Validation (Q16)

**Test case:** atlantis66 (the current `config.yaml`), 100 steps, dt = 1, fixed
`random_seed`.

`diagnose_run.py <output_dir>` makes PyGMT figures:

0. **Map sequence (the main check):** the plate configuration at each saved step, showing
   seafloor age, ridges, trenches with polarity, velocities and plate IDs. When runs are
   compared, each run gets a column at the same model times. Gate figures are presented
   through `show-me/`.
1. **Forsyth & Uyeda plot:** plate speed against trench fraction (the trend should be
   positive) and against continental fraction (the trend should be negative).
2. **Time series:** RMS speed, net rotation, maximum plate speed and step time.
3. **Map snapshots:** velocity arrows over the age grid, with trenches coloured by slab
   length.
4. **Regression test:** `force_model: simple` must reproduce the reference run exactly
   (rotations, plate IDs, tracers) over 10 steps. The reference was reset at 102d8b9,
   after the S2 tie-breaking rule. Before that, ties between equally near cells were
   broken by whatever the search method happened to do, and the model amplifies such
   differences into a different trajectory (see `show-me/2026-09-27-s2-faster-search/`).

Each gate shows the new mode side by side with `simple`.

## Stages (Q17)

Each stage ends at a review gate, and nothing moves to the next stage without approval.

| Stage | Contents | Gate |
|---|---|---|
| S0 | `diagnose_run.py`; preset and mode plumbing; mode log; save regression baseline; clean up warnings (`SettingWithCopyWarning`, the concat `FutureWarning` at run_SDP.py:603) | Figures for the current model; regression passes |
| S1 | Drag tensor with continent factor; speed calibration, cap and `speed_scale`; `net_rotation.csv` | F&U plot and time series compared with S0 |
| S2 | Numpy tracers; shared lookups; lazy physics grid and age raster | Step time; age grid unchanged; regression passes |
| S3 | Slab-length pull, depth-split suction, age polarity | F&U plot; slab-length maps |
| S4 | Ridge-push body force | Maps and time series |
| S5 | `breakoff` and `viscous` collisions | Collision sequence maps |
| S6 | Profiling pass | Step-time table |
| S7 | Update README and config documentation, with the sources below | Read-through |

## Branches (Q18)

- `review-fixes` holds the earlier review fixes and speed-ups, committed as 31832cb.
- `force-model-v1` branches from it, and this work goes there.
- Nothing is merged to `master` without approval.

## Found in S5

- **`viscous` is implicit** (user decision): the contact drag enters the same solve as basal
  and edge drag, as a coupling between the two plates on their relative convergence. It
  cannot overshoot, and sliding along the contact stays free. The strength is
  `collision_drag_factor: 10` × the edge drag. The plan's explicit push-back was dropped,
  because it would chatter like `resistance`.
- **`breakoff` on coupling slows the world down.** Velocity coupling ignores contact
  length, so a microcontinent on the largest ocean plate welded it to the supercontinent:
  RMS fell from 50 to 24 km/Myr over 200 Myr. `breakoff_contact: viscous` (drag ∝ contact
  length) holds RMS at about 39. Default pending the user's decision.
- **The supercontinent never breaks up** in 200 Myr of `physical`: nothing in v1 rifts a
  continent.
- **Cost grows late in physical runs** (2.75 s per step over the last 50 of 200, against
  1.56 s for simple). Seafloor tracers (151k) and slab tracers (127k) accumulate. S6 targets:
  thin the tracers, and reduce the slab-depth query to one point per slab column.

## Found in S4

- **Collision resistance chatters.** `collision: resistance` switches a push-back of
  weight 18 on at every converging continent–continent contact. The contact turns
  divergent, the push stops, and the contact closes again. This makes period-2
  oscillations of whole plates. It is the main source of the large step-to-step pole
  swings: in 100-Myr runs the share of plate area swinging more than 10° per step is
  10.2% (S3) and 22.5% (S3 + S4) with resistance, against 2.5% and 1.0% with `coupling`.
  S5's `viscous` mode (a push-back that can stop plates but never reverse them) targets
  this. Until then, compare force modes with `collision: coupling`.

## Decisions made during S3

- **Slabs at step 0:** each step-0 trench gets the slab column the model would have built,
  down to `initial_slab_depth_km: 1320`, so pull and suction are at full strength from the
  start. RMS still falls from about 50 to 44 km/Myr over the first 20 Myr as trenches move
  off those slabs.
- **Deep-slab suction weight:** `forces.slab_suction_deep: 3`, half of pull. Acting on both
  plates, that keeps the pair's total at today's 6.
- **Slab search:** the deepest slab tracer within 5° of the segment, on the overriding side
  (plus 1° across the trench line).
- **Polarity (`polarity: age`):** the continent rule is applied at the contact, not to the
  whole plate. A plate is continental at a contact if its own continent lies within
  `boundary_search_radius_deg` along at least half the contact length. It then overrides; if
  neither side is, the older seafloor subducts. The decision is kept, except that it flips when
  the downgoing side becomes continental and the overriding side is not; it never flips at a
  collision. The stored decision is final for every consumer (topography, earthquakes, tracers).
  With the whole-plate rule, age decided only 3 pairs in 100 Myr on atlantis66; at the contact
  it decides most ocean–ocean contacts.
- **Age grid:** each cell takes the nearest tracer *on its own plate*. Otherwise the
  overriding plate takes the downgoing plate's ages at trenches.

## Known limitations (accepted for v1)

- **Slab memory has no detachment.** A slab left behind by a trench that turned divergent
  pulls again at full strength if the boundary converges again while the slab tracers
  last (up to about 190 Myr). A new trench near an old slab also takes its pull, as
  does a nearby trench's slab within 5°. Breakoff after collision comes in S5.
- **Slab depth lags one step** (the slab tracers are updated after the forces).
- **Polarity is decided once per plate pair** (at first convergence, by age or by continent
  at the contact). Later the ages change, and they vary along a contact by about 18 Myr, so
  the older side subducts at only about half the ocean–ocean trench length after 50 Myr.
  Per-segment polarity is a v2 item (Q9c).

- **Ridges are not split into spreading segments and transforms.** Boundaries are continuous
  lines whose obliquity varies continuously: in the S0 run the length-weighted mean is about 40°,
  and about 40% of ridge and trench length is more oblique than 45°. The discrete
  ridge / trench / other classification therefore puts almost everything into ridge or trench.
  This keeps runs fast at low resolution. Splitting ridges into orthogonal spreading segments
  offset by transforms can be added later as a post-processing, appearance-only step on the
  saved boundaries.

## Known input issue

`config.yaml` pairs `atlantis3_rotations.rot`, which defines plates 101–901, with the
atlantis66 topologies and cratons, which use plates 111–777. This does no harm at
`reconstruction_time: 0`, because every rotation is then the identity. At any start
time above 0, those plates would silently get no rotation.

## Sources

- **Textbooks:**
  - Turcotte & Schubert, *Geodynamics* (3rd ed., 2014)
  - Kearey, Klepeis & Vine, *Global Tectonics* (3rd ed., 2009)
- **Drag and keels:** Solomon & Sleep (1974); Forsyth & Uyeda (1975); Chapple & Tullis
  (1977); Zhong (2001).
- **Slab pull and suction:** Carlson, Hilde & Uyeda (1983); Schellart (2004); Conrad &
  Lithgow-Bertelloni (2002, 2004).
- **Ridge push:** Lister (1975); Parsons & Richter (1980).
- **Coupled solutions:** Hager & O'Connell (1981); Ricard & Vigny (1989).
- **Net rotation:** Argus & Gordon (1991); Gripp & Gordon (2002).
