# Force model v2 — plan: rifting and plate splitting

Draft prepared with Claude Code. Decisions come from a design interview (Q1–Q15); the
question number is given in brackets.

## Goal

v1 never breaks a continent. In every 200-Myr `physical` run on atlantis66 (three tie
seeds), the supercontinent plate 222 stays whole and RMS speed drifts from 50 to 33–35 km/Myr.
The existing plate splitting only ever makes short-lived slivers: 5–8 new plates per 200 Myr
in physical, a median size of 0.04–0.10% of the sphere and a median life of 6–11 Myr. None
reaches 1% of the sphere and none is continental. Two things cause this:

- The cut is purely geometric (a thin neck), so the two halves need not move apart.
- Where plates overlap, the lower plate ID wins. A new plate always has the highest ID, so it
  loses area at every convergent contact.

v2 adds one splitting mechanism for **any** plate, continental or oceanic: a plate splits
where the forces would pull its two halves apart harder than the lithosphere can resist.
Every change is a mode flag, and `simple` keeps today's behaviour (Q14).

## Decisions

| # | Topic | Decision | Options kept |
|---|---|---|---|
| Q1 | Trigger | **Trial split:** solve each candidate cut's two halves against the plate as a whole (see Formulation) | `stress` prototype, a side experiment (a later stage) |
| Q2 | Candidate cuts | All three sources, each switchable: `sweep` (great circles), `weak_zones` (sutures etc., Q9), `necks` (rewritten on the sphere, without shapely) | any subset |
| Q3 | Split line | Least-cost path on the 1° physics grid within ±5° of the chosen cut. Cost is the strength field plus a little noise. The trial-split test is re-checked on the final line | — |
| Q4 | Resistance | The tension per km of cut is compared with the mean strength along the cut. Continental and oceanic strength differ | — |
| Q5 | Weakening | `rift_weakening: none \| cmb_heat` | default chosen from test runs |
| Q6 | Opening | Immediate: the cut is classified from motion, so a diverging cut is a ridge at once. For the first `rift_phase_myr` (20) it is labelled `rift` in the outputs (topography, earthquakes) | — |
| Q7 | Overlaps | `overlap: overriding`. At a trench the overriding plate keeps the overlap (stored polarity); at a collision the continent keeps it; otherwise the plate that held the cell last step. Done on the plate-ID raster, without shapely | `plate_id` (simple, lower ID wins) |
| Q8 | Screening | All plates above `min_child_area_deg2` × 2 are tested every step. A split needs the test passed on a similar cut for `rift_persist_myr` (3) consecutive steps. Both halves then get a cooldown | — |
| Q9 | Weak zones | **Suture tracers:** dropped on both plates along a continent–continent contact that has collided longer than `breakoff_delay_myr`. They are advected and never deleted, and written to `sutures_NNNN.geojson`. **Craton-ID lines:** where cratons of different original plates touch inside one plate. **Optional input file** of sutures for the starting world | — |
| Q10 | Welding | `collision_weld: off \| on`. After `weld_after_myr` (30) of continuous continental collision the smaller plate is absorbed; the suture tracers mark the join | default chosen from test runs |
| Q11 | Unstable shapes | Neck and ribbon cuts are only candidates and must pass the same force test. A short cut has a high pull per km, so necks split easily, but only where the halves separate | — |
| Q12 | Validation | atlantis66, 500 Myr; one run per development step, three tie seeds at each gate. The map sequence stays the main figure | — |
| Q13 | Strength | See Strength field | — |
| Q14 | `simple` | Rewrite the neck search and `resolve_overlaps(plate_id)` without shapely. Re-baseline once, after a statistical check against the existing three-seed `simple` ensemble | — |
| Q15 | Pattern defaults | See Strength field | — |

## Strength field (Q4, Q5, Q13, Q15)

This is the full field agreed in the interview. R2 starts with constant σc on continents and
the age-dependent ocean only; the other factors come later, one at a time (see Stages).

Strength σ is defined per cell on the 1° physics grid, in units of **full slab pull per unit
length** (the slab-pull weight, 6, per unit trench length). All distances are great-circle
distances on the sphere.

- **Continent:**
  σ = σc × max(interior, core) × suture × belts × weakening.
  - `interior`: rises linearly from 1 at the continent edge to 2 at 1000 km inland
    (`interior_strength_factor: 2`, `interior_strength_km: 1000`).
  - `core`: `craton_core_factor` (3) inside craton-core polygons from an optional input file
    (atlantis66 has none), carried with the plate. The larger of interior and core applies.
  - `suture`: 0.5 within 2° of a suture tracer or a craton-ID line
    (`suture_strength_factor: 0.5`, `suture_width_deg: 2`).
  - `belts`: 1 + 0.4·n, clipped to 0.3–1.7, where n is a unit-variance band-limited
    random pattern with features about 1500 km across (spherical harmonics, as for the
    terrain noise), fixed to the continent (`belt_amplitude: 0.4`, `belt_scale_km: 1500`).
    **This makes plate motions depend on `random_seed`,** which they don't today.
  - `weakening` (`cmb_heat` mode): 1 − 0.5 · heat / `cmb_max_heat`, i.e. strength is at most
    halved over the hottest mantle (`cmb_weakening: 0.5`).
- **Ocean:** σ = 3 σc × max(0.2, √(min(age, 80) / 80)) × weakening.
  - `oceanic_strength_ratio: 3`, `oceanic_strength_floor: 0.2`.
  - The √age saturation at 80 Myr matches slab pull (plate thickness grows as √age).
  - Old ocean is therefore 3× continent, 20-Myr ocean 1.5×, and ocean younger than about
    3 Myr 0.6×. Cuts along young ridge flanks are easier than through continents, which allows
    occasional ridge jumps.
- **σc:** a sweep of 0.25, 0.5 and 1.0 in 500-Myr runs (R2). The default is chosen from those
  runs.

Every factor has its own switch and can be set to 1 (off).

## Formulation

**Trial split (Q1, Q4).** Take a plate P with drag tensor D, torque τ and solved rotation ω
(from the global solve, including contact tensors). A candidate cut divides P's physics-grid
cells and boundary segments into halves 1 and 2. Then:

- D = D1 + D2 and τ = τ1 + τ2. Basal cells, body forces and boundary segments are assigned by
  side; the cut itself carries no edge force.
- The torque that half 2 must pass to half 1 to keep it moving with ω is τh = D1 ω − τ1.
- The cut has unit normal n̂(s), pointing from half 1 to half 2. A uniform normal traction t
  along the cut gives a torque on half 1 of t · M, where M = ∫ r × n̂ ds. The tension is then
  t = τh · M / |M|². Positive t means the halves are pulling apart.
- The cut **passes** if t > σ̄, the length-weighted mean strength along the cut.
- As a second check, the free halves, ωi = Di⁻¹ τi (with the same k), must open along at least
  half the cut length.

**Cost.** For the `sweep`, each plate's centroid c has 12 azimuths × 3 offsets (cut planes
through the origin with normal n cos δ + c sin δ, δ ∈ {−15°, 0, 15°} of arc). Cuts that
leave a half smaller than `min_child_area_deg2` are dropped. The per-candidate sums over cells
are one matrix product (candidates × cells) @ (cells × 9), so a whole step's test should cost
0.01–0.05 s. This is to be measured in R2.

**Split (Q3).**

- The least-cost path runs between the chosen cut's two end points on the plate edge. It is a
  Dijkstra search on the 1° grid (scipy `csgraph`), with great-circle edge lengths × mean cell
  strength, inside a ±5° corridor.
- Cells are assigned to a side by flood fill from the path.
- The test is re-run on the final line and the split is abandoned if it fails.
- The child gets a new plate ID and inherits the parent's rotation history, polarity decisions
  (parent, X) → (child, X), collision states along its part of the boundary, and its cratons.
- Tracers, seamounts and slab tracers take the plate ID of their cell.
- The cut's segments are labelled `rift` for `rift_phase_myr`.

**Craton-cut rule.** Today a cut that would separate cratons of the same plate is skipped. With
`rifting: trial` this rule is dropped; that is the point. A split craton keeps its original ID
on both sides, which later makes a craton-ID line if the halves are welded back together.

## Config schema

```yaml
force_modes:
  rifting: trial                 # none (simple: today's geometric neck splits) | trial
  rift_candidates: [sweep, necks, weak_zones]
  rift_weakening: cmb_heat       # none | cmb_heat
  overlap: overriding            # plate_id (simple) | overriding
  collision_weld: on             # off | on   (defaults decided when each is added)
params:
  rift_strength_continent: 0.5   # σc, units of full slab pull per unit length (R2 sweep 0.25/0.5/1.0)
  oceanic_strength_ratio: 3
  oceanic_strength_floor: 0.2
  interior_strength_factor: 2
  interior_strength_km: 1000
  craton_core_factor: 3
  craton_cores: null             # optional polygons file
  suture_strength_factor: 0.5
  suture_width_deg: 2
  sutures: null                  # optional initial sutures file
  belt_amplitude: 0.4            # 0 = off
  belt_scale_km: 1500
  cmb_weakening: 0.5
  rift_persist_myr: 3
  rift_phase_myr: 20
  rift_sweep_azimuths: 12
  rift_sweep_offsets_deg: [-15, 0, 15]
  rift_path_corridor_deg: 5
  weld_after_myr: 30
```

The `rifting: none` path keeps today's geometric neck fragmentation, rewritten without shapely
(Q14).

## Stages: start plain, add one thing at a time

The interview allowed a lot of options. **Build the plainest version that can rift first,
see what it does, then add the rest one at a time**, each with its own gate showing what it
changes. Every added feature stays optional.

Each stage ends with a gate: figures in `show-me/`, led by the map sequence, and a `simple`
regression check (bit-identical, except at R1's re-baseline). Runs are 500 Myr on atlantis66.

**Core (in order):**

| Stage | Content | Gate |
|---|---|---|
| R0 | Tracer thinning, so 500-Myr runs stay affordable. **Done:** `tracer_thinning: cell` (physical); tracers level off at ~70k instead of 200k at 500 Myr; tracer step 0.57 → 0.30 s late in a run; 1° age grid changes by < 1 Myr on average, but steps between young and old seafloor on one plate move by up to 0.5°. `show-me/2026-09-28-r0-tracer-thinning/` | physical 500 Myr: step time, tracer count, forces unchanged |
| R1 | Neck search and `resolve_overlaps(plate_id)` on the sphere, without shapely; re-baseline `simple`. **Done:** overlaps are resolved on the raster (highest ID burnt first), so shapely is gone from that step; `simple` is bit-identical, so there was no re-baseline. The neck search already worked on the raster. New `plate_rebuild` mode: `midpoint` (physical) removes the snapping to meridians and parallels (alignment 1.2 against 7.3); `rotated` was rejected because it makes plate poles jitter. `show-me/2026-09-28-r1-plate-rebuild/` | 3-seed 500-Myr physical runs per mode |
| R2 | **Plain trial split.** `sweep` candidates only, split along the **straight** great circle. Strength is one constant σc on continents; ocean is 3 σc × max(0.2, √(min(age,80)/80)), the one contrast asked for (easy to set to a constant). Persistence 3 Myr. No weakening, patterns, sutures, welding, rift label or path search. Overlaps stay `plate_id`. **Sweep done** (tie seed 0): σc 1.0 gives 1 rift in 500 Myr; σc 0.25 runs away (530 rifts, 182 plates); **σc 0.5 chosen as default for now** (to be tuned): 24 rifts, 222 rifts twice (275 and 406 Myr, both children alive at 500), largest continent 52% → 28%; 15 of 24 rifts are oceanic slivers of 0.1–1% of the sphere living 25–50 Myr; `rift_test` 0.2 s per step. `show-me/2026-09-28-r2-trial-split/` | first rifts; do new plates survive; cost per step; a quick σc sweep (0.25, 0.5, 1.0) |

**Then, one at a time, in an order decided by what R2 shows** (each is its own stage and gate):

- `overlap: overriding`, if R2's new plates are eaten at overlaps
- least-cost split path, instead of the straight cut
- `necks` as candidates, and ribbons
- `rift` label on young rifts (appearance only)
- suture tracers, craton-ID lines and `weak_zones` candidates
- `collision_weld`
- continental strength patterns, one at a time: sutures, interiors, craton cores, belts
- `rift_weakening: cmb_heat`
- the `stress` trigger prototype, a side experiment
- README, plan and config brought up to date

**Gate measures (Q12):**
- rifts per 100 Myr, and new plates surviving more than 50 Myr;
- the largest connected continent over time: continental cells joined across plate
  boundaries, as a share of all continent, which tracks break-up and re-assembly;
- RMS speed (it should stop drifting down), plate count and pole swings;
- tracer count and step time.

## Found in R0

- **RMS speed does not keep falling in 500-Myr physical runs.** It bottoms out at 27–32 km/Myr
  around 220–245 Myr and recovers to 42–43 km/Myr by 480 Myr, as continents 111, 333 and 444
  rearrange. Plate 222 is still whole at 500 Myr: that remains the rifting target.
- **Mixed seafloor ages on one plate** near the south pole on plate 114 (tracers of 20–190 Myr
  within 1°). This is pre-existing and not caused by thinning; left for later.

- **Boundaries snap to meridians and parallels** (`show-me/2026-09-28-boundary-orientation/`).
  Every step's raster rebuild (0.25° lon/lat raster, `polygonize`, then `toposimplify(0.4)` in
  degrees) snaps near-grid-aligned boundaries to exact meridians and parallels. 7–9× more
  boundary length lies within 1° of N–S or E–W than an isotropic set would have in physical, and
  4.4–4.8× in simple (slower plates stay snapped). The input topology has 1.8×. A no-force test
  with rigid rotations reproduces it, so it is the method. Thinning is not measurably worse
  (7.7 against 7.3; three physical runs without thinning give 7.3–9.4).
- **The body ridge-push force points mostly exactly N–S or E–W** (10–11× isotropic), with or
  without thinning. The age grid is piecewise constant (nearest tracer) and its gradient is taken
  along lon and lat. The effect on plate motion is not tested.
- The snapping is fixed in R1 (`plate_rebuild: midpoint`). The ridge-push direction is still
  open.
- **R1 found:**
  - `midpoint` boundaries are about 15% shorter in total, and there are slightly fewer plates
    (mean 8.7–9.5 against 10.1–10.9).
  - One `midpoint` run's alignment rose to about 4 in its last 40 Myr. Watch it in later
    gates.
  - Per-plate Forsyth–Uyeda correlations must leave out short-lived slivers (count plates
    living at least 50 Myr).

- **R2 found:**
  - Young ocean is weaker than continent under the agreed law (3σc√(age/80) < σc below about
    9 Myr; 0.6σc at the floor), so most rifts cut ridge flanks. This is kept for now (plain).
  - With σc 0.5, the continental correlation of Forsyth–Uyeda weakens (r −0.87 → −0.67).
  - σc 1.0 follows the no-rifting run almost exactly, so the trial mode leaves an unsplit
    world unchanged.
  - σc 0.5 over tie seeds 0, 1 and 2: 24, 28 and 21 rifts. 222 rifts at about 275 Myr in seeds 0
    and 1 but never in seed 2, so the supercontinent break-up is not yet robust. The early rifts
    (333, 115, 666, 555) repeat in all three.

## Open

- **Other shapely uses** (`update_plate_polygons` smoothing with topojson, `make_valid`, and
  `contains_xy`) are out of scope for v2 unless they block the rifting work.
- **Per-segment polarity** and **ridge/transform splitting** remain later items.
