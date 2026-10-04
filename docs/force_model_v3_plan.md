# Force model v3 — plan: new trenches, and recycling old ocean

Draft prepared with Claude Code, 2026-10-04. Status: **agreed scope** (decisions D1–D5 below);
stages still start plain and are gated one at a time.

## Where v2 left the model

v2 (closed 2026-10-04, `force_model_v2_plan.md`) made continents rift, stopped colliding
continents overlapping, and kept the plate count in check. The current model is
`configs/current_model.yaml`; three 1000-Myr runs are in the viewer and in
`show-me/2026-10-04-current-model/`.

What v2 cannot do: **start a new subduction zone**. Every trench in a run descends from the
starting map or from a polarity rule at a new contact. Two measures show the cost:

| Current model, 1000 Myr | Tie seed 0 | Tie seed 1 | Tie seed 2 | Earth today |
|---|---|---|---|---|
| Ocean older than 200 Myr | 25% | 28% | 25% | 0% |
| Oldest ocean | 1016 Myr | 1085 Myr | 1054 Myr | about 180 Myr (Mediterranean relics older) |
| Median ocean age | 126 Myr | 67 Myr | 71 Myr | about 60 Myr |
| Largest plate, max | 37% | 37% | 49% | Pacific, about 20% |

The Earth column is for orientation only, **not a target**: the oldest ocean today is about
180 Myr, but it may well have been much older at other times. What matters is that ocean is
never recycled unless it happens to sit next to a trench, so ocean of the starting map survives
for the whole run. It is also the strongest material in the rift test (strength 1.5, against
tensions of 0.2–0.4), so the plates made of it never break either. On Earth new trenches
start mostly at weaknesses near existing boundaries: passive margins, fracture zones and
transforms, often pushed by compression from elsewhere, or helped by a plume (Stern 2004;
Gurnis et al. 2004; Gerya et al. 2015). Old Atlantic-type margins are too strong to fail under
their own weight alone (Cloetingh et al. 1982), so something must weaken or push them.

## Goal

Add **subduction initiation**: a plate is cut along a weak line and one side starts to
subduct beneath the other. Like rifting, it is a mode flag that is off by default, and the
`simple` preset is untouched. Judge it on whole runs: whether old ocean can be recycled at all, the
size of the largest plate, the number of trenches and where they start. Present-day Earth
values are a reference, not something to tune towards.

## Stages: start plain, add one thing at a time

Every stage ends with a gate: a three-seed 500-Myr run of the current model plus the change
(D3), figures in `show-me/` led by the map sequence, and the gate measures below. The phase
closes with one three-seed 1000-Myr check.

| Stage | Content | Gate |
|---|---|---|
| S0 | **Gate tooling.** One command runs the three seeds of a config, writes the standard table (plates, slivers, largest plate, RMS, rifts, welds, hidden continent, ocean-age distribution), the map sequence, and exports the runs to the viewer. Runs are kept outside the scratchpad (a `runs/` directory in `config`, gitignored). | the current model reproduces its numbers |
| S1 | **Plain subduction initiation at passive margins** (D1, D2). Candidate lines: the continent–ocean edges inside one plate, where the ocean along the edge is older than `init_min_age_myr`. Trigger: a compression trial on the same footing as the rift test (the torque split across the candidate line, compression against a strength in compression); the line must pass on a similar position for `rift_persist_myr` steps. Polarity: the ocean subducts beneath the continent. | ocean-age distribution; largest plate; where and how often trenches start; plate count stays 10–25 |
| S2 | **Fracture zones and old ocean.** Add age steps inside a plate (fracture zones, ridge jumps) as candidate lines; the older side subducts. | as S1; are intra-oceanic arcs made? |
| S3 | **Back-arc rifting.** A weak zone on the overriding plate behind a retreating trench (`strength × factor` within a distance behind it), so rifts open there. | back-arc basins form behind rolling-back trenches and close again |
| S4 | **Plume weakening** (`rift_weakening: cmb_heat`, planned in v2 and not done): lower strength above hot spots of the CMB heat field, for rifting and initiation. | breakups near plume heads; does it change supercontinent lifetimes |

**Then the teaching outputs (D4)**, each its own small stage:

- **Teaching outputs:** ridges split into spreading segments and transforms, as a post-process
  of a chosen snapshot; viewer V2 grids (seafloor age, topography, divergence, plate ID) with a
  hover readout; a run info panel.

## Gate measures

- **Ocean-age distribution:** share of ocean older than 100, 150 and 200 Myr; oldest ocean;
  median; how much of the starting map's ocean survives. A diagnostic of recycling, not a
  target: Earth today (median about 60 Myr, oldest about 180 Myr) is one snapshot.
- **Largest plate over time**, and the plate-size distribution (slivers under 0.5% separately).
- **Trench starts:** number per 100 Myr, the setting (passive margin, fracture zone, other) and
  the age of the ocean at the start.
- Trench and ridge length; RMS speed and net rotation; continental rifts; hidden continent.
- Step time.

## Decisions (2026-10-04)

| # | Topic | Decision |
|---|---|---|
| D1 | S1 trigger | **Hybrid:** a rule picks the candidates (passive margins with ocean older than `init_min_age_myr`), the compression trial decides which fail |
| D2 | S1 location | **Passive margins only**; fracture zones and old ocean are S2 |
| D3 | Gate length | **500 Myr × three seeds** per stage; one 1000-Myr three-seed check when the phase closes |
| D4 | Scope | Subduction initiation (S0–S2), then **back-arc rifting and plume weakening** (S3, S4), then the **teaching outputs**. The v2 known issues are not in this phase |
| D5 | Earth values | Present-day Earth (ocean ages, plate sizes) is a reference, **not a target** to tune towards |

## Open (not in this phase)

- **v2 known issues:** hidden continent lost at a rift (rebuild continents by cutting the
  polygons, not from the raster); calibrating the stiff contact (the convergence left; the
  100-km ocean gap); per-segment trench polarity.

- Other shapely uses (`update_plate_polygons` smoothing with topojson, `make_valid`,
  `contains_xy`).
- `necks` as rift candidates, the `rift` label, craton-ID lines, the `stress` prototype and
  `overlap: overriding` (v2 items not done).
