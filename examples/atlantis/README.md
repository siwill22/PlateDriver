# Atlantis: a starting map

Atlantis is a synthetic Earth-like planet drawn in GPlates for the Southern Ocean Marine Geosciences course at
OUC Qingdao. `config.yaml` and `configs/current_model.yaml` start from it, so PlateDriver runs from a fresh clone.
The files come from the course's working folder (octopus, `Atlantis3/`), with the continent file renamed.

| File | What it holds |
|---|---|
| `atlantis66_topologies.gpml` | the plate boundaries (ridges, trenches with their polarity, transforms) |
| `atlantis66_topological_polygons.gpml` | the plate polygons built from those boundaries |
| `atlantis66_continents.geojson` | the continents, each with the plate it rides on (`PLATEID1`); `.gplates.xml` is GPlates' attribute mapping. In octopus this is `atlantis66_cratons.geojson` |
| `atlantis3_rotations.rot` | the rotation model the map was drawn with; with `reconstruction_time: 0` the plates start where they are drawn |

PlateDriver reads only the starting map. From there the plates move under the forces of the model, so a run is
not a reconstruction of Atlantis's history.
