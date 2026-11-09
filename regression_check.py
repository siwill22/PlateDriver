"""Compare two PlateDriver output folders file by file.

    python regression_check.py <reference_dir> <new_dir> [--tol 0]

Grids written with grid_output: compact (float32) against a full (float64) reference differ by
float32 rounding: up to about 5e-4 in topography (m), so compare them with --tol 1e-3.

Every output that exists in both folders is compared: NetCDF grids (all variables,
NaN-aware), CSVs, GeoJSONs (attributes and vertex coordinates) and rotation files.
Plots and run_info.yaml are skipped, and so are the wall-clock columns (seconds_*) of
the statistics CSVs.  Prints one line per file and exits with status 1
if anything differs by more than --tol (default 0: exact).

Refuses (status 2) a folder whose run_info.yaml records a slablet with uncommitted changes:
agreement with such a run could not be reproduced from a slablet commit.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd

SKIP_EXT = {'.png', '.yaml', '.log'}


def dirty_slablet(folder):
    """The slablet commit a run recorded, if it had uncommitted changes; otherwise None."""
    path = os.path.join(folder, 'run_info.yaml')
    if not os.path.exists(path):
        return None
    import yaml
    with open(path) as fh:
        info = yaml.safe_load(fh) or {}
    git = (info.get('slablet') or {}).get('git', '') if isinstance(info.get('slablet'), dict) else ''
    return git if str(git).endswith('-dirty') else None


def _max_diff(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return np.inf
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if (nan_a != nan_b).any():
        return np.inf
    return float(np.max(np.abs(a[~nan_a] - b[~nan_a]))) if (~nan_a).any() else 0.0


def _compare_frames(a, b):
    """Max difference between two DataFrames with the same columns (object columns exact)."""
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        return np.inf, 'columns/rows differ: {} x {} vs {} x {}'.format(
            len(a), list(a.columns), len(b), list(b.columns))
    worst = 0.0
    for col in a.columns:
        if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
            worst = max(worst, _max_diff(a[col], b[col]))
        elif not a[col].astype(str).equals(b[col].astype(str)):
            return np.inf, 'column {!r} differs'.format(col)
    return worst, ''


def compare_file(path_a, path_b):
    ext = os.path.splitext(path_a)[1]
    if ext == '.nc':
        da, db = xr.open_dataset(path_a), xr.open_dataset(path_b)
        if set(da.data_vars) != set(db.data_vars):
            return np.inf, 'variables differ'
        return max((_max_diff(da[v].values, db[v].values) for v in da.data_vars), default=0.0), ''
    if ext == '.csv':
        a, b = pd.read_csv(path_a), pd.read_csv(path_b)
        timing = [c for c in a.columns if c.startswith('seconds_')]
        return _compare_frames(a.drop(columns=timing), b.drop(columns=timing, errors='ignore'))
    if ext == '.geojson':
        ga, gb = gpd.read_file(path_a), gpd.read_file(path_b)
        worst, note = _compare_frames(pd.DataFrame(ga.drop(columns='geometry')),
                                      pd.DataFrame(gb.drop(columns='geometry')))
        if note:
            return worst, note
        ca, cb = ga.geometry.get_coordinates(), gb.geometry.get_coordinates()
        return max(worst, _max_diff(ca.values, cb.values)), ''
    if ext == '.rot':
        la, lb = open(path_a).read().split(), open(path_b).read().split()
        if len(la) != len(lb):
            return np.inf, 'different number of entries'
        worst = 0.0
        for x, y in zip(la, lb):
            try:
                worst = max(worst, abs(float(x) - float(y)))
            except ValueError:
                if x != y:
                    return np.inf, 'text differs: {!r} vs {!r}'.format(x, y)
        return worst, ''
    return 0.0, 'not compared'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('reference_dir')
    parser.add_argument('new_dir')
    parser.add_argument('--tol', type=float, default=0.0)
    args = parser.parse_args()

    dirty = {d: dirty_slablet(d) for d in (args.reference_dir, args.new_dir)}
    dirty = {d: g for d, g in dirty.items() if g}
    if dirty:
        for d, g in dirty.items():
            print('refused: {} ran on slablet {} (uncommitted changes)'.format(d, g))
        print('commit slablet and rerun; see decision 7 of slablet/docs/plans/2026-11-08-cross-repo-reconciliation.md')
        sys.exit(2)

    names_a = {f for f in os.listdir(args.reference_dir) if os.path.splitext(f)[1] not in SKIP_EXT}
    names_b = {f for f in os.listdir(args.new_dir) if os.path.splitext(f)[1] not in SKIP_EXT}
    common = sorted(names_a & names_b)
    failed = []
    for name in common:
        worst, note = compare_file(os.path.join(args.reference_dir, name),
                                   os.path.join(args.new_dir, name))
        status = 'OK  ' if worst <= args.tol else 'DIFF'
        if worst > args.tol:
            failed.append(name)
        print('{} {:32s} max diff {:<10.3g} {}'.format(status, name, worst, note))
    only_a, only_b = sorted(names_a - names_b), sorted(names_b - names_a)
    if only_a:
        print('only in reference:', ', '.join(only_a))
    if only_b:
        print('only in new:      ', ', '.join(only_b))
    print('\n{} of {} common files differ'.format(len(failed), len(common)))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
