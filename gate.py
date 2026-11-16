"""Run one gate: a config over several tie seeds, the standard table, the map sequence, and
(optionally) the viewer export.

    conda run -n pygmt17 python gate.py configs/current_model.yaml --name current_model \\
        [--seeds 0 1 2] [--until 500] [--export] [--replace-viewer]

Steps (each can be skipped):
  1. run   runs/<name>/tbN/ for every seed N, all seeds in parallel (run_SDP.py with the
           config, tie_break_seed N, total_time --until, output_dir runs/<name>/tbN).  A seed
           whose run already reached --until is not run again (--rerun forces it).
  2. table runs/<name>/gate_table.csv: one row per seed and period (--period, Myr) and one
           for the whole run, with the measures below; printed as a markdown table.
  3. figs  show-me/<date>-<name>/: diagnose_run.py's map sequence (one column per seed, every
           --until/10 Myr), time series and Forsyth–Uyeda plots, the table and a README.
  4. export (--export) the runs to the viewer as "<label>, tie seed N"; --replace-viewer first
           removes every other run from the viewer.

Measures (whole-run values are over all steps; area shares are of the sphere):
  plates, slivers (< 0.5%), plates > 2%      mean count per step
  largest plate                               mean and max share
  RMS speed                                   mean, km/Myr
  rifts, continental rifts, welds, accretions events (a continental rift leaves > 0.2% of
                                              the sphere of continent on both sides)
  hidden continent                            max share of continent hidden by overlap
  trench, ridge                               mean length, 1000 km
  ocean age at the period's end (saved step)  median, share older than 100 / 150 / 200 Myr,
                                              oldest, and the share that already existed at
                                              t = 0 (older than the elapsed time)
  step time                                   mean seconds per step
The gate stops if slablet has uncommitted changes (before running), or if any seed's run_info.yaml
records such a slablet (with --no-run): a gate result must be reproducible from commits.

Everything is from the run's own outputs (step_stats.csv, plate_stats.csv, rift_events.csv,
welds.csv, seafloor_age_NNNN.nc, continent_NNNN.nc or, before Nov 2026, craton_NNNN.nc); ocean areas are cos(lat)-weighted.
"""
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import xarray as xr
import yaml

REPO = os.path.dirname(os.path.abspath(__file__))
SLIVER = 0.005
BIG = 0.02
CONT_RIFT = 0.002


def seed_dir(runs_dir, name, seed):
    return os.path.join(runs_dir, name, 'tb%d' % seed)


def finished(d, until):
    f = os.path.join(d, 'step_stats.csv')
    if not os.path.exists(f):
        return False
    try:
        return float(pd.read_csv(f).time_myr.iloc[-1]) >= until
    except Exception:
        return False


def run_seeds(cfg_path, runs_dir, name, seeds, until, rerun):
    cfg = yaml.safe_load(open(cfg_path))
    procs = []
    for s in seeds:
        d = seed_dir(runs_dir, name, s)
        if finished(d, until) and not rerun:
            print('tb%d: already at %g Myr, not run again' % (s, until))
            continue
        os.makedirs(d, exist_ok=True)
        c = json.loads(json.dumps(cfg))           # deep copy
        c['output_dir'] = d
        c['params']['tie_break_seed'] = s
        c['run']['total_time'] = float(until)
        cpath = os.path.join(os.path.dirname(d), 'tb%d.yaml' % s)
        yaml.safe_dump(c, open(cpath, 'w'), sort_keys=False)
        log = open(os.path.join(os.path.dirname(d), 'tb%d.log' % s), 'w')
        print('tb%d: running to %g Myr (log %s)' % (s, until, log.name))
        procs.append((s, subprocess.Popen([sys.executable, os.path.join(REPO, 'run_SDP.py'), '--config', cpath],
                                          cwd=REPO, stdout=log, stderr=subprocess.STDOUT)))
    for s, p in procs:
        p.wait()
        print('tb%d: finished (exit %d)' % (s, p.returncode))
    bad = [s for s, p in procs if p.returncode != 0]
    if bad:
        raise SystemExit('runs failed for seeds %s; see the logs in %s' % (bad, os.path.join(runs_dir, name)))


def ocean_age(d, step, t):
    """Ocean-age measures at a saved step at model time t (None if the grids are missing)."""
    fa = os.path.join(d, 'seafloor_age_%04d.nc' % step)
    fc = [os.path.join(d, '%s_%04d.nc' % (n, step)) for n in ('continent', 'craton')]
    fc = [(f, n) for f, n in zip(fc, ('continent', 'craton')) if os.path.exists(f)]
    if not (os.path.exists(fa) and fc):
        return None
    a = xr.open_dataset(fa).z.values.astype(float)
    c = np.nan_to_num(xr.open_dataset(fc[0][0])[fc[0][1]].values.astype(float)) > 0
    lat = np.linspace(-90, 90, a.shape[0])
    w = np.cos(np.radians(lat))[:, None] * np.ones(a.shape)
    oc = ~c & np.isfinite(a)
    W = w[oc].sum()
    if W <= 0:
        return None
    order = np.argsort(a[oc])
    cw = np.cumsum(w[oc][order]) / W
    older = lambda x: 100 * w[oc & (a > x)].sum() / W
    return dict(median=float(a[oc][order][np.searchsorted(cw, 0.5)]), gt100=older(100), gt150=older(150),
                gt200=older(200), oldest=float(np.nanmax(a[oc])), from_t0=older(t))


def summarise(d, t0, t1):
    st = pd.read_csv(os.path.join(d, 'step_stats.csv'))
    pl = pd.read_csv(os.path.join(d, 'plate_stats.csv'))
    rf_path, wd_path = os.path.join(d, 'rift_events.csv'), os.path.join(d, 'welds.csv')
    rf = pd.read_csv(rf_path) if os.path.exists(rf_path) and os.path.getsize(rf_path) else pd.DataFrame(columns=['time_myr', 'parent', 'child'])
    wd = pd.read_csv(wd_path) if os.path.exists(wd_path) and os.path.getsize(wd_path) else pd.DataFrame(columns=['time_myr', 'accreted'])
    m = (st.time_myr >= t0) & (st.time_myr <= t1)
    s = st[m]
    p = pl[(pl.time_myr >= t0) & (pl.time_myr <= t1)]
    times = s.time_myr.values
    count = lambda q: q.groupby('time_myr').size().reindex(times, fill_value=0)
    big = p.groupby('time_myr').area_frac.max().reindex(times)
    hid = 100 * s.continent_hidden_km2 / (s.continent_area_km2 + s.continent_hidden_km2) \
        if 'continent_hidden_km2' in s else pd.Series(np.nan, index=s.index)
    pl['cont'] = pl.area_frac * pl.continental_frac
    cont = pl.set_index(['time_myr', 'plate_id']).cont
    r = rf[(rf.time_myr >= t0) & (rf.time_myr <= t1)]
    n_cont = sum(min(cont.get((t, a), 0), cont.get((t, b), 0)) > CONT_RIFT for t, a, b in zip(r.time_myr, r.parent, r.child))
    w = wd[(wd.time_myr >= t0) & (wd.time_myr <= t1)]
    row = {
        'plates': s.n_plates.mean(), 'slivers': count(p[p.area_frac < SLIVER]).mean(),
        'plates_gt2pct': count(p[p.area_frac >= BIG]).mean(),
        'largest_mean_pct': 100 * big.mean(), 'largest_max_pct': 100 * big.max(),
        'rms_km_myr': s.rms_speed_km_per_myr.mean(),
        'rifts': len(r), 'continental_rifts': int(n_cont), 'welds': len(w),
        'accretions': int(w['accreted'].astype(bool).sum()) if 'accreted' in w else 0,
        'hidden_max_pct': hid.max(),
        'trench_1000km': s.trench_km.mean() / 1e3, 'ridge_1000km': s.ridge_km.mean() / 1e3,
        'step_s': s.seconds_step.mean(),
    }
    saved = sorted(int(f[13:17]) for f in os.listdir(d) if f.startswith('seafloor_age_') and f.endswith('.nc'))
    end = max([x for x in saved if x <= t1], default=None)
    age = ocean_age(d, end, float(st.loc[st.step == end, 'time_myr'].iloc[0])) if end is not None else None
    if age:
        row.update({'age_at': end, 'age_median': age['median'], 'age_gt100_pct': age['gt100'],
                    'age_gt150_pct': age['gt150'], 'age_gt200_pct': age['gt200'], 'age_oldest': age['oldest'],
                    'age_from_t0_pct': age['from_t0']})
    return row


def gate_table(runs_dir, name, seeds, until, period):
    rows = []
    edges = list(np.arange(0, until, period)) + [until]
    for s in seeds:
        d = seed_dir(runs_dir, name, s)
        for a, b in zip(edges[:-1], edges[1:]):
            rows.append({'seed': s, 'period': '%g-%g' % (a, b), **summarise(d, a, b)})
        rows.append({'seed': s, 'period': 'whole run', **summarise(d, 0, until)})
    df = pd.DataFrame(rows)
    out = os.path.join(runs_dir, name, 'gate_table.csv')
    df.round(3).to_csv(out, index=False)
    return df, out


def markdown(df):
    show = [('seed', 'Seed', '%d'), ('period', 'Myr', '%s'), ('plates', 'Plates', '%.1f'),
            ('slivers', '< 0.5%', '%.1f'), ('largest_max_pct', 'Largest max %', '%.0f'),
            ('rms_km_myr', 'RMS', '%.0f'), ('rifts', 'Rifts', '%d'), ('continental_rifts', 'Cont. rifts', '%d'),
            ('welds', 'Welds', '%d'), ('hidden_max_pct', 'Hidden max %', '%.1f'),
            ('age_median', 'Ocean age median', '%.0f'), ('age_gt200_pct', '> 200 Myr %', '%.0f'),
            ('age_from_t0_pct', 'From t=0 %', '%.0f'), ('step_s', 's/step', '%.1f')]
    show = [x for x in show if x[0] in df]
    lines = ['| ' + ' | '.join(h for _, h, _ in show) + ' |', '|' + '---|' * len(show)]
    for _, r in df.iterrows():
        cells = []
        for k, _, f in show:
            v = r[k]
            cells.append('–' if v is None or (isinstance(v, float) and not np.isfinite(v)) else f % v)
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def figures(runs_dir, name, seeds, until, cfg_path, df, table_md):
    out = os.path.join(REPO, 'show-me', '%s-%s' % (datetime.date.today().isoformat(), name))
    os.makedirs(out, exist_ok=True)
    dirs = [seed_dir(runs_dir, name, s) for s in seeds]
    step = max(1, int(round(until / 10)))
    seq = [int(x) for x in range(0, int(until) + 1, step)]
    subprocess.run([sys.executable, os.path.join(REPO, 'diagnose_run.py'), *dirs,
                    '--labels', *['tie seed %d' % s for s in seeds], '--out', out,
                    '--sequence-steps', *map(str, seq), '--steps', str(int(until))], cwd=REPO, check=True)
    df.round(3).to_csv(os.path.join(out, 'gate_table.csv'), index=False)
    with open(os.path.join(out, 'README.md'), 'w') as f:
        f.write('# Gate: %s\n\nConfig `%s`, tie seeds %s, %g Myr; runs in `%s`.\n'
                'Made by `gate.py`; `sequence.png` is the main figure.\n\n%s\n'
                % (name, os.path.relpath(cfg_path, REPO), ', '.join(map(str, seeds)), until,
                   os.path.relpath(os.path.join(runs_dir, name), REPO), table_md))
    return out


def export(runs_dir, name, seeds, label, replace):
    data = os.path.join(REPO, 'viewer', 'public', 'data')
    if replace:
        keep = {'%s_tb%d' % (name, s) for s in seeds}
        rdir = os.path.join(data, 'runs')
        for d in os.listdir(rdir) if os.path.isdir(rdir) else []:
            if d not in keep:
                shutil.rmtree(os.path.join(rdir, d))
        with open(os.path.join(data, 'runs.json'), 'w') as f:
            json.dump([], f)
    for s in seeds:
        subprocess.run([sys.executable, os.path.join(REPO, 'viewer', 'prep', 'export_run.py'),
                        seed_dir(runs_dir, name, s), '--name', '%s_tb%d' % (name, s),
                        '--label', '%s, tie seed %d' % (label, s)], cwd=REPO, check=True,
                       stdout=subprocess.DEVNULL)
        print('exported tb%d to the viewer' % s)


def check_slablet_committed(runs_dir, name, seeds, no_run):
    """Stop the gate when its results would rest on uncommitted slablet code."""
    from regression_check import dirty_slablet
    if not no_run:
        try:
            import slablet
        except ImportError:
            return                       # PlateDriver without the mantle model: nothing to check
        git = str(slablet.version_info().get('git', ''))
        if git.endswith('-dirty'):
            raise SystemExit('gate stopped: slablet has uncommitted changes ({}); commit them first'.format(git))
    bad = [s for s in seeds if dirty_slablet(os.path.join(runs_dir, name, 'tb{}'.format(s)))]
    if bad:
        raise SystemExit('gate stopped: seeds {} ran on a slablet with uncommitted changes; rerun them'.format(bad))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('config')
    ap.add_argument('--name', help='gate name (default: the config file name)')
    ap.add_argument('--label', help='label in the viewer (default: the name)')
    ap.add_argument('--seeds', nargs='+', type=int, default=[0, 1, 2])
    ap.add_argument('--until', type=float, default=500.0, help='model time to run to (Myr)')
    ap.add_argument('--period', type=float, default=250.0, help='table period (Myr)')
    ap.add_argument('--runs-dir', default=os.path.join(REPO, 'runs'))
    ap.add_argument('--rerun', action='store_true', help='run seeds again even if finished')
    ap.add_argument('--no-run', action='store_true')
    ap.add_argument('--no-figs', action='store_true')
    ap.add_argument('--export', action='store_true', help='export the runs to the viewer')
    ap.add_argument('--replace-viewer', action='store_true', help='with --export: remove all other viewer runs')
    a = ap.parse_args()
    name = a.name or os.path.splitext(os.path.basename(a.config))[0]
    runs_dir = os.path.abspath(a.runs_dir)
    check_slablet_committed(runs_dir, name, a.seeds, a.no_run)
    if not a.no_run:
        run_seeds(os.path.abspath(a.config), runs_dir, name, a.seeds, a.until, a.rerun)
    df, path = gate_table(runs_dir, name, a.seeds, a.until, a.period)
    md = markdown(df)
    print('\n' + md + '\n\ntable: ' + path)
    if not a.no_figs:
        print('figures:', figures(runs_dir, name, a.seeds, a.until, os.path.abspath(a.config), df, md))
    if a.export:
        export(runs_dir, name, a.seeds, a.label or name, a.replace_viewer)


if __name__ == '__main__':
    main()
