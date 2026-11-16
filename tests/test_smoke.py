"""Smoke tests: 20 steps of each preset from config.yaml's starting map run to the end and
write their outputs. Skipped where the starting map's files are missing; the physical preset
is skipped without slablet, and simple then runs without the mantle outputs. About 2 min in all."""
import os
import subprocess
import sys

import pandas as pd
import pytest
import yaml

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
STEPS = 20


def _config(tmp_path, preset):
    with open(os.path.join(REPO, 'config.yaml')) as fh:
        cfg = yaml.safe_load(fh)
    files = cfg['input_files']
    needed = [files['rotation_model'], files.get('continents', files.get('cratons')), *files['dynamic_polygons']]
    missing = [f for f in needed if not os.path.exists(os.path.join(REPO, os.path.expanduser(f)))]
    if missing:
        pytest.skip('starting map not found: {}'.format(', '.join(missing)))
    cfg['force_model'] = preset
    cfg['output_dir'] = str(tmp_path / preset)
    cfg['run']['total_time'] = float(STEPS)
    cfg['params']['plot_interval'] = None
    cfg['params']['csv_save_interval'] = 10
    path = tmp_path / '{}.yaml'.format(preset)
    with open(path, 'w') as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    return path, tmp_path / preset


@pytest.mark.parametrize('preset', ['simple', 'physical'])
def test_a_short_run_finishes_and_writes_its_outputs(tmp_path, preset):
    if preset == 'physical':
        pytest.importorskip('slablet')
    config, out = _config(tmp_path, preset)
    run = subprocess.run([sys.executable, os.path.join(REPO, 'run_SDP.py'), '--config', str(config)],
                         cwd=REPO, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    steps = pd.read_csv(out / 'step_stats.csv')
    assert len(steps) == STEPS + 1
    assert steps['rms_speed_km_per_myr'].between(1, 500).all()
    for name in ('plate_id', 'seafloor_age', 'topo_bathy'):
        assert (out / '{}_{:04d}.nc'.format(name, STEPS)).exists()
    assert (out / 'rotations_{:04d}.rot'.format(STEPS)).exists()
