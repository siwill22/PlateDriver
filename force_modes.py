"""Force-model presets and mode flags.

config.yaml chooses a preset with ``force_model: simple | physical``; any key under
``force_modes`` overrides the preset for that one component.  ``simple`` is the
original model (boundary forces only, torque = angular velocity); ``physical`` is the
v1 force model described in force_model_v1_plan.md.
"""

MODE_CHOICES = {
    'drag':                ('none', 'area', 'area_edge'),
    'ridge_push':          ('boundary', 'body'),
    'slab_pull':           ('hard', 'ramp', 'slab_length'),
    'slab_suction':        ('ratio', 'depth_split'),
    'polarity':            ('id', 'age'),
    'collision':           ('resistance', 'viscous', 'coupling', 'breakoff'),
    'speed_normalisation': ('none', 'calibrate_once', 'per_step'),
    'tracer_thinning':     ('none', 'cell'),
    'plate_rebuild':       ('lonlat', 'midpoint', 'rotated'),
    'rifting':             ('none', 'trial'),
    'collision_weld':      ('off', 'on'),
    'sutures':             ('none', 'tracers'),
}

PRESETS = {
    'simple': {
        'drag':                'none',
        'ridge_push':          'boundary',
        'slab_pull':           'hard',
        'slab_suction':        'ratio',
        'polarity':            'id',
        'collision':           'resistance',
        'speed_normalisation': 'none',
        'tracer_thinning':     'none',
        'plate_rebuild':       'lonlat',
        'rifting':             'none',
        'collision_weld':      'off',
        'sutures':             'none',
    },
    'physical': {
        'drag':                'area_edge',
        'ridge_push':          'body',
        'slab_pull':           'slab_length',
        'slab_suction':        'depth_split',
        'polarity':            'age',
        'collision':           'breakoff',
        'speed_normalisation': 'calibrate_once',
        'tracer_thinning':     'cell',
        'plate_rebuild':       'midpoint',
        'rifting':             'trial',
        'collision_weld':      'off',     # physical default decided at the R3 gate
        'sutures':             'none',    # physical default decided at the R3 gate
    },
}

# Modes that exist in the code so far; the rest arrive stage by stage (see the plan).
IMPLEMENTED = {
    'drag':                {'none', 'area', 'area_edge'},
    'ridge_push':          {'boundary', 'body'},
    'slab_pull':           {'hard', 'ramp', 'slab_length'},
    'slab_suction':        {'ratio', 'depth_split'},
    'polarity':            {'id', 'age'},
    'collision':           {'resistance', 'viscous', 'coupling', 'breakoff'},
    'speed_normalisation': {'none', 'calibrate_once', 'per_step'},
    'tracer_thinning':     {'none', 'cell'},
    'plate_rebuild':       {'lonlat', 'midpoint', 'rotated'},
    'rifting':             {'none', 'trial'},
    'collision_weld':      {'off', 'on'},
    'sutures':             {'none', 'tracers'},
}

DEFAULT_PRESET = 'physical'


def resolve_force_modes(cfg):
    """Return (preset_name, {component: mode}) for a loaded config dict.

    Order: preset, then the legacy ``params.collision_mode`` (older configs), then
    ``force_modes`` overrides.  Raises SystemExit with a clear message for unknown
    presets, components or modes, and for modes not implemented yet.
    """
    preset = cfg.get('force_model', DEFAULT_PRESET)
    if preset not in PRESETS:
        raise SystemExit('force_model must be one of {}, not {!r}'.format(
            ', '.join(PRESETS), preset))
    modes = dict(PRESETS[preset])

    overrides = dict(cfg.get('force_modes') or {})
    legacy_collision = (cfg.get('params') or {}).get('collision_mode')
    if legacy_collision is not None and 'collision' not in overrides:
        overrides['collision'] = legacy_collision

    for component, mode in overrides.items():
        if component not in MODE_CHOICES:
            raise SystemExit('Unknown force_modes entry {!r}; expected one of: {}'.format(
                component, ', '.join(MODE_CHOICES)))
        if mode not in MODE_CHOICES[component]:
            raise SystemExit('force_modes.{} must be one of {}, not {!r}'.format(
                component, ', '.join(MODE_CHOICES[component]), mode))
        modes[component] = mode

    if modes['collision'] == 'viscous' and modes['drag'] == 'none':
        raise SystemExit('collision: viscous is a drag on colliding plates and needs '
                         'drag: area or area_edge (it has nothing to act against with drag: none).')
    missing = ['{}: {}'.format(c, m) for c, m in modes.items() if m not in IMPLEMENTED[c]]
    if missing:
        raise SystemExit('These force modes are not implemented yet:\n  {}\n'
                         'Use force_model: simple, or override them under force_modes.'.format(
                             '\n  '.join(missing)))
    return preset, modes


def describe(preset, modes):
    """One line per component, for the run log."""
    lines = ['force_model: {}'.format(preset)]
    for component, mode in modes.items():
        tag = '' if mode == PRESETS[preset][component] else '   (override)'
        lines.append('  {:20s} {}{}'.format(component, mode, tag))
    return '\n'.join(lines)
