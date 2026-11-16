"""Reproduce frozen analytical results in fresh output directories, without downloads."""
from pathlib import Path
import argparse
import importlib.metadata
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
import xarray as xr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repository', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--log-dir', type=Path, required=True)
    a = ap.parse_args()
    root, out, logs = a.repository.resolve(), a.output_dir.resolve(), a.log_dir.resolve()
    for path in [out, logs]:
        if path.exists():
            raise FileExistsError(path)
        if path.is_relative_to(root):
            raise ValueError('Output and log directories must be outside the immutable repository')
    if out == logs or out.is_relative_to(logs) or logs.is_relative_to(out):
        raise ValueError('Use separate output and log directories')
    out.mkdir(parents=True)
    logs.mkdir(parents=True)
    source = root / 'Source_Data/retained'
    original = root / 'code/original'
    views = out / 'input_views'
    views.mkdir()
    checked = []
    executed = []

    def link(target, destination):
        if not target.exists():
            raise FileNotFoundError(target)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.symlink_to(target, target_is_directory=target.is_dir())

    def run(batch, name, **kwargs):
        program = original / batch / name
        command = [sys.executable, '-I', '-B', str(program)]
        for key, value in kwargs.items():
            command.extend(['--' + key.replace('_', '-'), str(value)])
        log = logs / f'{len(executed) + 1:02d}_{program.stem}.log'
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
        with log.open('x') as f:
            result = subprocess.run(command, cwd=out, env=env, stdout=f, stderr=subprocess.STDOUT)
        executed.append({'script': str(program.relative_to(root)), 'returncode': result.returncode, 'log': str(log)})
        if result.returncode:
            raise RuntimeError(f'{name} failed ({result.returncode}); see {log}')
        print(f'PASS {name}', flush=True)

    def compare_csv(actual, expected):
        left = pd.read_csv(actual, float_precision='round_trip')
        right = pd.read_csv(expected, float_precision='round_trip')
        if list(left.columns) != list(right.columns) or left.shape != right.shape:
            raise ValueError(f'Table shape/column mismatch: {actual.name}')
        left = left.copy()
        right = right.copy()
        for col in left.columns:
            if col.endswith('_json'):
                for lval, rval in zip(left[col], right[col]):
                    np.testing.assert_allclose(json.loads(lval), json.loads(rval), rtol=1e-10, atol=1e-9)
                left[col] = ''
                right[col] = ''
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=False, rtol=1e-10, atol=1e-9)
        checked.append({'file': str(expected.relative_to(root)), 'rows': len(left), 'columns': len(left.columns), 'passed': True})

    def compare_folder(actual, expected):
        for f in sorted(actual.glob('*.csv')):
            compare_csv(f, expected / f.name)

    # A symlink layout satisfies the original programs' interfaces without altering them.
    geology_view = views / 'geological_inputs'
    link(source, geology_view / 'Source_Data')
    certificates = out / 'conditional_budget'
    geo_batch = 'geological_support_20261105T131309Z'
    run(geo_batch, 'certify_geological_budget.py', project_root=geology_view,
        config=source / 'geological_config.json', transfer_config=source / 'regional_transfer_config.json',
        output_dir=certificates)
    compare_folder(certificates, source / 'geological_support/conditional_budget')
    with xr.open_dataset(certificates / 'geology_conditioned_source_fields.nc') as generated:
        with xr.open_dataset(source / 'geological_support/conditional_budget/geology_conditioned_source_fields.nc') as frozen:
            xr.testing.assert_allclose(generated, frozen, rtol=1e-10, atol=1e-9)
    checked.append({'file': 'Source_Data/retained/geological_support/conditional_budget/geology_conditioned_source_fields.nc', 'passed': True})

    physical_view = views / 'physical_inputs'
    for name in ['frozen_source_config.json', 'F3_all_window_operator_pairs.csv']:
        link(source / name, physical_view / 'Source_Data' / name)
    for name in ['certified_histories.csv', 'source_only_bounds.csv']:
        link(certificates / name, physical_view / 'Source_Data' / name)
    link(source / 'regional_transfer_config.json', physical_view / 'regional_transfer_config.json')
    physical = out / 'physical_limits'
    physical_batch = 'physical_limits_20261106T004812Z'
    run(physical_batch, 'analyze_physical_limits.py', source_package=physical_view,
        config=source / 'physical_config.json', output_dir=physical)
    compare_folder(physical, source / 'physical_limits')
    run(physical_batch, 'derive_exposure_frontier.py', results=physical,
        config=source / 'physical_config.json', output_dir=physical / 'exposure')
    compare_folder(physical / 'exposure', source / 'physical_limits/exposure')
    run(physical_batch, 'verify_physical_limits.py', results=physical,
        source_package=physical_view, output=out / 'physical_independent_verification.json')
    run(physical_batch, 'verify_exposure_frontier.py', results=physical,
        output=out / 'exposure_independent_verification.json')

    analytical_view = views / 'analytical_inputs'
    link(source, analytical_view / 'Source_Data')
    evidence = out / 'evidence_closure'
    run('NG_evidence_closure_20261106T033806Z', 'analyze_evidence_dependencies.py',
        package=analytical_view, config=root / 'model_configs/NG_evidence_closure_20261106T033806Z.json',
        output_dir=evidence)
    compare_folder(evidence, source / 'evidence_closure')
    release = out / 'global_release_support'
    release_batch = 'global_release_support_20261106T040811Z'
    run(release_batch, 'analyze_release_support.py', package=analytical_view,
        config=root / 'model_configs/global_release_support_20261106T040811Z.json', output_dir=release)
    compare_folder(release, source / 'global_release_support')
    run(release_batch, 'verify_release_support.py', package=analytical_view, results=release,
        output=out / 'release_independent_verification.json')
    run(release_batch, 'derive_sampling_design.py', thresholds=release / 'heterogeneous_release_thresholds.csv',
        output=release / 'prospective_sampling_design.csv')
    compare_csv(release / 'prospective_sampling_design.csv', source / 'global_release_support/prospective_sampling_design.csv')

    regional = out / 'external_validation'
    run('NG_external_validation_20261106T045511Z', 'analyze_portable_station_table.py',
        stations=source / 'external_validation/stations.csv',
        support=root / 'model_configs/NG_external_validation_20261106T045511Z_station_support.json', output=regional)
    compare_folder(regional, source / 'external_validation')
    independent = {name: json.loads((out / name).read_text()) for name in [
        'physical_independent_verification.json', 'exposure_independent_verification.json', 'release_independent_verification.json']}
    result = {'passed': True, 'executed_scripts': executed, 'frozen_outputs_compared': checked,
              'independent_verification': independent,
              'csv_tolerance': {'rtol': 1e-10, 'atol': 1e-9},
              'python': sys.version, 'packages': {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy', 'xarray']},
              'scope': 'Processed-input analytical reproduction, not raw-data acquisition, climate simulation, pixel-identical figure reproduction or empirical global-source validation.',
              'downloads_performed': False, 'model_simulations_started': False,
              'repository_modified': False, 'public_upload_performed': False}
    with (out / 'reproduction_validation.json').open('x') as f:
        json.dump(result, f, indent=2)
    print(json.dumps({'passed': True, 'programs': len(executed), 'frozen_outputs_compared': len(checked),
                      'report': str(out / 'reproduction_validation.json')}, indent=2))


if __name__ == '__main__':
    main()
