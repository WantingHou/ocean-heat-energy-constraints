"""Validate a curated manuscript repository without modifying its files."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re


def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def cells(path):
    with path.open(newline='') as f:
        return list(csv.reader(f))[1:]


def check(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repository', type=Path, required=True)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    root = args.repository.resolve()
    if args.output and args.output.exists():
        raise FileExistsError(args.output)
    checks = {}
    listed = []
    for line in (root / 'metadata/files.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        target = root / name
        check(target.resolve().is_relative_to(root), f'Outside repository: {name}')
        check(target.is_file(), f'Missing file: {name}')
        check(hashlib.sha256(target.read_bytes()).hexdigest() == digest, f'Hash mismatch: {name}')
        listed.append(name)
    checks['base_manifest_files_verified'] = len(listed)
    checks['files_over_100_MiB'] = [str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.stat().st_size >= 100 * 1024**2]
    check(not checks['files_over_100_MiB'], 'GitHub file-size limit exceeded')
    forbidden = [str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and (':Zone.Identifier' in f.name or f.suffix in {'.aux', '.log', '.out'})]
    check(not forbidden, f'Unwanted generated/private files: {forbidden}')
    secret_hits = []
    pattern = re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|AKIA[A-Z0-9]{16}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----)')
    for f in root.rglob('*'):
        if f.is_file() and f.suffix.lower() in {'.py', '.md', '.csv', '.json', '.yml', '.yaml', '.tex', '.txt'}:
            if pattern.search(f.read_bytes()):
                secret_hits.append(str(f.relative_to(root)))
    check(not secret_hits, f'Potential secret material in {secret_hits}; values not printed')
    checks['known_secret_patterns_found'] = len(secret_hits)
    checks['secret_scan_scope'] = 'Common token/private-key signatures in text; not a guarantee against all sensitive content.'
    crosswalk = rows(root / 'metadata/table_crosswalk.csv')
    check(len(crosswalk) == 57, 'Expected 57 archived tables')
    check({int(r['previous_si_number']) for r in crosswalk} == set(range(1, 58)), 'Missing archive table IDs')
    for old, new in [(31, 2), (53, 3), (50, 4)]:
        check(cells(root / f'tables/archived_si/Archive_Table_{old:02d}.csv') == cells(root / f'tables/current_si/SI_Table_{new:02d}.csv'), f'Changed retained SI table {new}')
    checks['archived_tables'] = len(crosswalk)
    checks['current_table_rows'] = {str(n): len(cells(root / f'tables/current_si/SI_Table_{n:02d}.csv')) for n in range(1, 5)}
    check(checks['current_table_rows'] == {'1': 8, '2': 18, '3': 9, '4': 9}, 'Unexpected current SI table shape')
    data = root / 'Source_Data'
    allocations = {'Power': 'geological_power', 'Area': 'geological_area'}
    operators = {'Complete retention': 'direct', '100-year loss': 'loss_100yr', '30-year loss': 'loss_30yr'}
    support = {(r['allocation'], r['family'], r['loss']): r for r in rows(data / 'panels/F3d_condition_memberships.csv')}
    for row in cells(root / 'tables/current_si/SI_Table_02.csv'):
        allocation, family, operator, count, arange, brange, ids = row
        original = support[allocations[allocation], family, operators[operator]]
        check(int(count) == int(original['count']), f'Support count: {row}')
        if int(count):
            expected_a = f"{float(original['A_over_O_low']):.2f} to {float(original['A_over_O_high']):.2f}"
            expected_b = f"{float(original['B_over_O_low']):.2f} to {float(original['B_over_O_high']):.2f}"
            check(arange == expected_a and brange == expected_b and ids == original['display_ids'], f'Support endpoints: {row}')
        else:
            check(arange == brange == ids == 'Empty', 'Empty support has been changed to a number')
    checks['current_support_conditions_verified'] = len(support)
    paired = rows(data / 'paired_current_joint_support.csv')
    key = {(r['allocation'], int(r['candidate_index'])): r for r in paired}
    identities = {int(r['candidate_index']): r['display_id'] for r in rows(data / 'panels/F3c_primary_budgets.csv')}
    for row in cells(root / 'tables/current_si/SI_Table_04.csv'):
        allocation, candidate, observation, model, member, sign, value, label = row
        original = key[allocations[allocation], int(candidate)]
        check(observation == original['O_member'] and model == original['model'] and member == original['internal_member'], 'Changed observation/climate identity')
        check(sign == original['paired_extension_sign'] and label == identities[int(candidate)], 'Changed sign or budget label')
        check(abs(float(value) - float(original['B_1955_2023_ZJ'])) <= 0.0000005, 'Rounded B differs from underlying source')
    residual = 0.0
    for row in paired:
        for first, last in [(1955, 2014), (1955, 2023), (1992, 2023)]:
            suffix = f'_{first}_{last}_ZJ'
            residual = max(residual, abs(float(row['O' + suffix]) - sum(float(row[t + suffix]) for t in ['A', 'N', 'B', 'E'])))
    check(residual < 1e-10, f'Budget closure error {residual}')
    checks.update(primary_allocation_memberships=len(paired), unique_budget_identities=len(identities), max_budget_closure_error_ZJ=residual)
    regional = {(float(r['age_Ma']), r['summary']): r for r in rows(data / 'panels/F2c_paired_heat_flow.csv')}
    summary = {'Station mean': 'primary_station', 'Equal profiles': 'primary_profile_equal', 'Model': 'fixed_age_prediction'}
    # Resolve the model label from the source records instead of changing source data.
    model_labels = {r['summary'] for r in regional.values()} - {'primary_station', 'primary_profile_equal'}
    check(len(model_labels) == 1, 'Unexpected regional model summaries')
    summary['Model'] = model_labels.pop()
    regional_error = 0.0
    for row in cells(root / 'tables/current_si/SI_Table_03.csv'):
        age, label, conduction, deficit, total = row
        original = regional[float(age), summary[label]]
        for value, col in [(conduction, 'conduction_mW_m2'), (deficit, 'deficit_mW_m2'), (total, 'cooling_reference_mW_m2')]:
            check(abs(float(value) - float(original[col])) <= 0.0000005, 'Regional rounded value mismatch')
        actual = float(original['conduction_mW_m2']) + float(original['deficit_mW_m2'])
        regional_error = max(regional_error, abs(actual - 506.7 / math.sqrt(float(age))))
    check(regional_error < 1e-10, 'Regional budget mismatch')
    checks['regional_budget_rows_verified'] = len(regional)
    checks['max_regional_cooling_identity_error_mW_m2'] = regional_error
    figures = rows(root / 'metadata/figure_crosswalk.csv')
    check(len(figures) == 12, 'Expected twelve current figure references')
    for figure in figures:
        for path in figure['primary_data'].split(';'):
            check((data / path).exists(), f'Missing plotted data {path}')
    checks['current_figure_source_mappings_verified'] = len(figures)
    checks['inherited_display_metadata_is_historical'] = True
    checks['public_upload_performed'] = False
    checks['passed'] = True
    content = json.dumps(checks, indent=2) + '\n'
    if args.output:
        with args.output.open('x') as f:
            f.write(content)
    print(content)


if __name__ == '__main__':
    main()
