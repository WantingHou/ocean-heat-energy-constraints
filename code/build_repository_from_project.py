"""Create a manuscript-specific, create-only repository candidate from frozen files."""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import re
import sys

BASE = 'docs/manuscript/NG_integrated_editable_figures_20261109T031639Z'
CURRENT = 'docs/manuscript/Wnating_OHC/Wnating_OHC/manuscript'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text_cell(value):
    value = value.replace(r'\allowbreak{}', '').replace(r'\%', '%')
    value = value.replace(r'\_', '_').replace(r'\&', '&')
    value = value.replace(r'\textsuperscript', '^').replace(r'\textsubscript', '_')
    value = re.sub(r'\\(?:textbf|textrm|textit|ensuremath)\{([^{}]*)\}', r'\1', value)
    value = value.replace(r'\times', '×').replace(r'\leq', '≤').replace(r'\geq', '≥')
    return re.sub(r'\s+', ' ', value).strip()


def tables(path):
    source = path.read_text()
    out = []
    for m in re.finditer(r'\\section\*\{Supplementary Table (\d+)\. ([^\n]+)\}', source):
        end = source.index(r'\end{longtable}', m.end()) + len(r'\end{longtable}')
        fragment = source[m.start():end]
        header = fragment.split(r'\toprule', 1)[1].split(r'\midrule\endfirsthead', 1)[0]
        body = fragment.split(r'\bottomrule\endfoot', 1)[1].split(r'\end{longtable}', 1)[0]
        row_strings = [x.strip() for x in re.split(r'\\\\', header + body) if x.strip()]
        raw = [re.split(r'(?<!\\)&', x) for x in row_strings]
        rows = [[text_cell(c) for c in row] for row in raw]
        assert len(rows) > 1 and all(len(row) == len(rows[0]) for row in rows), m.group(1)
        out.append({'number': int(m.group(1)), 'title': m.group(2), 'rows': rows,
                    'fragment': fragment, 'line': source.count('\n', 0, m.start()) + 1})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    root, out = args.project_root.resolve(), args.output_dir.resolve()
    if out.exists():
        raise FileExistsError(out)
    current, base = root / CURRENT, root / BASE
    archived, selected = tables(base / 'Supplementary_Information_integrated_EN.tex'), tables(current / 'SI.tex')
    assert len(archived) == 57 and len(selected) == 4
    out.mkdir(parents=True)
    manifest = []

    def write(name, content):
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('x', encoding='utf-8', newline='') as f:
            f.write(content)

    def write_json(name, content):
        write(name, json.dumps(content, indent=2, ensure_ascii=False) + '\n')

    def write_csv(name, rows):
        buffer = io.StringIO(newline='')
        csv.writer(buffer).writerows(rows)
        write(name, buffer.getvalue())

    def preserve(source, name):
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as f:
            f.write(source.read_bytes())
        assert sha(target) == sha(source)
        manifest.append({'repository_path': name, 'project_source': str(source.relative_to(root)),
                         'bytes': source.stat().st_size, 'sha256': sha(source)})

    for f in sorted((base / 'Source_Data').rglob('*')):
        if f.is_file() and ':Zone.Identifier' not in f.name:
            preserve(f, 'Source_Data/' + f.relative_to(base / 'Source_Data').as_posix())

    crosswalk = [['archive_table_id', 'previous_si_number', 'title', 'current_si_table',
                  'classification', 'csv', 'verbatim_tex', 'source_line']]
    historical = {7, 8, 9, 11, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}
    retained = {31: 2, 53: 3, 50: 4}
    for t in sorted(archived, key=lambda x: x['number']):
        n = t['number']
        prefix = f'tables/archived_si/Archive_Table_{n:02d}'
        write_csv(prefix + '.csv', t['rows'])
        write(prefix + '.tex', '% Verbatim table fragment from preceding SI; not current numbering.\n' + t['fragment'] + '\n')
        state = 'retained_in_current_si' if n in retained else 'historical_condition' if n in historical else 'supporting_detail'
        if n == 54:
            state = 'predecessor_of_consolidated_evidence_audit'
        crosswalk.append([f'Archive_Table_{n:02d}', n, t['title'], retained.get(n, 1 if n == 54 else ''),
                          state, prefix + '.csv', prefix + '.tex', t['line']])
    write_csv('metadata/table_crosswalk.csv', crosswalk)
    for t in selected:
        prefix = f'tables/current_si/SI_Table_{t["number"]:02d}'
        write_csv(prefix + '.csv', t['rows'])
        write(prefix + '.tex', '% Verbatim current SI table fragment; manuscript text not included.\n' + t['fragment'] + '\n')
    numeric_matches = {}
    for old, new in retained.items():
        a = next(t for t in archived if t['number'] == old)['rows'][1:]
        b = next(t for t in selected if t['number'] == new)['rows'][1:]
        assert a == b, (old, new)
        numeric_matches[str(new)] = {'archive_number': old, 'rows_exactly_equal_after_tex_normalization': True, 'rows': len(b)}

    figure_rows = [['current_document', 'current_label', 'source_filename', 'sha256', 'matches_preceding_figure', 'primary_data']]
    figure_data = {
        'Figure_1_geological_support.pdf': 'geological_classes.nc;geological_background_cells.csv;age_sediment_area_histogram.csv;hydrothermal_catalogue_locations.csv;representativeness/',
        'Figure_2_paired_budget.pdf': 'SO240_stations_work_areas.csv;reorganization/F2b_violin_statistics.csv;main_SI_display/F2c_basin_area_ratios.csv;main_SI_display/F2c_allocation_power_intersections.csv',
        'Figure_3_energy_compatibility.pdf': 'retained/F2_observation_summary.csv;retained/F2_forcing_summary.csv;native_A_latitude_integral.csv;panels/F3c_primary_budgets.csv;main_SI_display/F3d_model_dependence_points.csv',
        'Figure_4_release_controls.pdf': 'physical_schematic_conditions.csv;figure4_display/F4b_plotted_values.csv;figure4_display/F4c_plotted_values.csv',
        'Figure_5_power_persistence_clear_labels.pdf': 'figure5_display/F5a_plotted_frontiers.csv;figure5_display/F5a_state_and_gap_display.csv;figure5_display/F5b_fixed_gain_curves.csv;figure5_display/F5b_cap_intersections.csv;figure5_display/F5b_saved_target_values.csv',
        'Supplementary_Figure_S2_allocation_sensitivity.pdf': 'retained/geological_support/conditional_budget/geology_conditioned_source_fields.nc;retained/integrated_display/selected_regional_diagnostic.csv',
        'Supplementary_Figure_S4_complete_history.pdf': 'retained/integrated_display/selected_global_annual_history.csv',
        'Supplementary_Figure_S5_time_depth_comparisons.pdf': 'independent_2025_O_A_N_comparisons.csv;time_depth_support.csv',
        'Supplementary_Figure_S6_regional_paired_budgets_clear_labels.pdf': 'panels/F2c_paired_heat_flow.csv',
        'Supplementary_Figure_S8_full_conditions.pdf': 'retained/evidence_closure/constraint_ablation.csv',
        'Supplementary_Figure_S9_condition_identities.pdf': 'panels/F3d_condition_memberships.csv;panels/F3_selected_memberships.csv',
        'Extended_Data_Figure_1_geological_strata.pdf': 'retained/global_release_support/stratum_support_envelopes.csv;retained/global_release_support/nested_regime_weights.csv',
    }
    for doc in ['manuscript', 'SI', 'Extended_Data']:
        text = (current / (doc + '.tex')).read_text()
        figures = re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text)
        for i, rel in enumerate(figures, 1):
            f = current / rel
            prior = base / 'figures' / f.name
            assert f.is_file() and sha(f) == sha(prior)
            label = ('Figure ' + re.search(r'Figure_(\d+)', f.name).group(1) if doc == 'manuscript'
                     else 'Supplementary Figure S' + str(i) if doc == 'SI' else 'Extended Data Figure ' + str(i))
            figure_rows.append([doc, label, f.name, sha(f), True, figure_data[f.name]])
    write_csv('metadata/figure_crosswalk.csv', figure_rows)
    write_json('metadata/manuscript_snapshot.json', {
        'authoritative_directory': CURRENT, 'preceding_source_data': BASE + '/Source_Data',
        'documents': {n: sha(current / n) for n in ['manuscript.tex', 'SI.tex', 'Extended_Data.tex', 'manuscript.pdf', 'SI.pdf', 'Extended_Data.pdf']},
        'current_si_table_count': 4, 'current_si_figure_count': 6, 'archived_table_count': 57,
        'retained_table_checks': numeric_matches, 'figure_hash_matches': len(figure_rows) - 1,
        'full_manuscript_and_figure_artwork_included': False,
        'repository_created': False, 'repository_url': None})

    scripts = {
        'geological_support_20261105T131309Z': ['certify_geological_budget.py', 'power_transfer_helpers.py', 'original_constraint_builder.py', 'integrate_geological_support.py', 'derive_evidence_diagnostics.py', 'verify_geological_support.py'],
        'physical_limits_20261106T004812Z': ['analyze_physical_limits.py', 'derive_exposure_frontier.py', 'verify_physical_limits.py', 'verify_exposure_frontier.py'],
        'NG_evidence_closure_20261106T033806Z': ['analyze_evidence_dependencies.py'],
        'global_release_support_20261106T040811Z': ['analyze_release_support.py', 'derive_sampling_design.py', 'verify_release_support.py'],
        'NG_external_validation_20261106T045511Z': ['analyze_portable_station_table.py', 'verify_external_validation.py'],
        'NG_regional_evidence_scope_20261108T095538Z': ['audit_existing_evidence.py'],
        'NG_main_SI_reorganization_20261108T144244Z': ['plot_publication.py', 'shared_plotting.py', 'colour_panels_publication.py', 'ppt_full.py'],
        'NG_integrated_editable_figures_20261109T031639Z': ['render_figure5_clear_labels.py', 'complete_figure_cleanup.py'],
    }
    for batch, names in scripts.items():
        for name in names:
            preserve(root / 'docs/manuscript' / batch / 'code' / name, f'code/original/{batch}/{name}')
    for name in ['NG_evidence_closure_20261106T033806Z.json', 'global_release_support_20261106T040811Z.json', 'NG_external_validation_20261106T045511Z_station_support.json']:
        preserve(root / 'model_configs' / name, 'model_configs/' + name)
    preserve(root / 'workflow/scripts/validate_ng_supporting_repository_20261116T074445Z.py', 'code/validate_repository.py')
    preserve(root / 'workflow/scripts/reproduce_ng_supporting_analyses_20261116T074445Z.py', 'code/reproduce_analyses.py')
    preserve(Path(__file__).resolve(), 'code/build_repository_from_project.py')

    write('README.md', '''# Ocean heat energy constraints — supporting data and code

Supporting materials for *Energy constraints on anthropogenic and seafloor contributions to modern ocean warming*.

This is a local repository candidate aligned with the author-selected manuscript and shortened Supplementary Information supplied on 16 November 2026. It is not yet a published repository or a DOI archive. `metadata/manuscript_snapshot.json` identifies the exact document hashes. No full manuscript, research plan, raw upstream archive, model installation, credentials, Windows metadata or compiled TeX log is included.

## Start here

- **Current SI**: four machine-readable tables in `tables/current_si/`; corresponding verbatim LaTeX table fragments retain formatting and units.
- **Supporting tables moved out of the SI**: all 57 preceding tables in `tables/archived_si/`. These use `Archive_Table_01`–`Archive_Table_57`, not the current SI numbering. `metadata/table_crosswalk.csv` distinguishes current matches, supporting detail and historical conditions. Older claims and definitions are not silently promoted to current results.
- **Figure and table source data**: `Source_Data/`, preserved byte-for-byte from the preceding integrated delivery. Final artwork hashes match for all five main figures, six current SI figures and one Extended Data figure.
- **Current figure numbering**: `metadata/figure_crosswalk.csv`. Current S1–S6 use older artwork filenames containing S2, S4, S5, S6, S8 and S9. The filename number is not the current caption number.
- **Reproduction**: `code/reproduce_analyses.py` regenerates analytical certificates, physical bounds, exposure frontiers, release-support calculations and the regional station comparison from supplied processed inputs. It does not run an ocean or climate model or download data.
- **Integrity and tabular checks**: `code/validate_repository.py` validates file hashes, current tables, matched budgets and source-data relationships.
- **Provenance**: `metadata/preserved_files.json` identifies source paths and SHA-256 hashes. `metadata/files.sha256` covers the repository content, excluding itself and the Git database.

## Reproduction commands

Use an existing Python environment with NumPy, pandas, SciPy and xarray. The original figure programs additionally require Matplotlib and Pillow. No installation is performed by these commands.

```bash
python3 -I code/validate_repository.py --repository .
python3 -I code/reproduce_analyses.py --repository . \\
  --output-dir /absolute/path/to/unused/analysis_output \\
  --log-dir /absolute/path/to/unused/analysis_logs
```

Both output directories must not exist. Original data and scripts are never overwritten. The second command uses frozen processed inputs; regeneration from original native geological grids or complete raw climate archives is a separate upstream stage and is not claimed by this command. The exact scope and remaining limitations are in `REPRODUCIBILITY.md`.

## Scientific interpretation

The primary source condition is geological-power allocation with S2_P3. Eligible-area allocation, 8/10 TW conditions and the expanded 14.4 MW condition are separate alternatives. Ordinary supply, additional anomaly power, gross axis-plus-flank power, source-only bounds, monthly timing relaxations and joint compatibility are distinct quantities. Fixed-mean benchmark percentages are not interchangeable with same-member paired ratios. Empty support and missing durations are not zero heat. Regional source allocations do not identify observed regional warming fractions. These materials do not establish actual global attribution or a universal measured vent-power cap.

## Numbering and historical metadata

Current SI Table 2 preserves preceding Table 31; Table 3 preserves Table 53; Table 4 preserves Table 50. Current Table 1 is the consolidated evidence audit, whose preceding counterpart is Table 54. Old table and figure references inside inherited files are provenance and must be resolved using the current crosswalk. In particular, inherited display manifests may describe earlier styling. The author-selected figure hashes and current caption numbering take precedence.

## Upstream data and rights

See `UPSTREAM_DATA.md`. Public availability does not by itself establish redistribution permission. Upstream datasets retain their own terms, attribution and citation requirements. No blanket license is assigned to third-party material. An author-approved license for original code/data and confirmation of redistribution terms remain release checks. See `PUBLICATION_CHECKLIST.md` before public upload.
''')
    write('REPRODUCIBILITY.md', '''# Reproduction scope

## Runnable without upstream downloads

The analytical runner uses original scripts without editing their scientific formulas. It creates an input layout with symlinks to supplied processed data, passes explicit paths, runs Python in isolated mode and stores logs outside the repository. It regenerates geological conditional certificates, physical monthly bounds, exposure frontiers, evidence-dependency diagnostics, heterogeneous release support and the SO240 regional comparison. All regenerated numerical CSV tables are compared to the frozen equivalents with declared floating-point tolerances. The original independent physical-kernel, scalar-dual exposure and heterogeneous-support verifiers are also executed.

The repository validator independently checks the current SI table row counts, current support counts and rounded ranges, nine primary allocation memberships, five budget identities, all three O=A+N+B+E windows and regional cooling-reference identities. A verification pass is numerical and accounting evidence, not empirical validation of a globally representative heat source.

## Original source code and upstream stages

`code/original/` preserves scientific and figure programs with their source batch identifiers. Some retained programs are archival entry points with original project-path requirements. Inclusion is not a claim that every archival program runs from the repository root. Native Seton/GlobSed integration and evidence-record grid lookups need the cited original grids, which are not redistributed here. The original climate/observation candidate generation also requires upstream archives. Frozen processed candidates, certificates, kernels and source configurations are supplied for the current analytical reproduction.

## Figures

`metadata/figure_crosswalk.csv` identifies all current plotted datasets. Final figure PDFs are SHA-256 identical to the identified preceding delivery. The repository does not include manuscript artwork or the supplied PowerPoint schematic. Original figure programs are retained for inspection, but full pixel-identical regeneration of every final figure is **not** claimed. Figure 4a is conceptual and its editable artwork requires a separate author-approved input. Figure programs also contain inherited path/layout assumptions. Numeric plotted coordinates are supplied, including the current Figure 4 and Figure 5 display exports. No schematic is treated as measured connectivity or temperature.

## Table exports

CSV tables preserve displayed values rather than manufacturing extra precision; full-precision calculations use `Source_Data/`. TeX commands for wrapping, emphasis and escaped punctuation are normalized in CSV exports. Mathematical notation not safely converted is retained as TeX. Exact table fragments accompany each CSV. Table fragments do not include all surrounding methodological prose; use the current manuscript and source-data definitions for interpretation.
''')
    write('UPSTREAM_DATA.md', '''# Upstream data references and redistribution review

These identifiers are transcribed from the selected manuscript. This packaging step has not independently fetched or relicensed upstream archives.

| Dataset | Cited identifier | Role |
|---|---|---|
| LargeEnsemble OHC reconstructions | https://doi.org/10.5281/zenodo.21376891 | Observational ensemble |
| IAPv4 | https://doi.org/10.5194/essd-16-3517-2024 | Gridded ocean context |
| DAMIP archive | https://doi.org/10.5281/zenodo.10121522 | Climate-response inputs |
| Seton 2020 crustal age | https://doi.org/10.1029/2020GC009214 | Geological screening |
| GlobSed v3 | https://doi.org/10.1029/2018GC008115; https://doi.org/10.1594/PANGAEA.982339 | Sediment screening; distinguish article and retained archive |
| Davies heat flow | https://doi.org/10.1002/ggge.20271 | Background diagnostic |
| InterRidge v3.4 | https://doi.org/10.1594/PANGAEA.917894 | Field-occurrence catalogue, not measured power |
| SO240 heat flow | https://doi.org/10.1594/PANGAEA.965467 | Regional station measurements |
| PB2002 | https://doi.org/10.1029/2001GC000252 | Plate-boundary context |
| ETOPO1 | NOAA public product, as cited in the manuscript | Geographic background |
| Local vent-power compilation | https://doi.org/10.1002/2015JB012245 | Separate local measurements and calibration groups |

The repository contains derived products and reused tabular subsets, including station/catalogue records. Their original terms and required attribution must be checked before public release; an open DOI is not proof of an open license. No raw native age/sediment grids, full climate archives, publisher PDFs or supplied figure artwork are included. The inherited data inventory files document additional lineage and evidence limitations. The SO240 erratum remains bibliographically identified but its correction text is not verified by this release preparation.
''')
    write('PUBLICATION_CHECKLIST.md', '''# Public release checklist

- [x] Use the author-selected November 16 manuscript and shortened SI as the authority.
- [x] Keep existing project files unchanged; do not stage the entire research workspace.
- [x] Export four current SI tables and preserve all 57 preceding tables under archival IDs.
- [x] Supply current figure-number crosswalk and source-file hashes.
- [ ] Review `metadata/validation.json` and `metadata/reproduction_validation.json` once generated; do not claim tests ran before those records exist.
- [ ] Confirm GitHub owner/organization and repository name. Suggested name: `ocean-heat-energy-constraints`.
- [ ] Authenticate through an approved GitHub client. The preparation environment initially had no `gh` command. Do not store tokens in this project.
- [ ] Confirm upstream redistribution terms for reused records and choose an author-approved original-code/data license. No license is assumed here.
- [ ] Create the requested PUBLIC repository and upload only this curated repository, after reviewing its exact file manifest.
- [ ] Verify the remote commit and public access. Create a versioned release; an archival DOI is recommended but not fabricated.
- [ ] Insert the verified repository URL in manuscript Data availability and replace `https://github.com/xxx` in Code availability. Add the repository link to the SI reading guide if approved. Preserve the existing manuscripts by creating dated linked versions unless in-place edits are explicitly approved.
- [ ] Compile linked document versions into a fresh output location, check links, references and layout, and retain logs under `results/logs/`.

No GitHub URL has been invented. The current manuscript's statement that code is already public is not verified until the upload and access check are complete. The public repository need not include the complete unpublished manuscript.
''')
    write('.gitignore', '__pycache__/\n*.py[cod]\n.DS_Store\n*:Zone.Identifier*\n*.aux\n*.log\n*.out\n')
    write_json('metadata/preserved_files.json', manifest)
    # A file/column inventory is a structural dictionary, not an inferred scientific glossary.
    inventory = [['path', 'bytes', 'sha256', 'rows_excluding_header', 'columns_json']]
    for f in sorted((out / 'Source_Data').rglob('*')):
        if not f.is_file():
            continue
        count, columns = '', []
        if f.suffix == '.csv':
            with f.open(newline='') as stream:
                reader = csv.reader(stream)
                columns = next(reader, [])
                count = sum(1 for _ in reader)
        inventory.append([f.relative_to(out).as_posix(), f.stat().st_size, sha(f), count, json.dumps(columns, ensure_ascii=False)])
    write_csv('metadata/source_data_inventory.csv', inventory)
    write('metadata/files.sha256', ''.join(f'{sha(f)}  {f.relative_to(out).as_posix()}\n' for f in sorted(out.rglob('*')) if f.is_file()))
    print(json.dumps({'repository': str(out), 'source_data_files': len(inventory) - 1,
                      'preserved_files': len(manifest), 'archived_tables': len(archived),
                      'current_tables': len(selected), 'figure_hash_matches': len(figure_rows) - 1}, indent=2))


if __name__ == '__main__':
    main()
