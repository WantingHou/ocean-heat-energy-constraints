"""Add complete table context and verified reproduction records without overwriting."""
from pathlib import Path
import argparse
import csv
import importlib.metadata
import json
import re
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-root', type=Path, required=True)
    ap.add_argument('--repository', type=Path, required=True)
    ap.add_argument('--reproduction-dir', type=Path, required=True)
    a = ap.parse_args()
    root, repo, reproduced = a.project_root.resolve(), a.repository.resolve(), a.reproduction_dir.resolve()
    source = root / 'docs/manuscript/NG_integrated_editable_figures_20261109T031639Z/Supplementary_Information_integrated_EN.tex'
    document = source.read_text()
    outputs = ['metadata/validation.json', 'metadata/reproduction_validation.json', 'metadata/table_context.json',
               'tables/Supporting_Tables.tex', 'requirements-tested.txt', 'RELEASE_STATUS.md']
    for rel in outputs:
        if (repo / rel).exists():
            raise FileExistsError(repo / rel)

    def write(rel, content):
        with (repo / rel).open('x', encoding='utf-8') as f:
            f.write(content)

    subprocess.run([sys.executable, '-I', '-B', str(repo / 'code/validate_repository.py'), '--repository', str(repo),
                    '--output', str(repo / 'metadata/validation.json')], check=True)
    report = json.loads((reproduced / 'reproduction_validation.json').read_text())
    if report['passed'] is not True:
        raise ValueError('Reproduction did not pass')
    # Public-facing metadata keeps reproducibility results but omits local absolute log locations.
    for run in report['executed_scripts']:
        run['log'] = 'external_reproduction_logs/' + Path(run['log']).name
    write('metadata/reproduction_validation.json', json.dumps(report, indent=2) + '\n')
    context = []
    parts = []
    sections = list(re.finditer(r'\\section\*\{([^\n]+)\}', document))
    for i, match in enumerate(sections):
        table = re.fullmatch(r'Supplementary Table (\d+)\. (.+)', match.group(1))
        if not table:
            continue
        number = int(table.group(1))
        end = sections[i + 1].start() if i + 1 < len(sections) else document.index(r'\end{document}')
        body = document[match.end():end].strip()
        assert body.count(r'\begin{longtable}') == body.count(r'\end{longtable}') == 1
        before, rest = body.split(r'\begin{longtable}', 1)
        after = rest.split(r'\end{longtable}', 1)[1]
        context.append({'archive_table': f'Archive_Table_{number:02d}', 'previous_title': table.group(1),
                        'source_line': document.count('\n', 0, match.start()) + 1,
                        'before_table_tex': before, 'after_table_tex': after})
        parts.append((number, r'\section*{Repository Table A' + f'{number:02d}. ' + table.group(2) + '}\n' + body + '\n'))
    assert len(parts) == 57
    preamble = r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern,textcomp,amsmath,geometry,graphicx,longtable,booktabs,hyperref,fancyhdr,lineno,array}
\geometry{left=14mm,right=14mm,top=14mm,bottom=14mm}
\hypersetup{hidelinks,pdftitle={Ocean heat supporting tables}}
\urlstyle{same}
\pagestyle{fancy}\fancyhf{}\fancyfoot[C]{\thepage}\renewcommand{\headrulewidth}{0pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{6pt}\renewcommand{\baselinestretch}{1.13}
\setlength{\emergencystretch}{3em}\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.20}
\begin{document}
{\Large\bfseries Repository supporting tables\par}
Energy constraints on anthropogenic and seafloor contributions to modern ocean warming\par
\section*{Scope and numbering}
This companion preserves all 57 numerical and evidence tables from the preceding integrated Supplementary Information, including the immediately surrounding table notes and the original supplementary reference list. A01--A57 are repository identifiers, not the four table numbers of the current shortened Supplementary Information. No values have been recomputed for this document.\par
Current SI Tables 2, 3 and 4 correspond to A31, A53 and A50. Current SI Table 1 is the consolidated evidence audit with predecessor A54. Current CSV tables and the full numbering crosswalk are supplied separately. References to other table or figure numbers inside the preserved text follow the preceding SI numbering. The current manuscript and the repository figure crosswalk govern present figure labels.\par
Earlier-condition comparisons are preserved for traceability. They do not replace current source conditions or identify actual global attribution. Tables A07--A09, A11, A15 and A17--A26 belong to historical conditions; A10 is an auxiliary extended-window comparison. Distinguish additional from gross power, fixed-mean from paired percentages and separate necessary minima from joint attainable histories. Missing or empty results do not mean zero heat.\par
'''
    refs = next((i, m) for i, m in enumerate(sections) if m.group(1) == 'Supplementary References')
    i, m = refs
    reference_block = document[m.start():sections[i + 1].start()]
    write('tables/Supporting_Tables.tex', preamble + '\n'.join(text for _, text in sorted(parts)) + '\n' + reference_block + '\n\\end{document}\n')
    write('metadata/table_context.json', json.dumps(sorted(context, key=lambda x: x['archive_table']), indent=2, ensure_ascii=False) + '\n')
    packages = ['numpy', 'pandas', 'scipy', 'xarray', 'matplotlib', 'Pillow']
    write('requirements-tested.txt', '# Versions present during validation; no installation was performed.\n' + '\n'.join(f'{name}=={importlib.metadata.version(name)}' for name in packages) + '\n')
    write('RELEASE_STATUS.md', '''# Release preparation status

**Local preparation verified; no public GitHub repository has been created or uploaded.**

## Completed checks

- 229 preserved Source Data files; no original project data were overwritten.
- 57 preceding SI tables exported as CSV and verbatim table fragments; four current SI tables exported separately.
- `tables/Supporting_Tables.tex` supplies complete compilable table bodies, surrounding notes and the preceding SI reference list. Per-table `.tex` files are raw archival fragments, not standalone compilable documents. A PDF is included only after successful compilation.
- Current SI Tables 2/3/4 match preceding Tables 31/53/50; current Table 1 is the consolidated evidence audit.
- Twelve current figure hashes match the identified source delivery. Full artwork is not redistributed here.
- Current SI membership/range checks, three-window budget closure, regional cooling identities and SHA-256 checks pass. See `metadata/validation.json`.
- Ten original analytical/verifier programs ran successfully. Thirty-four frozen CSV/NetCDF products matched at the declared tolerances. See `metadata/reproduction_validation.json`.
- The independently checked exposure frontier has 4,986 points and 72 inversion judgements. These checks establish numerical consistency, not empirical global attribution.

## Publication still pending

The owner/organization, authentication and remote repository URL have not been supplied. Suggested repository name is `ocean-heat-energy-constraints`; the author selected a public repository. The original manuscript still contains its placeholder URL. No speculative URL or DOI has been inserted.

Author-approved licensing and attribution review remain pending, especially for reused third-party records. No license has been invented. Original manuscript, SI and Research Plan files remain unchanged. An actual verified GitHub URL is required before creating linked manuscript versions.

`metadata/files.sha256` describes the initial immutable package. Additional verification records and this document are added afterward. `metadata/release_files.sha256`, when present, is the complete release-candidate manifest excluding itself. Inherited source configurations, README files and display manifests retain their original historical wording and path assumptions; the current numbering crosswalk and this release status take precedence for packaging.
''')
    print(json.dumps({'added': outputs, 'full_table_contexts': len(context), 'reproduction_passed': report['passed']}, indent=2))


if __name__ == '__main__':
    main()
