# Ocean heat energy constraints — supporting data and code

Supporting materials for *Energy constraints on anthropogenic and seafloor contributions to modern ocean warming*.

This repository contains supporting materials aligned with the author-selected manuscript and shortened Supplementary Information supplied on 16 November 2026. Repository URL: https://github.com/WantingHou/ocean-heat-energy-constraints. No archival DOI has been assigned. `metadata/manuscript_snapshot.json` identifies the exact source-document hashes. No full manuscript, research plan, raw upstream archive, model installation, credentials, Windows metadata or compiled TeX log is included.

## Start here

- **Current SI**: four machine-readable tables in `tables/current_si/`; corresponding verbatim LaTeX table fragments retain formatting and units.
- **Supporting tables moved out of the SI**: all 57 preceding tables in `tables/archived_si/`. These use `Archive_Table_01`–`Archive_Table_57`, not the current SI numbering. `metadata/table_crosswalk.csv` distinguishes current matches, supporting detail and historical conditions. Older claims and definitions are not silently promoted to current results.
- **Figure and table source data**: `Source_Data/`, preserved byte-for-byte from the preceding integrated delivery. Final artwork hashes match for all five main figures, six current SI figures and one Extended Data figure.
- **Current figure numbering**: `metadata/figure_crosswalk.csv`. Current S1–S6 use older artwork filenames containing S2, S4, S5, S6, S8 and S9. The filename number is not the current caption number.
- **Reproduction**: `code/reproduce_analyses.py` regenerates analytical certificates, physical bounds, exposure frontiers, release-support calculations and the regional station comparison from supplied processed inputs. It does not run an ocean or climate model or download data.
- **Integrity and tabular checks**: `code/validate_repository.py` validates file hashes, current tables, matched budgets and source-data relationships.
- **Provenance**: `metadata/preserved_files.json` identifies original scientific-data/code paths and SHA-256 hashes. `metadata/files.sha256` covers all release content except the two checksum manifests and the Git database; `metadata/release_files.sha256` additionally covers `metadata/files.sha256`.
- **Readable supporting tables**: [Supporting_Tables.pdf](tables/Supporting_Tables.pdf) contains all 57 archived tables, their accompanying notes and references. See `RELEASE_STATUS.md` for verification scope.

## Reproduction commands

Use an existing Python environment with NumPy, pandas, SciPy and xarray. The original figure programs additionally require Matplotlib and Pillow. No installation is performed by these commands.

```bash
python3 -I code/validate_repository.py --repository .
python3 -I code/reproduce_analyses.py --repository . \
  --output-dir /absolute/path/to/unused/analysis_output \
  --log-dir /absolute/path/to/unused/analysis_logs
```

Both output directories must not exist. Original data and scripts are never overwritten. The second command uses frozen processed inputs; regeneration from original native geological grids or complete raw climate archives is a separate upstream stage and is not claimed by this command. The exact scope and remaining limitations are in `REPRODUCIBILITY.md`.

## Scientific interpretation

The primary source condition is geological-power allocation with S2_P3. Eligible-area allocation, 8/10 TW conditions and the expanded 14.4 MW condition are separate alternatives. Ordinary supply, additional anomaly power, gross axis-plus-flank power, source-only bounds, monthly timing relaxations and joint compatibility are distinct quantities. Fixed-mean benchmark percentages are not interchangeable with same-member paired ratios. Empty support and missing durations are not zero heat. Regional source allocations do not identify observed regional warming fractions. These materials do not establish actual global attribution or a universal measured vent-power cap.

## Numbering and historical metadata

Current SI Table 2 preserves preceding Table 31; Table 3 preserves Table 53; Table 4 preserves Table 50. Current Table 1 is the consolidated evidence audit, whose preceding counterpart is Table 54. Old table and figure references inside inherited files are provenance and must be resolved using the current crosswalk. In particular, inherited display manifests may describe earlier styling. The author-selected figure hashes and current caption numbering take precedence.

## Upstream data and rights

See `UPSTREAM_DATA.md`. Public availability does not by itself establish redistribution permission. Upstream datasets retain their own terms, attribution and citation requirements. No blanket license is assigned to third-party material. An author-approved license for original code/data and confirmation of redistribution terms remain release checks. See `PUBLICATION_CHECKLIST.md` before public upload.
