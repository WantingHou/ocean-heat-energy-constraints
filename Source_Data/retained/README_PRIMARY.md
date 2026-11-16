# Field-calibrated primary data

Use `field_power/certified/` for the new S2_P3 global results:

- all_candidate_classification.csv: every original S2 complete discrete tuple, including rejections and required epoch anomalies.
- certified_histories.csv: 23 complete source certificates, with true paired observation/climate identities and unchanged required heat contrasts.
- certified_bounds.csv: all energy and same-member ratio endpoints over three windows and three operators, with attained witness IDs.
- source_only_bounds.csv: source-only extrema, distinct from joint-budget bounds.
- reference_power_envelopes.csv and proof.json: correlated-reference and temporal power-box necessity/sufficiency certificate.

Use `field_power/spatial/` for the selected source-dependent figures:

- primary_mapped_fields.nc: unchanged native-derived A, updated B and Bh, original wet area, and constructed-certificate spread; actual 1970–2014 window.
- primary_basin_histories.csv and primary_regional_diagnostic.csv: signed regional energies and complete-history identities.
- primary_annual_history.csv: selected GISS A, B and observation member 33321111, referenced to 1955; not annual two-source closure.
- selected_witness.json and verification.json: candidate 92589 and independent original-matrix checks.

`field_power/condition_package.json` records S1, S2 and S2_P3 and the field-evidence distinction. The other inherited data, including main_window/ and necessary_conditions/, refer to archived S1/S2 analyses, not S2_P3. They remain for explicit stress comparisons and supplementary figures. Do not replace their set identities or mix their source histories with the primary figures.
