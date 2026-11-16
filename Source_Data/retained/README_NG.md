# Selected NG data

This document supplements the inherited input documentation. The authoritative figure-number mapping is ../figure_manifest_NG.csv.

- `main_window/stored_history_requirements.csv`: 2,102 budget histories and two separate source-only witnesses; complete-history inventory, power/flow and date-specific stocks. Only budget records have observation denominators.
- `main_window/epoch_power_audit.csv`: three modern epochs for each of those 2,104 records.
- `main_window/selected_reservoir_trajectories.csv`: actual/ordinary old stocks and young stocks at year end, plus annual-mean ocean B states, for identified witnesses.
- `main_window/saved_support_requirement_envelopes.csv`: direct-retention finite-support threshold summaries; each minimum keeps its own witness and companion coordinates. No support is NA, not zero demand.
- `necessary_conditions/necessary_conditions.csv`: analytical S2/direct necessary contact and flow floors for explicit target heat. The first three use percentages of mean O, not same-member ratios.
- `necessary_conditions/necessary_condition_witness_checks.csv`: comparisons between the local-reference analytical floors and saved nonnegative S2/direct records, including the source-only upper witness.

The analytical derivation preserves correlated references, full young reservation and normal-counterpart feasibility. It relaxes some constraints and therefore supplies necessary, not sufficient, conditions. It does not run the archived optimizer or create a new source ensemble. The preliminary parent-directory audit with first-window ratios is excluded.

Other files are unchanged inherited inputs for global, spatial and approximation-sensitivity results. Keep A, B, Bh and positive retained Ch distinct. The 1970–2014 spatial diagnostic is not a 1955–2023 observed partition. Geographic layers and external reuse conditions remain documented in ../map_sources.csv. All interpretation and reproduction details are in the selected manuscripts and ../REPRODUCTION_SELECTED.md.
