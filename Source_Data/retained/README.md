# Spatial Source Data

## Selected supports
Global attribution remains1955–2023,0–6000m and the original908-member product. Spatial A and source comparisons use1970–2014,0–6000m;1955model subsets did not complete their bounded transfer and are not used. No post2014spatial pattern extension is implemented. ACCESS-CM2 r1 historical/hist-nat is the sole available spatial subset, not the original model ensemble mean and not a climate model retained in the compatible S1/S2 support.

The selected annual operator averages12monthly column states equally, matching the frozen source ocean_states annual mean. A preliminary duration-weighted reduction was corrected only for this operator alignment. No source or global climate result was changed. The original model-global contrast remains a distinct support in integral_comparison.csv, not a target for normalization.

## Arrays and units
- native_anthropogenic_fields.nc: original300×360curvilinear A field, latitude/longitude in degrees, native area_m2, volume_0_6000_m3 and basin_index. A_ZJ is per native cell; A_GJ_m2 is column energy divided by wet-cell area. Native wet thickness is derived from reference cell volume/area and capped at6000m; rho1035kg/m3 and cp4000J/kg/K define the original heat diagnostic convention.
- mapped_fields.nc: conservative2degree native-cell-centroid sums of energy, area and volume, basin-area fractions, mixed-resolution native-A-plus-basin-B, and the sole heat-capacity-mixedB sensitivity. This preserves nativeA texture rather than inventing it.
- mapped_basin_fields.nc: selected main Fig1fields. A and B are both six-basin means, with fractional basin areas retained in boundary bins. Each *_ZJ is bin energy; *_GJ_m2 is energy/actual summed wet area. main A+B is not a matched original compatible tuple.
- saved_history_basin_support.csv: all1207direct histories ×4windows ×6basins, including original candidate/member/source IDs; B, Bh, mixed and conductive-donor components; region area and heat capacity. Negative responses are retained. Candidate source variables remain in frozen_direct_candidates.csv.
- saved_history_basin_components_ZJ.csv: seven explicit terms in ZJ. Young routing is split by actual axis/mixed discharge. Paired recipient/donor terms use the same conduction-region weights, cancelling in every region. Native class contributions use their own source-ledger weights.
- regional_support_endpoints.csv: regionwise extrema, each with its actual attaining candidate. These describe saved-history support, not exact regional extrema or maps that can be stitched into a global upper bound.
- selected_spatial_histories.csv: original global1955–2023B-energy-low/high candidates for each set; map integrals use1970–2014from the same source histories.
- regional_diagnostic.csv: independent nativeA, selectedsourceB/Bh, sum and alternate mixing by basin. O_ZJ and O_minus_A_minus_B_ZJ are deliberately missing, not zeros.
- native_endpoint_monthly_global_diagnostics.csv:48monthly absolute reference heat diagnostics and explicit annual weights; only endpoint differences define warming.
- basin_mask_assignment_audit.csv: every native wet-cell location/area, category and distance to nearest valid retained NOAA mask node. Polar and Mediterranean overrides are explicit. Unassigned marginal waters, notably the Red Sea assigned to Indian, can be far from the provider’s basin-integration mask; this is a geographical convention, not observed support.

B=Bh+mixed+conduction_donor. Bh is background-relative axis plus flank response including routing received by those classes; Ch is positive retained input and is not substituted. No23.03%global full-B endpoint is relabelled as a Bh upper.

The capacity-mixed alternative uses the same globalB and native rho*cp*volume, giving uniform water-column temperature increment but variable GJ/m2. Neither allocation is an observed transport response or a physical spatial envelope. MainB retains the original near-full-depth global comparison condition, including the assumption that omitted deeper response does not materially alter the scale.

## Actual missing information
No usable0–6000mspatial observation endpoint was obtained. The additional code/data capsule from the same Zenodo record contains918configuration global time series, not the selected908members' gridded endpoint histories; its map files are decadal spread or trend-uncertainty summaries, not observed warming fields. None of these capsule values replaces the selected908global ensemble. The IAP official metadata pages and cited author portal did not yield numeric fields within the bounded access pass. No O or residual heatmap, error envelope, or member matching is manufactured.

Original input paths/hashes are in spatial_input_inventory.csv; acquisition records and access failures accompany provenance/. Fields remain editable CSV/NetCDF rather than an opaque object or executable archive. All figures redraw offline from these files.
