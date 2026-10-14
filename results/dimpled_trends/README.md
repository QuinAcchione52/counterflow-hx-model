# Dimpled effective roughness: historical vs revised fits

This comparison uses the model-results snapshot `da28165`. It analyzes saved
results and does not rerun the solver or alter experimental observations.

Download **Dimpled_eD_comparison_and_trends.xlsx** for a filterable 40-case
comparison, group summaries and exploratory Pearson/Spearman coefficients.
**Dimpled_eD_comparison_and_trends.pdf** contains reading notes and six plots.
Individual figures are also available as numbered PNGs.

## Which fits are being compared?

- Historical: direct-Colebrook thermal extension, no smooth-trained multiplier;
  e/D up to 1 was allowed. All 40 had fitted HTC, only four with e/D <= 0.05.
- New: revised DH-relative-roughness extension with the smooth-trained empirical
  correction; e/D restricted to 0..0.05. Twenty-four matched. Thirteen targets
  were below the corrected smooth baseline and three above the .05 endpoint.
- This includes code/formulation and baseline-calibration changes, not a
  controlled test of one change. Historical 'DH' means the earlier model using
  a Colebrook replacement, not the original smooth-only published correlation.
- All-case screening data are at 20 segments. Blank new e/D is **not zero**.
  Bounds symbols in figure 01 are flags, not fitted roughness values. No
  numerical fit beyond the bound or negative e/D is invented.

## Main observations

For the **same 24 paired cases**, median old e/D is **0.25257**, versus new
median **0.01726**. The median of the individual percentage reductions is
**93.83%**. This is distinct from comparing medians over different samples.

| Source group | Cases | Historical e/D <= .05 | New matched | Below smooth | Above cap |
|---|---:|---:|---:|---:|---:|
| m12/q2 | 11 | 0 | 7 | 4 | 0 |
| m12/q3 | 10 | 0 | 5 | 4 | 1 |
| m15/q2 | 9 | 2 | 7 | 1 | 1 |
| m15/q3 | 10 | 2 | 5 | 4 | 1 |

The source m12/m15 labels are used as group identifiers. Plots use reported
actual mass flux (about 24 and 30 kg/m² s), not the label numbers as mass flux.
The heat-flux predictor is the **nominal** 2000/3000 W/m² experimental group,
not exact measured duty. Numerical flux is preserved separately as an output.

In the new matched subset, Spearman associations with fitted e/D are:

- Nominal heat-flux group: **+0.476** (higher e/D tends to occur at 3 kW/m²).
- Mean CO2 temperature: **+0.400**; Pearson is only **+0.016**, so there is no
  useful simple linear temperature relationship across this subset.
- Reported mass flux: **-0.234**, a weaker association.

The historical all-40 temperature Spearman correlation is **+0.594** and
mass-flux correlation **-0.368**. The workbook additionally compares historical
fits for only the same 24 cases to distinguish sample-selection effects.

## Interpretation limits

These are **exploratory fitted-parameter associations, not validated predictive
correlations or physical roughness laws**. The 16 unavailable new fits are
systematically selected by model reachability. Correlations on the remaining
24 can be biased; no complete-40 new-e/D correlation is justified.

The smooth correction already depends on temperature relative to the
pseudocritical point, mass flux and nominal flux group. New e/D associations
therefore partly reflect that correction and remaining baseline/operating-
condition error. The groups have different temperature coverage. A fitted e/D
cannot be attributed solely to geometry, and the .05 cutoff is a tested bound,
not proof of physically credible dimples. HTC matching does not establish duty
or pressure-drop agreement. Reported experimental uncertainties are unchanged.

To reproduce with the solver environment (requires SciPy, CoolProp, pandas,
matplotlib and xlsxwriter):

```bash
python analysis/compare_dimpled_roughness.py
```

Sources: `results/Spencer-vs-Numerical-HTC-All-40-Cases.csv`,
`results/smooth_baseline/Dimpled-transfer.summary.csv`, `Dimpled_cases.csv`,
`experimental_targets.csv`, and `results/smooth_baseline/model.json`.

## Exceptions and case transitions

Of the original four cases within the .05 bound, only two remain matched:
May24_4 and May24_2. May27_4 and May27_2 now have targets below the corrected
smooth baseline. The new set therefore contains **22 newly reachable cases +
2 retained cases = 24**; it is not simply the original four plus 20.

Within the 24 paired fits, e/D decreases for 23 cases. The exception is
**May24_2**, increasing from 0.040781 to 0.049887. A correction that improves the
overall count does not improve every individual case.

For the matched subset, median e/D is **0.01314** in the nominal 2 kW/m² group
(n=14), versus **0.02368** in the 3 kW/m² group (n=10). This group comparison is
conditional on the same selection and correction-dependence limitations above.
