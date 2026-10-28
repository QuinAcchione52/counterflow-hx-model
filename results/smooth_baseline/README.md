# Optional empirical smooth baseline: HEOS screening results

The correction is separate from published Dang–Hihara and from roughness. It
multiplies local Nu inside the coupled wall iteration; hydraulic friction is
not multiplied. Defaults remain uncorrected. The mean CO2 temperature, flow and
pressure come from Spencer; the CO2 inlet is inferred and shared water
conditions are assumptions. Nominal flux group is a regression predictor and
is not an imposed or independently matched heat flux.

- `training_factors.csv`: 38 individual fitted smooth factors; these are not
  predictions. All use e/D=0 and the original smooth DH friction expression.
- `heldout_predictions.xlsx` / `.csv`: each factor is predicted by regression
  trained on the other 37 cases, then used in a forward coupled solve. The
  held-out target HTC is excluded from that regression. Mean absolute HTC error
  decreases from 46.73% to 14.22%; 9/38 are within 5%, worst absolute error 43.46%.
- `model.json`: all-38 fitted regression for optional transfer to dimpled cases.
  Fixed geometry/water assumptions and property backend are enforced. Range
  extrapolation must be explicit and is labeled. Min/max ranges do not certify
  interpolation quality or account for feature combinations.
- `Dimpled-transfer.xlsx` / `.summary.csv`: smooth-trained correction held fixed;
  individual e/D fits bounded to 0..0.05. Unbracketed cases are not forced to
  match and no negative or enlarged roughness is used.

All-case screening uses 20 segments. Three representative full-model smooth
predictions checked at 100 segments changed HTC by less than 0.002%. This checks
numerical resolution in those cases, not validity of the correlation. Closely
related cases from a single campaign mean leave-one-case-out results may be
optimistic for a new campaign. No experimental uncertainty was modified.
Neither fitted HTC agreement nor a small fitted e/D proves physical roughness,
correct heat flux, pressure-drop agreement, or the source of the baseline error.

The completed transfer screen fits **24/40 dimpled cases** within e/D=0..0.05.
The other 16 are unbracketed: 13 targets fall below the corrected e/D=0
prediction, and three exceed the e/D=0.05 prediction. This count describes HTC
calibration reachability; 0.05 is a tested bound, not an established physical
roughness limit for dimples. Dimpled HTC targets are used to fit individual e/D,
so this is not held-out validation of dimple predictions.

Examples at 20 segments: May19_3 e/D=0.010524 and May27_10 e/D=0.015321.
Retain the spreadsheet heat-flux columns: matching HTC does not imply matching
the nominal 2/3 kW/m² experimental group or the actual measured per-case duty.

May24_1 is flagged as outside a recorded smooth-training range. The transfer
batch explicitly allowed this baseline extrapolation; it is not silently
included as an in-range prediction. For example, the matched May19_3 HTC still
has predicted flux about 1436 W/m² versus the nominal 2000 W/m² group; May27_10
predicts about 1412 W/m² versus its nominal 3000 group. These comparisons are
warnings about incomplete operating-condition reproduction, not comparisons to
exact measured per-case heat flux (which is unavailable in the summary tables).

`Dimpled_fine_grid_checks.csv` confirms the two example fits at 100 segments:
May19_3 e/D=0.01052448, HTC=1862.00; May27_10 e/D=0.01531948,
HTC=2642.00 W/(m² K). These remain conditional empirical fits.

## Complete 40-case workbook

`All_40_Dimpled_Corrected_Results.xlsx` collects all 40 targets, both bounded
endpoint HTC predictions and all 40 smooth-trained multipliers in one main
sheet. It also includes the full solver summary, source experimental data,
assumed inputs, correction metadata, smooth held-out validation, two fine-grid
checks and the 24 available matched axial profiles (32 sheets total).

For 16 unfitted cases, final fitted e/D, duty, outlet, fitted HTC and pressure
outputs remain blank. The workbook does not invent a selected state for an
unbracketed case. Dimpled data use the final correction trained on all smooth
cases; they are not held-out dimple predictions. Rebuild with
`python analysis/export_dimpled_corrected_workbook.py` (pandas, xlsxwriter and
openpyxl required; no solver rerun).
