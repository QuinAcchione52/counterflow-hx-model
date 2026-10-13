# counterflow-hx-model

Counterflow CO₂/water heat exchanger with a Dang–Hihara CO₂ correlation and an
annular water-side model. Calibration compares the **experimental-definition
average HTC**, not the arithmetic mean of the local correlation values.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python CounterHX_full.py --match-bulk --reference experimental_targets.csv
python -m pytest -q
```

The default property backend is CoolProp **HEOS**. To use your licensed REFPROP
installation, set `HX_PROPERTY_BACKEND=REFPROP` and configure CoolProp's REFPROP
path. Backend errors are not silently replaced with another backend. Results
record the backend; cloud results have not been validated against REFPROP.

`Test.csv` provides the inlet conditions; `D (mm)` means **CO₂ inner diameter**.
The thesis specifies 0.625 inch OD and 0.049 inch wall: ID **13.3858 mm**, wall
**1.2446 mm**, OD **15.875 mm**. The technical note rounds the wall to 1.2 mm;
this input uses the thesis's inch dimensions. Water gap is a radial gap.
The user confirmed that the CO₂ inlet is a guess chosen to reproduce the paper's mean bulk temperature. `target_mean_CO2_C=35.23` and `target_htc_W_m2K=2642`
are supplied for May27_10. When `target_mean_CO2_C` is present, the runner infers the inlet by default
instead of fixing the guess. Use `--fixed-inlet` only when you want to hold the
supplied inlet fixed. `--match-bulk` explicitly enables inference from the
reference table for other inputs. The water inlet is held at the supplied 20.7 °C; its experimental
provenance remains unverified.

Outputs contain case profiles and a `Summary` sheet, plus a companion summary
CSV. Pressure-drop columns are Pa except the explicitly labeled psi column.
Heat transfer `Q` is negative for CO₂ cooling; `Q_leak` is positive for heat
leaving water to ambient. The linear leak coefficients represent whole-section
heat **into water**, using inlet water temperature minus ambient temperature,
and are distributed once over the segments. Confirm that sign convention with
the original experimental leakage fit.

## Adjust or fit e/D

```bash
python CounterHX_full.py --roughness 0.01 --output sweep-001.xlsx
python CounterHX_full.py --fit --match-bulk --reference experimental_targets.csv --ep-max 0.05
# Alternatively, a target for a single-case input:
python CounterHX_full.py --fit --match-bulk --bulk-temp 35.23 --target-htc 2642
```

For multiple cases, supply their measured inlet conditions in the input CSV and
add `target_htc_W_m2K` per case, or join `experimental_targets.csv` with
`--reference`. That reference contains the 40 **dimpled** runs from thesis
Table 13 (printed p. 56); it does not contain all inlet conditions needed to run
the solver. Each Case must match exactly. `ep_max (e/Dh)` controls the profile;
use `Constant` to calibrate a single uniform effective e/D per case.

The fit solves for experimental-definition HTC:

`h = CO₂ duty / (π × inner diameter × length × LMTD_bulk_to_inner_wall)`

Endpoint wall temperatures are calculated at both ends. A separate
length-weighted local HTC mean is also reported. The fit requires a sign-changing
bracket and a converged solution, then verifies HTC within 0.1%. It rejects
unbracketed targets before opening the output file. Default bounds are 0..0.05;
Bounds above 0.05 require explicit `--allow-extrapolation` when fitting.
The numerical guard is 1.0. Values beyond the usual rough-pipe range are only
empirical fitting knobs, **not physical roughness or a validated Colebrook
extension to dimpled tubes**. To investigate a numerical match explicitly:

```bash
python CounterHX_full.py --fit --match-bulk --ep-max 1 --allow-extrapolation
```

`--match-bulk` solves for CO₂ inlet temperature separately at every trial e/D,
so the endpoint mean stays at the measured bulk temperature throughout fitting.
Outputs label that inlet `inferred_CO2_inlet_C` and `CO2_inlet_is_inferred`.
The CSV inlet is used as a fixed boundary with `--fixed-inlet`. In inverse
mode the solver uses a numerical initial guess based on the target mean; the
input CSV is not rewritten.
The tool does not add an independent HTC enhancement multiplier.

## May27_10 diagnosis

Using HEOS, 100 segments, the corrected geometry, supplied inlet conditions,
and e/D=0.001 gives approximately:

| Quantity | Model | Thesis Table 13 |
|---|---:|---:|
| Experimental-definition HTC, W/(m² K) | 849 | 2642 |
| Length-mean local HTC, W/(m² K) | 868 | 3170 (experimental sampled mean) |
| Mean endpoint CO₂ temperature, °C | 36.12 | 35.23 |
| CO₂ pressure drop, psi | 0.000361 | 0.0604 |
| CO₂ heat flux, W/m² | 1391 | nominal 3000 test group |

The old code reproduces roughly 543 W/(m² K) local HTC with the outer diameter
entered as ID and its forced-turbulent water branch. Correcting geometry and
endpoint treatment does **not** remove the full experimental discrepancy.
A 20/40/100-segment check changes the corrected average HTC by under 0.02%.
Both fluid energy residuals are below 1e-5 W.

An e/D sweep at 20 segments gives about 821 (smooth), 849 (0.001), 1374 (0.05),
and 1641 (0.1) W/(m² K). **2642 is not bracketed** even by 0..0.1 for this model
and input. Do not interpret a large effective roughness as measured roughness,
and do not tune HTC while ignoring duty, mean temperature, or pressure drop.

The water inlet Reynolds number is about 1269. The previous code clamped it to
2300 and applied a turbulent correlation. The new baseline integrates the
fully developed laminar annular velocity/energy equations for a uniform inner
wall heat flux and insulated outer wall; 2300..4000 uses linear transition to
Gnielinski. This removes the artificial turbulent floor, but **fully developed,
forced laminar flow is not a validated description of the experimental water
annulus**. Thermal entrance effects, buoyancy/mixed convection, and bath inlet
mixing need investigation. The laminar model is a stated baseline, not a fitted
water correlation. CO₂ friction uses film Reynolds number to preserve the
Dang–Hihara convention, with Colebrook replacing its smooth friction formula.
The papers explicitly identify extrapolation at these low CO₂ mass fluxes.

A diagnostic combination of the guessed 36.5 °C inlet plus
the thesis's 35.23 °C endpoint mean implies an outlet of 33.96 °C. HEOS predicts
about 245 W / 5.83 kW/m² for those endpoints, versus the nominal 3 kW/m² group.
That is not an experimental duty estimate, because the inlet is guessed.
The group label is not a substitute for the measured duty. Without raw endpoint
and wall temperatures, the reference's case label alone cannot establish that
all CSV inputs describe that experiment.

`Constant=True` in the solver is a diagnostic imposed-flux calculation; it does
not enforce the water correlation's conjugate wall balance and is not used by
the calibration CLI. `whole_section=True` is rejected because the old branch
contained placeholder wall temperatures and effectiveness.

Matching HTC alone is calibration, not validation of the exchanger model.
After the user clarified that the inlet is guessed, a bulk-matched sweep gives:

| Effective e/D | Inferred CO₂ inlet, °C | Spencer-definition HTC, W/(m² K) |
|---|---:|---:|
| 0 | 35.4413 | 732 |
| 0.001 | 35.4421 | 756 |
| 0.01 | 35.4467 | 917 |
| 0.05 | 35.4534 | 1291 |
| 0.1 | 35.4565 | 1600 |

All rows hold mean CO₂ temperature at 35.23 °C. The change in mean temperature
is no longer confounded with roughness. A match beyond this range must be
explicitly reported as extrapolated and checked against the other observables.


The bulk/HTC joint fit at 20 segments reaches **2642 W/(m² K)** with
**effective e/D ≈ 0.363687** and **inferred inlet ≈ 35.462078 °C**. The predicted
mean temperature is 35.230000 °C. This is a numerical fit in the extrapolated
Colebrook extension, not a physically credible tube roughness: it implies an
absolute roughness of about 4.87 mm in a 13.39 mm bore. Predicted heat flux is
about 1412 W/m² and pressure drop about 0.00269 psi, so other experimental
observables still do not match. The numerical local HTC mean is about 2665,
versus the experimental sampled mean of 3170. The two HTC averages should not
be interchanged. All numerical results in this diagnosis use HEOS.


Fine-grid confirmation at **100 segments** gives HTC **2642.176 W/(m² K)**
(error **0.00665%**) at the same fitted e/D, with mean bulk temperature
**35.2299995 °C** and inferred inlet **35.4620801 °C**. Energy residuals are
below 1e-6 W. See `validation/May27_10_HEOS.csv` for the summary. Nine automated
tests cover the laminar limiting solution, experimental HTC definition, inlet
wall correlation, counterflow boundaries/energy conservation, fitting behavior,
and inlet inference. A real bounded fit to 1000 W/(m² K) at 40 segments also
returned e/D ≈ 0.00883879 and reproduced the target.


For the supplied May27_10 CSV, changing `ep_max (e/Dh)` and running
`python CounterHX_full.py` now keeps the mean bulk temperature at 35.23 °C
by default. Alternatively use `--roughness 0.01` to override e/D without
rewriting the CSV. `--fixed-inlet` disables this inverse-temperature step for
cases where the inlet is a measured condition. Automatic fitting still requires
`--fit`; large extrapolated fitting bounds require `--allow-extrapolation`.

## Dimpled cases, one at a time or as a batch

`Dimpled_cases.csv` contains all 40 dimpled runs from thesis Table 13, with each
case's published CO2 mass flow, pressure, mean bulk temperature, and HTC target.
The mass flow and pressure retain the table's rounding. Water inlet **20.7 °C**,
water flow **0.087 kg/s**, ambient **22 °C**, geometry, and the existing leakage
fit are shared assumptions inherited from the supplied May27_10 input. They
are not per-case measurements from the thesis. The nominal heat-flux group is
recorded separately from model-predicted flux.

```bash
# First published dimpled case; fit only within the conventional test bounds:
python calibrate_dimpled.py --case May19_3 --output May19_3-fit.xlsx

# Explicitly allow an effective calibration beyond those bounds:
python calibrate_dimpled.py --case May19_3 --ep-max 1 --allow-extrapolation --output May19_3-fit.xlsx

# All 40, with independent cases in separate processes:
python calibrate_dimpled.py --workers 3 --ep-max 1 --allow-extrapolation
```

The runner checks the zero-roughness and e/D=0.05 predictions before trying an
extended bracket. A root that matches HTC within 0.1% is labeled either
`fit_within_0_to_005` or `extrapolated_fit`; the first label describes the tested
bounds, not proof of physical roughness. Unbracketed targets and solver errors
are reported separately and do not abort the other cases. The output contains
a Summary, Inputs_assumptions, and successfully fitted profiles; a partial CSV
is refreshed as cases finish. Default 20 segments is a screening grid; confirm
important cases at 100 segments before interpreting differences between cases.

The unknown CO2 inlet is now solved jointly with the water outlet inside each
frozen-U thermal solve. This avoids repeatedly converging the whole model for
trial inlet temperatures. It preserves the same energy balance and endpoint
mean constraint. It has been checked against the previous nested solver.

## Water-side sensitivity

Water conditions can affect the inferred CO2 inlet, heat duty, wall temperature,
and therefore the HTC correlation. They must not be tuned freely alongside e/D
and then interpreted as independently identified physical parameters.

For actual temperature/flow sensitivity, edit `T_H2O (C)` or `m_H2O (kg/s)` in
a copy of the input CSV. An optional `water_htc_multiplier` column or
`--water-htc-multiplier` flag is available for **diagnostic sensitivity only**;
its default is 1.0 and it does not establish a new validated water correlation.
It scales water HTC and adjusts its correlation wall-flux iteration consistently.
It never directly scales CO2 HTC. No water multiplier is fitted automatically.

```bash
python CounterHX_full.py --roughness 0.05 --water-htc-multiplier 2 --output Water-sensitivity.xlsx
```

At the May27_10 target mean of 35.23 °C and e/D=0.05, the tested water changes
have not raised the predicted CO2 HTC to 2642 W/(m² K). The baseline gives about
1291, a 2x water HTC gives 1207, and a 4x water HTC gives 1087. Increasing water
HTC increases duty but also changes the CO2 wall/bulk properties; it is not a
direct way to increase the CO2-side correlation coefficient. The default water
model and conditions therefore remain unchanged for the dimpled-case batch.
A physically supported entrance/mixed-convection water model remains a separate
model-development question; arbitrary multipliers are not a correction for it.


### Completed dimpled screening run

All 40 cases fitted their published HTC and bulk-temperature targets with the
shared assumptions above. Four point fits lie within 0..0.05: May24_4
(0.0197812), May27_4 (0.0172585), May24_2
(0.04078125), and May27_2
(0.02593258). The other 36 require extrapolated effective e/D.
The overall screening range is about 0.01726..0.64462; this variation is in an
effective calibration parameter for the same physical tube, not changing
measured surface roughness. See `validation/Dimpled_cases_HEOS_screening.csv`.

The first listed dimpled case, May19_3, was refitted at 100 segments:
e/D=0.28503290, HTC=1862.0000 W/(m² K), mean CO2 temperature=35.75 °C,
inferred inlet=36.269523 °C. The e/D differs by 0.0733% from the 20-segment fit.
May24_4 was also checked at 100 segments: e/D=0.01977386, HTC≈332 W/(m² K);
the e/D difference is -0.0372%. These are representative fine-grid checks,
not a claim that every case has been run at 100 segments.

The summary records the reported HTC uncertainty band separately from the
point-fit result. It is not treated as an independently established confidence
interval, nor are the point targets changed to obtain smaller roughness values.
Water-side sensitivity results are in `validation/May27_10_water_sensitivity_HEOS.csv`.

## Smooth-tube baseline comparison

`Smooth_cases.csv` contains the 38 smooth runs from thesis Table 12, with source
labels preserved. The comparison fixes e/D=0 and **does not fit HTC**. It uses
the same assumed water inlet of 20.7 °C, water flow of 0.087 kg/s, baseline water
correlation, and leakage assumptions. The CO2 inlet is inferred to match the
reported mean bulk temperature. This is a conditional base-model comparison,
not validation using independently measured inlet conditions.

```bash
python validate_smooth.py --workers 3 --segments 100
```

All 38 cases completed at 100 segments and satisfied energy/bulk constraints.
One case is within 5% HTC (11-Apr: 217.15 predicted vs 212 reported, +2.43%).
The mean absolute per-case percentage error is **46.96%**; 33 cases underpredict.
The signed mean error is -42.83%. None are within the reported per-case HTC
uncertainty bands under these shared assumptions. The 11-Apr case is borderline:
2.43% error versus a reported 2.42% band, with rounded source inputs.

A separate CO2-only check bypassed the water model and imposed the nominal
2 or 3 kW/m² heat-flux group, inferred inlet/outlet from CO2 enthalpy balance and
the reported mean, and evaluated the same zero-roughness CO2 correlation.
Its mean absolute HTC percentage error is **45.14%**. For smooth peak May11_2:
HTC is 622.42 vs 2268 W/(m² K), while the nominal-flux-implied inlet is 36.0186 °C,
close to the paper's reported 36.01 °C (printed p.37). This is a diagnostic using
nominal heat flux, not a substitute for exact per-case measured duty.

The mismatch therefore occurs even without dimples and is not removed by
bypassing the assumed water-side model in that diagnostic. It is not proof of
experimental error or a definitive validation of the code. Thesis section 7.2 /
Table 8 (printed pp.38–39) already reports 40–62% smooth-tube deviations for Liu,
which it says agrees better than Dang–Hihara, and identifies mass flux below the
correlations' validation ranges. These experiments also use a larger tube and
lower nominal flux than Dang–Hihara's cited database.

The thesis procedure (printed p.29) describes prechilling water for low-temperature
runs and adjusting water temperature to reach the flux target. Accordingly,
constant water temperature across every experimental case is an approximation.
See `results/Smooth-tube-HTC-comparison.xlsx` for both checks and profiles;
CSV values are in `validation/Smooth_zero_roughness_HEOS.csv` and
`validation/Smooth_CO2_only_nominal_flux_HEOS.csv`.

To reproduce the water-bypassed diagnostic:

```bash
python smooth_nominal_flux_check.py --segments 100
```
