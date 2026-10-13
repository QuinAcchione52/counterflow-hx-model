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
The given inlet remains an initial guess; the input CSV is not rewritten.
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
below 1e-6 W. See `validation/May27_10_HEOS.csv` for the summary. Ten automated
tests cover the laminar limiting solution, experimental HTC definition, inlet
wall correlation, counterflow boundaries/energy conservation, fitting behavior,
inlet inference, and independent hydraulic calibration. A real bounded fit to 1000 W/(m² K) at 40 segments also
returned e/D ≈ 0.00883879 and reproduced the target.

## Separate pressure-drop fit

```bash
python pressure_calibration.py --roughness 0.3636871236009118 --segments 100
```

This holds the thermal effective e/D fixed, infers the inlet to match the
reference mean bulk temperature, and then independently fits viscous pressure
loss to `pressure_drop_psi` in `experimental_targets.csv`. The signed
acceleration contribution is not multiplied. It outputs `Pressure-fit.xlsx`
and `Pressure-fit.summary.csv`; `--target-pressure-psi` can override a single
case's target. The calibration is a one-way postprocessing fit on the unchanged
thermal profile, not a new coupled solution or an experimentally established
friction law. Absolute-pressure feedback is omitted as in the original solver;
416.44 Pa is only about 0.00515% of the 8.094 MPa operating pressure here.

For May27_10 on the HTC-calibrated 100-segment profile:

| Quantity | Value |
|---|---:|
| Unscaled viscous drop | 0.00272173 psi |
| Acceleration contribution (cooling pressure recovery) | -0.0000342463 psi |
| Baseline total pressure drop | 0.00268748 psi |
| Experimental total pressure drop | 0.0604 psi |
| Fitted viscous multiplier | 22.20435 |
| Calibrated total pressure drop | 0.0604 psi |
| Unchanged HTC | 2642.176 W/(m² K) |

The script also independently inverts Colebrook for a **pressure-only**
equivalent e/D on the frozen thermal profile. Matching 0.0604 psi requires
**e/D ≈ 2.26324**, corresponding to about **30.30 mm** of equivalent roughness
in a 13.39 mm bore. This is geometrically impossible as wall roughness. The
extended 0..3 inversion is only a diagnostic to expose that contradiction;
it is never fed back into the Dang–Hihara HTC correlation. The multiplier
calibration and pressure-only e/D inversion are alternative ways of fitting
the same pressure target; do not apply both simultaneously.

A separate coupled thermal run with e/D=1 and the mean bulk matched predicts
only 0.0085014 psi and HTC about 4897 W/(m² K). Thus one common roughness
parameter cannot reproduce the reported HTC and pressure for this model and
these supplied conditions. The pressure fit has separate hydraulic and thermal
parameters and cannot be called a validated model simply because both targets
can be forced to match.

The thesis says the differential transducer measures the test-section drop,
but the actual tap-to-tap wetted path and all restrictions need confirmation.
Check reducers, connection tubes, fittings, valves, and zero/calibration records
before interpreting the missing loss. A Darcy/Fanning convention difference
(factor four) cannot explain the remaining factor 22 on this HTC-fitted profile.
The thesis also reports less than 4% pressure difference between smooth and
dimpled tubes, so assigning the entire measured loss to dimple roughness is
not supported by that observation. See `validation/May27_10_pressure_HEOS.csv`.
