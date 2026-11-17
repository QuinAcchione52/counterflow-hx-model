# Dimpled CO2 pressure-drop audit: all 40 cases

Download `Pressure_drop_audit_all_40.xlsx` for the comparison, inputs, method
notes and all 40 thermal-state profiles. `Pressure_comparison_all_40.pdf` / PNG
shows experimental-vs-predicted pressure drop and the ratio vs bulk temperature.
`All_40_pressure_audit.csv` is the full machine-readable comparison.

## Findings

- All 40 experimental targets lie outside the predictions from hydraulic
  e/D=0..0.05 on the selected fixed thermal states. No pressure fit is forced.
- At hydraulic e/D=.05, the experimental/model ratio is **69.54–218.09**;
  mean absolute percentage error is **99.12%**.
- Experimental table values correspond to **341–692 Pa**. Modeled e/D=.05
  total losses are **2.23–8.34 Pa**.
- Reintegrating saved profiles independently reproduces all 24 existing
  pressure outputs to within **1.9e-13 Pa**. Darcy friction, signed acceleration
  and psi conversion do not show an arithmetic discrepancy in this check.
- Cooling acceleration is negative here (a pressure recovery), only about
  **-0.370 to -0.015 Pa**. Its sign is appropriate and it cannot explain the gap.
- Assigning the entire measured drop to this 1m tube would require apparent
  Darcy factors **4.85–16.43**, versus modeled mean f **.0728–.0758** at e/D=.05.
  These are diagnostics, not measured tube friction factors or validated fits.

For May27_10, measured table ΔP is .0604 psi, versus .000773 psi at hydraulic
e/D=.05. The saved HTC-fitted roughness gives .000487 psi. For May19_3, the
corresponding figures are .1004, .000589 and .000330 psi. The .05 sensitivity
is therefore a different comparison from the previously quoted pressure errors
at each case's HTC-fitted e/D.

## Method and limits

The existing thermal solver and HTC calibration are **unchanged**. Twenty-four
saved HTC-fitted profiles are reused. The 16 unmatched cases are solved at
thermal e/D=0, with their existing smooth-trained multiplier and reported mean
temperature. These are labeled smooth-baseline states, not fitted dimple HTC
solutions. Each hydraulic roughness sensitivity uses the same thermal states.
Pressure properties use HEOS at the input pressure, just as the current HX.

The independent calculation uses:

- Darcy friction: integral `f_D G²/(2 rho D) dx`, with trapezoidal endpoint f/rho.
- Signed acceleration: `G² (1/rho_out - 1/rho_in)`.
- `1 psi = 6894.757293168 Pa`; the archived solver uses rounded 6894.757.

No pressure multiplier, dimple form-drag correlation, extra K or sensor offset
is introduced. Apparent f, equivalent additional K, equivalent extra pipe length
and equivalent CO2 static head are **alternative scale diagnostics only**.
For example the residual corresponds to about **5.6–21.2 cm of CO2 static head**
at inlet density. That does not demonstrate an actual tap elevation error:
impulse-line fluids, routing and zeroing must be known. Added equivalent lengths
of 66–216m are not inferred actual tap spans. No adjustment is made from these
columns.

## Source verification required before calibration

The original thesis/technical-note PDFs are not available in the current
workspace. Pressure targets and psi units here are from the previously
extracted Table 13 data (`experimental_targets.csv`), **not freshly verified
against the paper**. Pressure-tap span, restrictions, elevation, impulse-line
routing, transducer range/units, zero offset and pressure uncertainty remain
unverified. HTC uncertainty must not be used as pressure uncertainty.

Recheck the source pressure units and tube/tap geometry first; distinguish
straight-tube friction, dimple form drag, minor losses and measurement effects.
A Darcy/Fanning conversion alone is a factor of four and cannot account for
this 70–218 ratio. The measured ΔP is at most .00865% of inlet pressure; ignoring
pressure feedback is unlikely to explain this large discrepancy, although a
fully pressure-coupled property solve has not been performed.

Only after the measurement boundary is established should any missing-loss
model be calibrated, with held-out tests and parameters separate from the HTC
multiplier and thermal e/D. A match obtained by fitting an arbitrary loss
coefficient would not establish the physical source of the discrepancy.

Reproduce with the Python solver dependencies plus matplotlib/openpyxl:

```bash
python analysis/audit_dimpled_pressure.py
python -m pytest tests/test_pressure_audit.py -q
```

The script checks saved-profile agreement and the test module checks an
independent uniform-pipe Darcy reference, cooling acceleration sign and grid
independence, and invalid-state rejection. No existing HTC results are changed.
