# First correlation audit: equations, roughness profiles, and solver bookkeeping

The equation check uses the user's uploaded Spencer thesis, Eqs.4–10 (printed
pp.9–10). This is a check against that transcription, not independent verification
of the original Dang–Hihara article. Direct primary-reference lookup was blocked
by the execution environment's network proxy; no verification was disabled.

## Findings and changes

1. **Step with one step silently ignored e/D.** `linspace(0,ep,1)` returned zero.
   A one-step profile now means the requested constant e/D. Multiple steps use
   physical normalized x rather than node-index partitions. Changing the grid
   therefore does not move a step boundary.
2. **Custom e/D arrays were silently reduced to their first value.** Arrays now
   require one finite nonnegative e/D per node and are used without truncation.
   Invalid step counts, negative roughness, and malformed arrays are rejected.
3. **The previous zero-roughness calculation was not the stated DH smooth
   friction equation.** `dang_hihara` uses `(1.82*log10(Re_f)-1.64)^(-2)` with
   e/D=0. The new default `dh_roughness` multiplies this smooth friction by
   `Colebrook(Re_f,e/D)/Colebrook(Re_f,0)`. It recovers the smooth equation
   exactly at zero and introduces a local relative roughness enhancement.
   This ratio extension is empirical, not a validated published rough-tube DH
   correlation. `colebrook` retains the previous direct replacement for comparison.
4. **Invalid property states were clipped without changing the stated inputs.**
   Correlations now reject invalid states instead of silently replacing bulk or
   wall temperatures. Water pressure is always Pa at the correlation boundary;
   it is no longer guessed to be atm from its magnitude. The outer-wall state
   is evaluated at its actual calculated temperature.
5. **Water wall iteration was fixed at three iterations.** It now tests
   temperature convergence. The fully developed laminar water baseline and
   turbulent/transition model have not otherwise been replaced or calibrated.
6. **Thermal and hydraulic friction were conflated.** Thermal DH friction uses
   film Reynolds number as transcribed; mechanical pressure friction now uses
   bulk Reynolds number with Colebrook (or laminar/transition handling). Segment
   viscous pressure integrates the endpoint `f/rho` values. These two friction
   factors are exported separately. This is not a pressure calibration.
7. **Fluid phase labels used incorrect numeric enum assignments.** They now use
   CoolProp's named phase constants; actual bulk phase is included in the output.
8. **Correlation branch selection was hidden in the final HTC.** Output now
   includes film Reynolds number, bulk/integrated heat capacity, selected
   Prandtl number and branch, film conductivity, and the local wall-flux residual.

9. **Energy closure did not guarantee local wall-flux closure.** With a stepped
   roughness profile, averaging segment duty onto nodes could disagree with
   the node conductance and wall temperatures. Coupled node flux now uses the
   local conductance and bulk temperature difference; segment duty remains
   integrated separately for the energy balance. Final CO2/water wall-flux and
   wall-conduction residuals are exported and water-wall closure is required
   for convergence. A three-step test checks both sides and conduction.
10. **An arbitrary initial 3 kW/m² duty could create a false freezing state.**
    Coupled runs now start at zero duty and build the flux from the solution.
    Prescribed-flux mode retains its requested flux. Actual temperatures remain
    unclipped; this changes only the numerical initial guess.

The CO2 heat-transfer inputs are explicitly K, Pa, m, kg/(m² s), and W/m² into
CO2. Positive q (heating) is rejected because the selected model is for cooling.
Bulk Re<3000 is rejected as a conservative turbulent-regime guard; this is not
an assertion that all larger Reynolds numbers are within DH's validation data.
Geometry and low mass-flux screening flags are based on Spencer's cited
1–6 mm tube database and approximately 200 kg/(m² s) experiments, not a claim
of full published applicability bounds.

## What the checks establish

- Equal wall/bulk, zero-flux DH is checked with independent CoolProp properties
  against the transcribed equations.
- The roughness extension has the exact original smooth limit.
- Colebrook is checked against the standard Re=100000 smooth Darcy factor
  (about 0.0179898) and its implicit equation residual.
- Physical profile boundaries survive grid changes; one step equals constant;
  custom arrays preserve every node value; sine remains nonnegative.
- Counterflow boundary conditions, the joint bulk-temperature inverse solve,
  and both energy balances remain tested.
- Constant, Step, Linear, and Sine profiles were exercised in the actual HX.
  For the example with one step, Constant and Step return the same solution.

At 40 segments under the assumed water conditions, May27_10 at e/D=0 gives
740.46 W/(m² K) with the stated DH smooth expression versus 731.98 with the
previous Colebrook replacement. At e/D=0.05, the normalized DH extension gives
1304.15 versus 1291.30. May19_3 smooth gives 705.64 versus 698.53. These changes
are about 1%, not an explanation of the earlier large experimental discrepancy.

The historical all-case tables were produced using the previous direct
Colebrook thermal friction and earlier pressure bookkeeping. They remain
archived results; their exact fitted e/D values and pressure values should not
be presented as newly generated results from the revised formulation. The
`colebrook` option preserves the earlier thermal formulation, but corrected
hydraulic friction means pressure results need not be identical.

## Remaining physical limitations

DH was developed for smooth, smaller tubes and substantially higher CO2 mass
fluxes than these Spencer cases. Roughness enhancement alone is not a model of
buoyancy, dimple-induced secondary flows, or altered wetted geometry. The water
model remains a fully developed forced-convection baseline; thermal entrance
and mixed-convection effects have not been independently validated.

For nonuniform profiles, the experimental endpoint-LMTD HTC and the
length-mean local HTC are different. Profiles that return to smooth at both
ends (e.g. Sine) can enhance interior heat transfer while changing the endpoint
metric much less. Compare duty, local HTC, wall/bulk profiles, and both average
metrics; do not force a single endpoint metric to represent every roughness
profile's performance.
