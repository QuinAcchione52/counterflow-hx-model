# Pressure comparison: what Spencer's documents establish

Sources: *Summer 2026 Marinac_sCO2_Thesis_v10.pdf* and
*Enhanced_sCO2_Tubes_Technical_Note.pdf*, uploaded by the user. The repeated
uploads are byte-identical to the originals. This review concerns the pressure
comparison; it is not an independent audit of every thesis result.

## Measurement boundary

- Thesis printed pp. 18–19: the heated/cooled test tube is 1000 mm long,
  0.625 inch OD, 0.049 inch wall. This gives a 13.3858 mm ID.
- Thesis printed p. 18: 1000 mm entrance length and 500 mm exit length are
  installed around the test section. Their existence does not prove the pressure
  taps span all 2500 mm.
- Thesis printed p. 23: absolute pressure sensors are at the entrance and exit;
  a GP50 model 216 differential transducer measures test-section pressure drop.
- Technical note Fig. 1 (PDF p. 6): the dashed DPT sensing lines connect just
  upstream/downstream of the drawn test section. Pump, flow meter, coils, and
  relief valve are outside these drawn taps. The diagram gives neither actual
  tap spacing nor a dimensioned connection/reducer geometry.
- Technical note pressure reduction (PDF p. 10) and thesis printed p. 31:
  published pressure drop is measured directly by the differential sensor;
  subtracting absolute pressure readings is only a consistency check.

**Conclusion:** using the total loop pressure drop is not supported by the
schematic. Exact tap-to-tap hydraulic geometry remains unknown. Connection
losses cannot be excluded, but cannot be asserted to be the cause.

## Sensor calibration and uncertainty

Thesis printed p. 30, Eq. 13 / Table 5 gives the pressure conversion as
`pressure = m*V + b`; the DPT entry has slope 1 and offset +0.00509.
The reported pressure units in the surrounding text are psi. That offset is
about 8.4% of the May27_10 value; it cannot alone explain the factor-22 mismatch.
Published reduced pressure values should not have the listed offset subtracted
again without the original raw-signal reduction records.

Both instrumentation tables quote differential-pressure accuracy ±0.05%, but
these documents do not establish the exact installed sensor range or whether
this percentage refers to full scale or reading. Consequently an absolute
pressure uncertainty in psi cannot be calculated from that percentage alone.
Raw DPT voltages/readings, the installed range/output configuration, and
zero-pressure records at operating conditions are required for that check.

Possible unequal hydrostatic heads in the impulse lines are another item to
check, not a proven error. The fitted thermal inlet state gives density about
433 kg/m³; approximately 0.098 m of uncompensated CO₂ head could produce the
entire 0.0604 psi reading. The horizontal test-section description and schematic
do not establish the impulse-line elevations, fluid temperatures, or whether
that head has already been compensated.

## Experimental trend

Thesis printed pp. 47–48, section 7.6, and technical note section 4.4 report
less than 4% difference between smooth and dimpled pressure drops. The thesis
reports little observed temperature dependence over 30–40 °C and approximately
0.01 psi increase with the higher mass-flow setting. These observations are not
consistent with attributing the entire missing loss solely to dimples. They do
not establish a sensor error or identify a specific missing component.

May27_10, thesis Table 13, printed p. 56:
pressure 8094 kPa, rounded CO₂ flow 0.0043 kg/s, mean temperature 35.23 °C,
average HTC 2642 W/(m² K), local sampled mean HTC 3170 W/(m² K), pressure drop
0.0604 psi. The user confirmed the 36.5 °C CSV inlet is guessed, not measured;
model comparisons therefore infer the inlet to hold the published mean fixed.

## Separate fits

At the HTC-calibrated effective e/D=0.363687, 100-segment HEOS solution:

| Pressure term | psi |
|---|---:|
| Viscous loss | 0.00272172966 |
| Acceleration (cooling recovery) | -0.00003424634 |
| Total model loss | 0.00268748332 |
| Experiment | 0.0604 |

Matching the experiment by scaling **only viscous loss** requires a multiplier
22.204353. The thermal profile, HTC of 2642.176 W/(m² K), inferred inlet of
35.462080 °C, and signed acceleration term remain unchanged. This is a separate
empirical hydraulic calibration on a frozen thermal solution, not a prediction
of the reported pressure loss from known geometry.

Independently fitting a hydraulic Colebrook e/D on that same frozen thermal
solution gives e/D=2.263244, equivalent roughness 30.2953 mm. That exceeds the
13.3858 mm bore and cannot describe real wall roughness. This extrapolated
inversion is a diagnostic only and is not fed back into the HTC correlation.
The two pressure fits are alternatives and are not applied together.

Even a coupled thermal run using e/D=1, with the bulk mean matched, gives only
0.0085014 psi and HTC about 4897 W/(m² K). Thus increasing the common thermal
roughness does not preserve the HTC match while recovering pressure.

The missing loss corresponds to an additional equivalent K around 370 using
the modeled inlet dynamic pressure. This K is not a measured fitting coefficient
and is not evidence of a particular valve or reducer. A Darcy/Fanning factor
of four would also be insufficient to explain this difference.

## What would resolve the physical interpretation

A dimensioned or annotated view of both DPT process taps, the exact wetted
path between them (including any diameter changes), and raw DPT calibration /
zero records would let us test a physical loss model against the experiment.
Until then, the hydraulic multiplier can reproduce the number but should not
be interpreted as tube roughness or used as a validated cross-case correlation.
