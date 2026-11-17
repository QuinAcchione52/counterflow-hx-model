# MATLAB plots for the 40 dimpled cases

Download `plot_dimpled_results.m` and run it in MATLAB. This is the supplied
script with Figures 7 and 8 appended; the original figures are retained.

Keep these two workbooks in Downloads (or select them when prompted):

- `All_40_Dimpled_Corrected_Results.xlsx`
- `Spencer-vs-Numerical-HTC-All-40-Cases.xlsx`

The additional figures reproduce the supplied temperature PDF: four source
classes, reported bulk temperature on x, fitted e/D on y, red original fits,
blue available revised fits, and green dashed pressure-dependent pseudocritical
reference lines. Figure 7 has all 40 original fits; Figure 8 shows only the 24
revised fits, with y limits 0..0.01. Nineteen values are above this display
range; five are visible. Off-scale counts are labeled per panel. Unfitted revised values remain absent, not zero.

The four HEOS Cp-maximum references are embedded so CoolProp is not required in
MATLAB. Each is calculated at the mean pressure of all cases in its class;
these are group references, not the fixed critical temperature or the exact
pseudocritical temperature of every individual run. A pressure check refuses
silent reuse with a different dataset.

PDF and PNG versions are exported to `Dimpled_Temperature_Plots` alongside the
selected corrected workbook. Uses `readtable` with preserved headers,
`xline`, `sgtitle`, and `exportgraphics` (MATLAB R2020b or newer recommended).

Workbook columns, case joins, class counts and mean pressures were checked
against the saved workbooks. MATLAB/Octave is unavailable in the execution
workspace, so MATLAB execution and rendering were not tested here.
