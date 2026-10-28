"""Export all saved dimpled results with explicit fit/validation labels."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/smooth_baseline/All_40_Dimpled_Corrected_Results.xlsx'
base=ROOT/'results/smooth_baseline'
r=pd.read_csv(base/'Dimpled-transfer.summary.csv')
t=pd.read_csv(ROOT/'results/dimpled_trends/All_40_old_vs_new_eD.csv').set_index('Case').loc[r.Case].reset_index()
i=pd.read_csv(ROOT/'Dimpled_cases.csv').set_index('Case').loc[r.Case].reset_index()
e=pd.read_csv(ROOT/'experimental_targets.csv').set_index('Case').loc[r.Case].reset_index()
model=json.loads((base/'model.json').read_text())
s=r.copy();s['CO2_baseline_correction_multiplier']=t.smooth_baseline_multiplier.to_numpy()
s['CO2_baseline_is_calibrated']=True
s['CO2_heat_friction_model']=s.CO2_heat_friction_model.fillna('dh_roughness')
s['property_backend']=s.property_backend.fillna(model['backend'])
status=t.new_status.map({'matched_0_to_005':'HTC fitted within e/D 0–0.05','target_below_smooth':'Unfitted: target below corrected smooth HTC','target_above_005_cap':'Unfitted: target above e/D=0.05 HTC'})
main=pd.DataFrame({
 'Case':r.Case,'Source test group':i.test_group,'Fit status':status,
 'Spencer target HTC (W/m² K)':r.target_htc_W_m2K,
 'Numerical HTC at e/D=0 (W/m² K)':r.smooth_model_htc_W_m2K,
 'Numerical HTC at e/D=0.05 (W/m² K)':r.htc_at_eD_005_W_m2K,
 'Numerical fitted HTC (W/m² K)':r.htc_spencer_W_m2K,
 'Fitted effective e/D':r.effective_e_over_D,
 'Smooth-trained Nu multiplier':t.smooth_baseline_multiplier,
 'Baseline range extrapolated':t.baseline_extrapolated,
 'Fitted HTC error (%)':100*r.relative_htc_error,
 'Reported experimental HTC uncertainty (%)':r.experimental_htc_uncertainty_percent,
 'Spencer mean CO₂ temperature (°C)':r.target_mean_CO2_C,
 'Numerical mean CO₂ temperature (°C)':r.CO2_mean_endpoint_C,
 'Inferred CO₂ inlet (°C)':r.inferred_CO2_inlet_C,
 'Numerical CO₂ outlet (°C)':r.CO2_out_C,
 'Experimental CO₂ mass flow (kg/s)':r.experimental_CO2_flow_kg_s,
 'Reported CO₂ mass flux (kg/m² s)':t.mass_flux_kg_m2_s,
 'Computed CO₂ mass flux (kg/m² s)':t.computed_mass_flux_kg_m2_s,
 'Experimental inlet pressure (MPa)':r.experimental_pressure_MPa,
 'Assumed water inlet (°C)':r.assumed_water_inlet_C,
 'Assumed water flow (kg/s)':r.assumed_water_flow_kg_s,
 'Numerical water outlet (°C)':r.water_out_C,
 'Numerical CO₂ heat duty (W)':r.duty_CO2_W,
 'Numerical CO₂ heat flux (W/m²)':r.heat_flux_CO2_W_m2,
 'Nominal experimental heat-flux group (W/m²)':r.nominal_heat_flux_W_m2,
 'Flux difference from nominal (%)':100*(r.predicted_to_nominal_flux_ratio-1),
 'Numerical CO₂ pressure drop (psi)':r.CO2_pressure_drop_psi,
 'Experimental CO₂ pressure drop (psi)':e.pressure_drop_psi,
 'Numerical pressure-drop error (%)':100*(r.CO2_pressure_drop_psi/e.pressure_drop_psi-1),
 'Numerical length-mean local HTC (W/m² K)':r.htc_local_length_mean_W_m2K,
 'Bulk-to-wall LMTD (K)':r.bulk_wall_lmtd_K,
 'CO₂ energy residual (W)':r.CO2_energy_residual_W,
 'Water energy residual (W)':r.water_energy_residual_W,
 'Segments':r.segments,'Property backend':s.property_backend})
assert len(main)==40 and main.Case.is_unique
assert main['Smooth-trained Nu multiplier'].notna().all()
assert main['Fitted effective e/D'].notna().sum()==24
notes=[
 ['Purpose','All 40 dimpled cases after applying the final smooth-trained empirical Nu correction; saved results, not a new solver run.'],
 ['Important validation distinction','The multiplier relationship was tested by leaving out each smooth case (9/38 within 5%), then trained on ALL 38 smooth cases for dimpled transfer. Dimpled outputs are NOT held-out dimple predictions.'],
 ['How dimpled matching was performed','The smooth-trained multiplier is calculated from operating conditions and held fixed. Only e/D is fitted to each dimpled HTC target, bounded to 0..0.05.'],
 ['Cases available','40 targets and endpoint HTC predictions; 24 matched profiles. Thirteen targets fall below corrected smooth HTC; three exceed the HTC at e/D=.05.'],
 ['Blank output cells','An unfitted case has no selected final e/D or fitted profile. Duty, outlet, fitted HTC and pressure cells are deliberately blank, not zero. Both endpoint HTC predictions remain available for every case.'],
 ['All 40 multipliers','Multiplier values for unbracketed cases are taken from the previously exported operating-condition correction, rather than leaving them blank as in the raw batch summary.'],
 ['Experimental vs assumed','CO₂ pressure, flow, mean temperature and target HTC come from Spencer. Water inlet 20.7°C and flow .087 kg/s are shared assumptions; CO₂ inlet is inferred. Original CSV inlet values are guesses, not measured inlets.'],
 ['Nominal heat-flux group','2/3 kW/m² is the nominal test group, not exact measured per-case duty. The multiplier uses this group as a predictor; duty is still solved by the HX.'],
 ['Fit is not physical validation','Matching HTC does not validate physical roughness, pressure drop or heat flux. The .05 limit is a chosen bound, not a validated dimple-roughness limit.'],
 ['Extrapolation','May24_1 is flagged outside the smooth calibration ranges; the transfer batch explicitly allowed baseline range extrapolation.'],
 ['Uncertainty','Reported experimental HTC uncertainty is unchanged. Prediction error is a separate metric.'],
 ['Grid and properties','Main batch uses 20 segments and HEOS. Two separately saved 100-segment fits are on Fine_grid_checks; do not substitute those e/D values into the 20-segment batch without labeling.'],
 ['Profile tabs','The 24 tabs beginning with case names preserve the original matched axial profiles. No profiles are invented for the 16 unfitted cases.'],
 ['Sources','Model-results snapshot da28165; Dimpled-transfer.xlsx/.summary.csv; Dimpled_cases.csv; experimental_targets.csv; All_40_old_vs_new_eD.csv; model.json.'],
]
metadata=[]
for key,value in model.items():metadata.append([key,json.dumps(value,ensure_ascii=False) if isinstance(value,(dict,list)) else value])
source=pd.ExcelFile(base/'Dimpled-transfer.xlsx')
frames={'Read_me':pd.DataFrame(notes,columns=['Topic','Explanation']),
 'All_40_cases':main,'All_solver_results':s,'Experimental_data':e,'Inputs_and_assumptions':i,
 'Correction_model':pd.DataFrame(metadata,columns=['Setting','Value']),
 'Smooth_heldout_validation':pd.read_csv(base/'heldout_predictions.csv'),
 'Fine_grid_checks':pd.read_csv(base/'Dimpled_fine_grid_checks.csv')}
for name in source.sheet_names:
 if name not in ['Summary','Inputs_assumptions']:frames[name]=pd.read_excel(source,sheet_name=name)
with pd.ExcelWriter(OUT,engine='xlsxwriter') as writer:
 wb=writer.book
 header=wb.add_format({'bold':True,'font_color':'white','bg_color':'#142B45','text_wrap':True,'valign':'top','border':0})
 number=wb.add_format({'num_format':'0.000000'});wrap=wb.add_format({'text_wrap':True,'valign':'top'})
 for name,frame in frames.items():
  frame.to_excel(writer,sheet_name=name,index=False);ws=writer.sheets[name];ws.freeze_panes(1,1);ws.autofilter(0,0,len(frame),len(frame.columns)-1)
  ws.set_row(0,46,header)
  for j,column in enumerate(frame.columns):
   ws.write(0,j,column,header)
   ws.set_column(j,j, max(17,min(27,len(str(column))*.8)),number if pd.api.types.is_numeric_dtype(frame[column]) else None)
  if name in ['Read_me','Correction_model']:
   ws.set_column(0,0,32,wrap);ws.set_column(1,1,112,wrap)
   for row in range(1,len(frame)+1):ws.set_row(row,62 if name=='Read_me' else 52)
  if name=='All_40_cases':
   ws.set_column(0,0,16);ws.set_column(2,2,46,wrap)
   ws.conditional_format(1,2,40,2,{'type':'text','criteria':'containing','value':'Unfitted','format':wb.add_format({'bg_color':'#FCE8D5','text_wrap':True})})
   ws.conditional_format(1,2,40,2,{'type':'text','criteria':'begins with','value':'HTC fitted','format':wb.add_format({'bg_color':'#DBF0E9','text_wrap':True})})
   ws.set_landscape();ws.fit_to_pages(1,0);ws.repeat_rows(0)
print(f'Saved {OUT} ({len(frames)} sheets; all 40 cases, 24 axial profiles).')
