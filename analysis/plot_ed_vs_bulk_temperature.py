"""Replace run index with reported bulk temperature, separated by source class."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/hx-temperature-plot')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr,pearsonr
from scipy.optimize import minimize_scalar
import CoolProp.CoolProp as CP
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/dimpled_trends'
d=pd.read_csv(OUT/'All_40_old_vs_new_eD.csv')
assert len(d)==40 and d.Case.is_unique
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
@lru_cache(maxsize=256)
def tpc_C(pressure_MPa):
 state=CP.AbstractState('HEOS','CO2')
 if pressure_MPa*1e6<=state.p_critical():
  raise ValueError('Pseudocritical reference requires pressure above the critical pressure')
 def negative_cp(T):
  state.update(CP.PT_INPUTS,pressure_MPa*1e6,float(T))
  return -state.cpmass()
 optimum=minimize_scalar(negative_cp,bounds=(298.,320.),method='bounded',options={'xatol':1e-8})
 if not optimum.success:raise RuntimeError('Heat-capacity maximum search did not converge')
 return float(optimum.x-273.15)
references={}
for group,r in d.groupby('test_group',sort=True):
 pressure=float(r['P (Mpa)'].mean())
 references[group]={'source_group':group,'reference_pressure_mean_MPa':pressure,'pseudocritical_temperature_C':tpc_C(pressure),'pressure_min_MPa':float(r['P (Mpa)'].min()),'pressure_max_MPa':float(r['P (Mpa)'].max()),'case_Tpc_min_C':tpc_C(float(r['P (Mpa)'].min())),'case_Tpc_max_C':tpc_C(float(r['P (Mpa)'].max())),'backend':'HEOS','definition':'Cp maximum at group arithmetic-mean inlet pressure'}
figures=[];stats=[]
for compare in [False,True]:
 fig,axs=plt.subplots(2,2,figsize=(12.5,8),sharex=True,sharey=True)
 for ax,(group,r) in zip(axs.ravel(),d.groupby('test_group',sort=True)):
  if not compare:
   ax.scatter(r.bulk_CO2_C,r.old_fitted_eD,color='red',s=49,label=f'Original fits (n={len(r)})',zorder=3)
  if compare:
   p=r.dropna(subset=['new_fitted_eD'])
   ax.scatter(p.bulk_CO2_C,p.new_fitted_eD,color='blue',s=42,label=f'Corrected fits (n={len(p)})',zorder=4)
  if not compare:
   ax.axhline(.05,color='#666666',ls='--',lw=1,label='Reference e/D = 0.05')
  else:
   above=int((p.new_fitted_eD>.1).sum())
   if above:
    ax.text(.98,.05,f'{above} of {len(p)} blue values above 0.1 (off-scale)',transform=ax.transAxes,ha='right',va='bottom',fontsize=8,color='#444444')
  ref=references[group]
  ax.axvline(ref['pseudocritical_temperature_C'],color='#00834B',ls='--',lw=1.6,zorder=2,
   label=f"Tpc = {ref['pseudocritical_temperature_C']:.2f} °C at {ref['reference_pressure_mean_MPa']:.3f} MPa")
  gmin,gmax=r.mass_flux_kg_m2_s.min(),r.mass_flux_kg_m2_s.max()
  q=r.nominal_heat_flux_W_m2.iloc[0]/1000
  ax.set_title(f'{group} | G = {gmin:g}–{gmax:g} kg/(m²·s) | nominal q″ = {q:g} kW/m²',fontsize=11,pad=10)
  ax.set(xlim=(25,39.5),ylim=(0,.1) if compare else (0,.70),xlabel='Reported mean CO₂ bulk temperature (°C)',ylabel='Fitted effective e/D')
  ax.set_xticks([26,28,30,32,34,36,38]);ax.set_yticks(np.linspace(0,.1,6) if compare else np.arange(0,.71,.1));ax.grid(alpha=.18)
  ax.legend(loc='upper left',fontsize=8,framealpha=.9)
  if not compare:
   peak=r.loc[r.old_fitted_eD.idxmax()]
   stats.append({'source_group':group,'n_original':len(r),'n_corrected':int(r.new_fitted_eD.notna().sum()),'Spearman_temperature_vs_original_eD':float(spearmanr(r.bulk_CO2_C,r.old_fitted_eD).statistic),'Pearson_temperature_vs_original_eD':float(pearsonr(r.bulk_CO2_C,r.old_fitted_eD).statistic),'largest_original_eD':peak.old_fitted_eD,'temperature_at_largest_eD_C':peak.bulk_CO2_C,'peak_case':peak.Case})
 fig.suptitle('Fitted e/D vs bulk CO₂ temperature for each operating-condition class',fontsize=16,y=.98)
 subtitle='Original 40 DH-based fits' if not compare else 'Revised smooth-corrected fits only (blue; 24 matched cases); y-axis: 0–0.1'
 fig.text(.5,.935,subtitle,ha='center',fontsize=11)
 note='Effective e/D is fitted to HTC, not measured roughness. Points show individual runs; no trend law is fitted.'
 if compare:note+=' All 24 fitted values are within the displayed range; 16 revised cases are unfitted.'
 fig.text(.5,.039,note,ha='center',fontsize=8.5)
 fig.text(.5,.018,'Green dashed line: CO₂ Cp maximum at each class mean inlet pressure (HEOS); individual runs have pressure-dependent Tpc.',ha='center',fontsize=8.5)
 fig.tight_layout(rect=(0,.055,1,.91))
 stem='eD_vs_Bulk_Temperature_Corrected_Blue_0_to_01' if compare else 'eD_vs_Bulk_Temperature_Original_40'
 for extension in ['png','pdf']:fig.savefig(OUT/(stem+'.'+extension),dpi=220,bbox_inches='tight')
 figures.append(fig)
with PdfPages(OUT/'eD_vs_Bulk_Temperature_by_Class.pdf') as pdf:
 for fig in figures:pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
pd.DataFrame(list(references.values())).to_csv(OUT/'Pseudocritical_temperature_references.csv',index=False)
pd.DataFrame(stats).to_csv(OUT/'Temperature_by_class_statistics.csv',index=False)
print(pd.DataFrame(stats).to_string(index=False))
