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
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/dimpled_trends'
d=pd.read_csv(OUT/'All_40_old_vs_new_eD.csv')
assert len(d)==40 and d.Case.is_unique
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
figures=[];stats=[]
for compare in [False,True]:
 fig,axs=plt.subplots(2,2,figsize=(12.5,8),sharex=True,sharey=True)
 for ax,(group,r) in zip(axs.ravel(),d.groupby('test_group',sort=True)):
  ax.scatter(r.bulk_CO2_C,r.old_fitted_eD,color='red',s=49,label=f'Original fits (n={len(r)})',zorder=3)
  if compare:
   p=r.dropna(subset=['new_fitted_eD'])
   ax.scatter(p.bulk_CO2_C,p.new_fitted_eD,color='blue',s=42,label=f'Corrected fits (n={len(p)})',zorder=4)
  ax.axhline(.05,color='#666666',ls='--',lw=1,label='Reference e/D = 0.05')
  gmin,gmax=r.mass_flux_kg_m2_s.min(),r.mass_flux_kg_m2_s.max()
  q=r.nominal_heat_flux_W_m2.iloc[0]/1000
  ax.set_title(f'{group} | G = {gmin:g}–{gmax:g} kg/(m²·s) | nominal q″ = {q:g} kW/m²',fontsize=11,pad=10)
  ax.set(xlim=(25,39.5),ylim=(0,.70),xlabel='Reported mean CO₂ bulk temperature (°C)',ylabel='Fitted effective e/D')
  ax.set_xticks([26,28,30,32,34,36,38]);ax.set_yticks(np.arange(0,.71,.1));ax.grid(alpha=.18)
  ax.legend(loc='upper left',fontsize=8,framealpha=.9)
  if not compare:
   peak=r.loc[r.old_fitted_eD.idxmax()]
   stats.append({'source_group':group,'n_original':len(r),'n_corrected':int(r.new_fitted_eD.notna().sum()),'Spearman_temperature_vs_original_eD':float(spearmanr(r.bulk_CO2_C,r.old_fitted_eD).statistic),'Pearson_temperature_vs_original_eD':float(pearsonr(r.bulk_CO2_C,r.old_fitted_eD).statistic),'largest_original_eD':peak.old_fitted_eD,'temperature_at_largest_eD_C':peak.bulk_CO2_C,'peak_case':peak.Case})
 fig.suptitle('Fitted e/D vs bulk CO₂ temperature for each operating-condition class',fontsize=16,y=.98)
 subtitle='Original 40 DH-based fits' if not compare else 'Original 40 DH-based fits (red) and revised smooth-corrected fits (blue; 24 matched cases)'
 fig.text(.5,.935,subtitle,ha='center',fontsize=11)
 note='Effective e/D is fitted to HTC, not measured roughness. Points show individual runs; no trend law is fitted.'
 if compare:note+=' The 16 unfitted revised cases have no blue e/D point.'
 fig.text(.5,.025,note,ha='center',fontsize=9)
 fig.tight_layout(rect=(0,.055,1,.91))
 stem='eD_vs_Bulk_Temperature_Original_and_Corrected' if compare else 'eD_vs_Bulk_Temperature_Original_40'
 for extension in ['png','pdf']:fig.savefig(OUT/(stem+'.'+extension),dpi=220,bbox_inches='tight')
 figures.append(fig)
with PdfPages(OUT/'eD_vs_Bulk_Temperature_by_Class.pdf') as pdf:
 for fig in figures:pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
pd.DataFrame(stats).to_csv(OUT/'Temperature_by_class_statistics.csv',index=False)
print(pd.DataFrame(stats).to_string(index=False))
