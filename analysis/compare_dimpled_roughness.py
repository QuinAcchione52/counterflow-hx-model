"""Compare archived and revised dimpled e/D fits; describe operating-condition trends.
No solver fits or experimental observations are changed by this analysis.
"""
from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/hx-trends-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr,pearsonr
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from smooth_baseline import pseudocritical_temperature, correction, load_model
OUT=ROOT/'results/dimpled_trends';OUT.mkdir(parents=True,exist_ok=True)
old=pd.read_csv(ROOT/'results/Spencer-vs-Numerical-HTC-All-40-Cases.csv')
new=pd.read_csv(ROOT/'results/smooth_baseline/Dimpled-transfer.summary.csv')
inputs=pd.read_csv(ROOT/'Dimpled_cases.csv')
ref=pd.read_csv(ROOT/'experimental_targets.csv')
assert all(d.Case.is_unique for d in [old,new,inputs,ref])
assert set(old.Case)==set(new.Case)==set(inputs.Case)==set(ref.Case)
d=pd.DataFrame({'Case':old.Case}).merge(inputs,on='Case',validate='one_to_one').merge(
 ref[['Case','mass_flux_kg_m2_s']],on='Case',validate='one_to_one')
a=old.set_index('Case').loc[d.Case];b=new.set_index('Case').loc[d.Case]
for name,values in {
 'old_fitted_eD':a['Fitted effective e/D'],
 'new_fitted_eD':b.effective_e_over_D,
 'new_smooth_HTC_W_m2K':b.smooth_model_htc_W_m2K,
 'new_cap005_HTC_W_m2K':b.htc_at_eD_005_W_m2K,
 'new_fitted_HTC_W_m2K':b.htc_spencer_W_m2K,
 'old_fitted_HTC_W_m2K':a['Numerical HTC (W/m² K)'],
 'old_numerical_flux_W_m2':a['Numerical heat flux (W/m²)'],
 'new_numerical_flux_W_m2':b.heat_flux_CO2_W_m2,
 'smooth_baseline_multiplier':b.CO2_baseline_correction_multiplier,
 'baseline_extrapolated':b.baseline_model_extrapolated,
 }.items():d[name]=values.to_numpy()
baseline_model=load_model(ROOT/'results/smooth_baseline/model.json')
for index,row in d.iterrows():
 factor,outside=correction(row,baseline_model,allow_extrapolation=True)
 d.loc[index,'smooth_baseline_multiplier']=factor
 d.loc[index,'baseline_extrapolated']=outside
d['computed_mass_flux_kg_m2_s']=d['m_CO2 (kg/s)']/(np.pi*(d['D (mm)']/1000)**2/4)
d['bulk_CO2_C']=d.target_mean_CO2_C
d['delta_Tpc_K']=[r.bulk_CO2_C-(pseudocritical_temperature(r['P (Mpa)']*1e6)-273.15) for _,r in d.iterrows()]
d['mass_flow_group']=d.test_group.str.split('/').str[0]
d['new_status']=np.where(d.new_fitted_eD.notna(),'matched_0_to_005',np.where(d.target_htc_W_m2K<d.new_smooth_HTC_W_m2K,'target_below_smooth','target_above_005_cap'))
d['old_within_005']=d.old_fitted_eD.le(.05)
d['old_minus_new_eD']=d.old_fitted_eD-d.new_fitted_eD
d['eD_reduction_percent']=100*(1-d.new_fitted_eD/d.old_fitted_eD)
d['new_minus_old_HTC_W_m2K']=d.new_fitted_HTC_W_m2K-d.old_fitted_HTC_W_m2K
paired=d[d.new_fitted_eD.notna()].copy()
assert len(d)==40 and len(paired)==24 and d.new_fitted_eD.isna().sum()==16
assert (paired.new_fitted_eD.between(0,.05)).all()
# Correlations are descriptive and conditional on fit availability; no causal claims.
factors={'mass_flux_kg_m2_s':'Reported mass flux','nominal_heat_flux_W_m2':'Nominal flux group','bulk_CO2_C':'Mean CO2 temperature','delta_Tpc_K':'Temperature relative to Tpc','P (Mpa)':'Inlet pressure'}
rows=[]
for label,subset,ycol in [('Historical all 40',d,'old_fitted_eD'),('Historical paired 24',paired,'old_fitted_eD'),('New matched 24',paired,'new_fitted_eD')]:
 for col,name in factors.items():
  valid=subset[[col,ycol]].dropna();rho=spearmanr(valid[col],valid[ycol]).statistic;r=pearsonr(valid[col],valid[ycol]).statistic
  rows.append({'dataset':label,'factor':name,'input_column':col,'n':len(valid),'Spearman_rho':rho,'Pearson_r':r})
corr=pd.DataFrame(rows)
records=[]
for group,sub in d.groupby('test_group',sort=True):
 p=sub.dropna(subset=['new_fitted_eD'])
 records.append({'test_group':group,'cases':len(sub),'old_within_005':int(sub.old_within_005.sum()),'new_matched':len(p),
  'new_below_smooth':int(sub.new_status.eq('target_below_smooth').sum()),'new_above_cap':int(sub.new_status.eq('target_above_005_cap').sum()),
  'mass_flux_min':sub.mass_flux_kg_m2_s.min(),'mass_flux_max':sub.mass_flux_kg_m2_s.max(),
  'bulk_T_min_C':sub.bulk_CO2_C.min(),'bulk_T_max_C':sub.bulk_CO2_C.max(),
  'old_eD_median_all':sub.old_fitted_eD.median(),'old_eD_median_paired':p.old_fitted_eD.median(),
  'new_eD_median_matched':p.new_fitted_eD.median()})
groups=pd.DataFrame(records)
transitions=pd.crosstab(d.old_within_005,d.new_status).reset_index().rename(columns={'old_within_005':'historical_eD_within_005'})
transitions.to_csv(OUT/'Fit_status_transitions.csv',index=False)
bins=[]
for label,sub,y in [('Historical all 40',d,'old_fitted_eD'),('New matched 24',paired,'new_fitted_eD')]:
 for group,part in sub.groupby('test_group'):
  for tempbin,chunk in part.groupby(pd.cut(part.bulk_CO2_C,[20,30,34,37,45]),observed=True):
   bins.append({'dataset':label,'test_group':group,'temperature_bin_C':str(tempbin),'n':len(chunk),'median_eD':chunk[y].median()})
bins=pd.DataFrame(bins)
d.to_csv(OUT/'All_40_old_vs_new_eD.csv',index=False)
groups.to_csv(OUT/'Group_summary.csv',index=False);corr.to_csv(OUT/'Exploratory_correlations.csv',index=False)
notes=[
 ('Historical model','Archived direct-Colebrook thermal extension with no smooth-trained multiplier; historical fit bounds allowed e/D up to 1.'),
 ('New model','Revised DH-relative-roughness extension plus smooth-trained baseline multiplier; fitted e/D restricted to 0..0.05.'),
 ('Not a controlled single-change test','Comparison includes implementation changes AND empirical baseline correction.'),
 ('Missing new fits','13 targets below smooth baseline; 3 above e/D=.05 endpoint. Blank e/D is not zero. No fits outside tested bound are invented.'),
 ('Selection bias','New e/D trends use only 24 matched cases. Missing cases depend on target and operating conditions; these are not representative complete-sample correlations.'),
 ('Flux predictor','Heat flux on the trend x-axis is the nominal experimental 2000/3000 group, not exact measured flux. Numerical flux is a solver output.'),
 ('Mass flux','Plots use source-table reported mass flux; geometry-derived mass flux is also exported. m12/m15 are source group labels, not assumed exact mass fluxes.'),
 ('Built-in dependence','The smooth multiplier already uses temperature relative to Tpc, mass flux and nominal flux. New e/D trends partly reflect that correction and its residual errors.'),
 ('Interpretation','e/D values are fitted effective parameters. The 0.05 cap is a chosen tested bound, not proven physical roughness.'),
 ('Statistics','Pearson and Spearman coefficients are exploratory associations, not validated laws, causal effects or independent predictive validation.'),
 ('Sources','Model results snapshot da28165; results/Spencer-vs-Numerical-HTC-All-40-Cases.csv; results/smooth_baseline/Dimpled-transfer.summary.csv; Dimpled_cases.csv; experimental_targets.csv.'),
]
with pd.ExcelWriter(OUT/'Dimpled_eD_comparison_and_trends.xlsx',engine='xlsxwriter') as writer:
 for name,frame in [('All_40_comparison',d),('Group_summary',groups),('Fit_status_transitions',transitions),('Correlations',corr),('Temperature_bins',bins),('Read_me',pd.DataFrame(notes,columns=['Item','Explanation']))]:
  frame.to_excel(writer,sheet_name=name,index=False);ws=writer.sheets[name];ws.freeze_panes(1,1);ws.autofilter(0,0,len(frame),len(frame.columns)-1)
  for i,c in enumerate(frame.columns):ws.set_column(i,i,min(34,max(15,len(c)+2)))
  if name=='Read_me':ws.set_column(0,0,30);ws.set_column(1,1,110,writer.book.add_format({'text_wrap':True}));[ws.set_row(i,48) for i in range(1,len(frame)+1)]
# Standalone publication-style charts and a multi-page PDF report.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
colors={'m12/q2':'#2478b4','m12/q3':'#d97c25','m15/q2':'#008c85','m15/q3':'#9b5eaa'}
figures=[]
def output(fig,name):
 fig.savefig(OUT/(name+'.png'),bbox_inches='tight',facecolor='white');figures.append(fig)
fig,axs=plt.subplots(2,1,figsize=(14,8),sharex=True,gridspec_kw={'height_ratios':[1,1]});x=np.arange(40)
axs[0].bar(x,d.old_fitted_eD,color=[colors[g] for g in d.test_group]);axs[0].axhline(.05,color='black',ls='--',lw=1,label='Tested e/D=.05 bound');axs[0].set_ylabel('Historical fitted e/D');axs[0].legend()
for status,marker,color in [('matched_0_to_005','o','#008c85'),('target_below_smooth','v','#d97c25'),('target_above_005_cap','^','#be4f50')]:
 mask=d.new_status.eq(status);y=d.loc[mask,'new_fitted_eD'] if status=='matched_0_to_005' else np.full(mask.sum(),0 if status=='target_below_smooth' else .05)
 axs[1].scatter(x[mask],y,marker=marker,color=color,s=42,label={'matched_0_to_005':'Fitted (24)','target_below_smooth':'Target below smooth (13): symbol at bound, NOT a fit','target_above_005_cap':'Target above cap (3): symbol at bound, NOT a fit'}[status])
axs[1].axhline(.05,color='black',ls='--',lw=1);axs[1].set_ylim(-.006,.061);axs[1].set_ylabel('New fitted e/D or bound flag');axs[1].set_xticks(x,d.Case,rotation=90,fontsize=8);axs[1].legend(fontsize=8,loc='upper right')
fig.suptitle('All 40 dimpled cases: historical vs revised + smooth-corrected calibration',fontsize=15);fig.tight_layout();output(fig,'01_all_40_comparison')
fig,ax=plt.subplots(figsize=(8,6))
for g,part in paired.groupby('test_group'):ax.scatter(part.old_fitted_eD,part.new_fitted_eD,label=g,color=colors[g],s=50)
ax.axhline(.05,color='gray',ls='--');ax.set(xlabel='Historical fitted e/D',ylabel='New fitted e/D',title='Paired comparison: 24 matched cases only');ax.legend(title='Source group');fig.tight_layout();output(fig,'02_paired_eD')
fig,axs=plt.subplots(2,3,figsize=(14,8))
for row,(subset,ycol,label) in enumerate([(d,'old_fitted_eD','Historical: all 40'),(paired,'new_fitted_eD','New: matched 24 only')]):
 for col,(factor,xlabel) in enumerate([('mass_flux_kg_m2_s','Reported mass flux (kg/m²·s)'),('nominal_heat_flux_W_m2','Nominal experimental heat-flux group (W/m²)'),('bulk_CO2_C','Reported mean CO₂ temperature (°C)')]):
  ax=axs[row,col]
  for g,part in subset.groupby('test_group'):
   xx=part[factor].to_numpy().copy()
   if factor=='nominal_heat_flux_W_m2':xx=xx+np.linspace(-32,32,len(xx))
   ax.scatter(xx,part[ycol],color=colors[g],label=g,s=32,alpha=.85)
  ax.set(xlabel=xlabel,ylabel='Required fitted e/D',title=label)
  if col==1:ax.set_xticks([2000,3000]);ax.set_xlim(1800,3200)
  if row==0 and col==0:ax.legend(fontsize=8,title='Source group')
fig.suptitle('Required e/D vs operating conditions (nominal-flux points jittered for visibility)',fontsize=14)
fig.tight_layout();output(fig,'03_operating_condition_trends')
fig,axs=plt.subplots(2,2,figsize=(12,8))
for ax,(g,part) in zip(axs.ravel(),d.groupby('test_group')):
 ax.scatter(part.bulk_CO2_C,part.old_fitted_eD,color='#327ec0',label='Historical e/D',s=42)
 p=part.dropna(subset=['new_fitted_eD']);ax.scatter(p.bulk_CO2_C,p.new_fitted_eD,color='#008c85',label='New fitted e/D',s=42)
 ax.axhline(.05,color='gray',ls='--');ax.set(xlabel='Reported mean CO₂ temperature (°C)',ylabel='Required fitted e/D',title=f'{g}: new fits {len(p)}/{len(part)}');ax.legend(fontsize=8)
fig.suptitle('Temperature trends separated by source mass-flow / nominal-flux group',fontsize=14);fig.tight_layout();output(fig,'04_temperature_by_group')
fig,ax=plt.subplots(figsize=(9,5));matrix=corr.pivot(index='factor',columns='dataset',values='Spearman_rho')[['Historical all 40','Historical paired 24','New matched 24']]
im=ax.imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto');ax.set_xticks(range(3),matrix.columns);ax.set_yticks(range(len(matrix)),matrix.index)
for i in range(len(matrix)):
 for j in range(3):ax.text(j,i,f'{matrix.iloc[i,j]:+.2f}',ha='center',va='center',color='white' if abs(matrix.iloc[i,j])>.55 else 'black')
fig.colorbar(im,ax=ax,label='Spearman rank correlation');ax.set_title('Exploratory associations with fitted e/D; no causal interpretation');fig.tight_layout();output(fig,'05_correlation_heatmap')
fig,ax=plt.subplots(figsize=(9,5));xx=np.arange(len(groups));w=.23
ax.bar(xx-w,groups.old_within_005,w,label='Historical ≤ .05',color='#327ec0');ax.bar(xx,groups.new_matched,w,label='New matched ≤ .05',color='#008c85');ax.bar(xx+w,groups.new_below_smooth+groups.new_above_cap,w,label='New unfitted',color='#d97c25')
ax.set_xticks(xx,groups.test_group);ax.set(ylabel='Cases',title='Fit reachability by source group');ax.legend();fig.tight_layout();output(fig,'06_group_reachability')
with PdfPages(OUT/'Dimpled_eD_comparison_and_trends.pdf') as pdf:
 import textwrap
 for page,items in enumerate([notes[:5],notes[5:]],1):
  fig=plt.figure(figsize=(11.7,8.3));fig.text(.07,.91,'Dimpled e/D comparison and operating trends',fontsize=20,weight='bold')
  text='40 archived fits compared with 24 revised bounded fits — reading notes '+str(page)+'/2\n\n'+ '\n\n'.join(k+': '+v for k,v in items)
  fig.text(.07,.84,'\n'.join(textwrap.fill(line,112) for line in text.split('\n')),fontsize=11,va='top',linespacing=1.45);pdf.savefig(fig);plt.close(fig)
 for fig in figures:pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
summary={'cases':40,'paired_cases':24,'old_eD_median_all':float(d.old_fitted_eD.median()),'old_eD_median_paired':float(paired.old_fitted_eD.median()),'new_eD_median_matched':float(paired.new_fitted_eD.median()),'median_paired_reduction_percent':float(paired.eD_reduction_percent.median()),'new_below_smooth':13,'new_above_cap':3}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));print(groups.to_string(index=False));print(corr.to_string(index=False))
