"""Audit dimpled tube pressure losses without fitting an HTC multiplier or pipe loss.
Saved thermal profiles are reused; unmatched cases use the corrected smooth
thermal baseline. Hydraulic roughness sensitivity is evaluated at fixed states.
"""
from pathlib import Path
import sys,os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/hx-pressure-mpl')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from Correlations import hydraulic_friction
from CounterHX_full import match_bulk_case
from concurrent.futures import ProcessPoolExecutor,as_completed
PSI_PA=6894.757293168
OUT=ROOT/'results/pressure_audit'


def pressure_components(x,rho,mu,G,D,eD):
    """SI Darcy losses; constant mass flux; positive drop in flow direction."""
    x,rho,mu=[np.asarray(v,dtype=float) for v in (x,rho,mu)]
    if x.ndim!=1 or len(x)<2 or x.shape!=rho.shape or x.shape!=mu.shape:
        raise ValueError('Matching 1D position/density/viscosity arrays required')
    if not all(np.isfinite(v).all() for v in [x,rho,mu]) or not np.isfinite([G,D]).all() or np.any(np.diff(x)<=0) or min(rho.min(),mu.min(),G,D)<=0:
        raise ValueError('Finite increasing positions and positive fluid properties, G, D required')
    if not np.isfinite(eD) or not 0<=eD<=1:raise ValueError('e/D must be finite in 0..1')
    f=np.array([hydraulic_friction(float(re),eD) for re in G*D/mu])
    visc=np.diff(x)/D*G*G/4*(f[:-1]/rho[:-1]+f[1:]/rho[1:])
    acc=G*G*np.diff(1/rho)
    return float(visc.sum()),float(acc.sum()),f


def state_profile(values,profile=None):
    row=pd.Series(values)
    if profile is None:
        row['CO2_baseline_correction_multiplier']=row.smooth_baseline_multiplier
        profile=match_bulk_case(row,0.,row.target_mean_CO2_C,20)
        kind='Corrected smooth thermal baseline; no fitted dimple HTC'
    else:
        profile=pd.DataFrame(profile);kind='Saved dimple HTC-fitted thermal profile'
    state=CP.AbstractState('HEOS','CO2');rho=[];mu=[]
    for T in profile.CO2b:
        state.update(CP.PT_INPUTS,row['P (Mpa)']*1e6,float(T));rho.append(state.rhomass());mu.append(state.viscosity())
    D=row['D (mm)']/1000.;G=row['m_CO2 (kg/s)']/(np.pi*D*D/4);x=profile.pos.to_numpy()
    rho=np.array(rho);mu=np.array(mu)
    v0,acc,f0=pressure_components(x,rho,mu,G,D,0.)
    v5,acc5,f5=pressure_components(x,rho,mu,G,D,.05)
    target=float(row.pressure_drop_psi)*PSI_PA
    # The coefficient multiplying a spatially constant apparent Darcy f.
    coefficient=float((np.diff(x)/D*G*G/4*(1/rho[:-1]+1/rho[1:])).sum())
    residual=target-v5-acc;dynamic=G*G/(2*rho[0])
    result={'Case':row.Case,'test_group':row.test_group,'thermal_state_source':kind,
      'target_mean_CO2_C':row.target_mean_CO2_C,'CO2_inlet_pressure_MPa':row['P (Mpa)'],
      'CO2_mass_flow_kg_s':row['m_CO2 (kg/s)'],'computed_mass_flux_kg_m2_s':G,
      'diameter_m':D,'length_m':x[-1]-x[0],'segments':len(x)-1,
      'thermal_eD':row.new_fitted_eD,'thermal_smooth_multiplier':row.smooth_baseline_multiplier,
      'target_pressure_psi':row.pressure_drop_psi,'target_pressure_Pa':target,
      'rho_inlet_kg_m3':rho[0],'rho_outlet_kg_m3':rho[-1],'Re_inlet':G*D/mu[0],
      'smooth_friction_Pa':v0,'acceleration_Pa':acc,'smooth_total_Pa':v0+acc,
      'smooth_total_psi':(v0+acc)/PSI_PA,'eD005_friction_Pa':v5,'eD005_total_Pa':v5+acc,
      'eD005_total_psi':(v5+acc)/PSI_PA,
      'smooth_error_percent':100*((v0+acc)/target-1),
      'eD005_error_percent':100*((v5+acc)/target-1),
      'target_to_eD005_ratio':target/(v5+acc),
      'target_bracketed_0_to_005':min(v0+acc,v5+acc)<=target<=max(v0+acc,v5+acc),
      'required_apparent_Darcy_f_if_all_drop_is_tube_friction':(target-acc)/coefficient,
      'modeled_Darcy_f_eD005_mean':float(np.mean(f5)),
      'unexplained_drop_at_eD005_Pa':residual,
      'equivalent_additional_K_inlet_reference_NOT_fitted':residual/dynamic,
      'equivalent_additional_pipe_length_m_at_eD005_NOT_measured':residual/(v5/(x[-1]-x[0])),
      'equivalent_CO2_static_head_m_NOT_measured':residual/(rho[0]*9.80665),
      'tube_span_and_pressure_units_verified_against_original_pdf':False,
      'property_backend':'HEOS'}
    if pd.notna(row.new_fitted_eD):
        vf,af,ff=pressure_components(x,rho,mu,G,D,float(row.new_fitted_eD))
        result['HTC_fit_pressure_psi']=(vf+af)/PSI_PA
        saved=float(profile.CO2_dp_cumulative_psi.iloc[-1])
        result['saved_pressure_psi']=saved
        result['recomputed_minus_saved_Pa']=(vf+af)-saved*6894.757
        if not np.isclose(vf+af,saved*6894.757,rtol=1e-10,atol=1e-9):
            raise AssertionError('Independent pressure integration disagrees with saved solver')
    profile=profile.copy();profile['audit_rho_bulk_kg_m3']=rho;profile['audit_mu_bulk_Pa_s']=mu
    profile['audit_f_smooth']=f0;profile['audit_f_eD005']=f5
    return result,profile


def main():
    if os.environ.get('HX_PROPERTY_BACKEND','HEOS')!='HEOS':raise ValueError('Audit requires HEOS')
    OUT.mkdir(parents=True,exist_ok=True)
    inputs=pd.read_csv(ROOT/'Dimpled_cases.csv')
    trends=pd.read_csv(ROOT/'results/dimpled_trends/All_40_old_vs_new_eD.csv')
    targets=pd.read_csv(ROOT/'experimental_targets.csv')
    data=inputs.merge(trends[['Case','smooth_baseline_multiplier','new_fitted_eD']],on='Case',validate='one_to_one').merge(targets[['Case','pressure_drop_psi']],on='Case',validate='one_to_one')
    workbook=pd.ExcelFile(ROOT/'results/smooth_baseline/Dimpled-transfer.xlsx')
    results=[];profiles={}
    with ProcessPoolExecutor(max_workers=3) as pool:
        tasks=[]
        for _,row in data.iterrows():
            profile=pd.read_excel(workbook,sheet_name=row.Case).to_dict('list') if row.Case in workbook.sheet_names else None
            tasks.append(pool.submit(state_profile,row.to_dict(),profile))
        for task in as_completed(tasks):
            result,profile=task.result();results.append(result);profiles[result['Case']]=profile
            print(len(results),result['Case'],result['eD005_error_percent'],flush=True)
    summary=pd.DataFrame(results).set_index('Case').loc[data.Case].reset_index()
    summary.to_csv(OUT/'All_40_pressure_audit.csv',index=False)
    notes=[
      ('Scope','Audit only: no pressure multiplier or minor-loss coefficient is fitted. No HTC targets or thermal fits are changed.'),
      ('Thermal states','Reuse 24 saved HTC-fitted profiles. For 16 unmatched cases solve the smooth thermal baseline with the existing smooth-trained multiplier and inferred mean. These are not HTC fits.'),
      ('Hydraulic sensitivity','Hydraulic e/D=0 and .05 are evaluated on the same fixed thermal profile for each case. This isolates pressure sensitivity; it is not a recoupled thermal roughness sweep.'),
      ('Equation','Darcy: integral f G²/(2 rho D) dx, trapezoid endpoint f/rho. Acceleration: G²(1/rho_out - 1/rho_in), signed. 1 psi = 6894.757293168 Pa.'),
      ('Measurement scope','Original reports are unavailable in current workspace. Experimental psi values are from the previously extracted table, not newly source-verified. Tap span, pressure-tap elevation, fittings and sensor zero remain unknown.'),
      ('Equivalent losses','Apparent f, added K, pipe length and static-head columns express the missing pressure scale under separate hypothetical explanations. They are not measured/fitted hardware parameters or diagnoses.'),
      ('Missing physics','No dimple-specific form-drag correlation or unverified fitting/tap loss is added. Dimple geometry is not equivalent sand-grain roughness.'),
      ('Pressure feedback','The HX currently evaluates properties at fixed input pressure; pressure loss is postprocessed. Small experimental ΔP relative to 8 MPa makes this unlikely to explain the large ratio, but a coupled pressure solve is not performed.'),
      ('Uncertainty','No per-case pressure uncertainty or zero-offset correction was available; HTC uncertainty is not substituted.'),
      ('Interpretation','A target outside hydraulic e/D=0..0.05 predictions is not forced to match. .05 is a tested bound, not validated dimple geometry.')]
    with pd.ExcelWriter(OUT/'Pressure_drop_audit_all_40.xlsx',engine='xlsxwriter') as writer:
        for name,frame in [('Read_me',pd.DataFrame(notes,columns=['Item','Explanation'])),('All_40_pressure',summary),('Inputs',data)]:
            frame.to_excel(writer,sheet_name=name,index=False);ws=writer.sheets[name];ws.freeze_panes(1,1);ws.autofilter(0,0,len(frame),len(frame.columns)-1);ws.set_column(0,len(frame.columns)-1,23)
            if name=='Read_me':ws.set_column(0,0,26);ws.set_column(1,1,110,writer.book.add_format({'text_wrap':True}));[ws.set_row(j,65) for j in range(1,len(frame)+1)]
        for case,profile in profiles.items():profile.to_excel(writer,sheet_name=case[:31],index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for group,sub in summary.groupby('test_group'):
        axs[0].scatter(sub.target_pressure_psi,sub.eD005_total_psi,label=group)
        axs[1].scatter(sub.target_mean_CO2_C,sub.target_to_eD005_ratio,label=group)
    lo=min(summary.eD005_total_psi.min(),summary.target_pressure_psi.min())*.7;hi=summary.target_pressure_psi.max()*1.4
    axs[0].plot([lo,hi],[lo,hi],'k--',lw=1,label='Exact agreement');axs[0].set(xscale='log',yscale='log',xlabel='Spencer ΔP (psi), extracted table',ylabel='Model ΔP (psi), hydraulic e/D=.05',title='40 cases: fixed thermal-state pressure sensitivity')
    axs[1].set(xlabel='Reported mean CO₂ temperature (°C)',ylabel='Experimental / model ΔP at hydraulic e/D=.05',title='Remaining pressure-drop ratio');axs[0].legend(fontsize=8);axs[1].legend(fontsize=8)
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(OUT/f'Pressure_comparison_all_40.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    statistics={'cases':len(summary),'targets_bracketed_0_to_005':int(summary.target_bracketed_0_to_005.sum()),'eD005_MAE_percent':float(summary.eD005_error_percent.abs().mean()),'eD005_experiment_to_model_ratio_min':float(summary.target_to_eD005_ratio.min()),'eD005_experiment_to_model_ratio_max':float(summary.target_to_eD005_ratio.max()),'saved_profile_max_reintegration_error_Pa':float(summary.recomputed_minus_saved_Pa.abs().max())}
    (OUT/'summary.json').write_text(json.dumps(statistics,indent=2)+'\n');print(json.dumps(statistics,indent=2),flush=True)

if __name__=='__main__':main()
