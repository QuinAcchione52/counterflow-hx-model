"""Fit dimpled cases independently, keeping assumed water conditions explicit."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
from scipy.optimize import brentq
from CounterHX_full import match_bulk_case


def calibrate_one(row_dict, segments, max_ep):
    row=pd.Series(row_dict)
    target=float(row.target_htc_W_m2K)
    cache={}
    def evaluate(ep):
        if ep not in cache:
            cache[ep]=match_bulk_case(row,ep,float(row.target_mean_CO2_C),segments)
        return float(cache[ep].attrs['summary']['htc_spencer_W_m2K'])-target
    summary={'Case':row.Case,'target_htc_W_m2K':target,
             'target_mean_CO2_C':row.target_mean_CO2_C,
             'experimental_CO2_flow_kg_s':row['m_CO2 (kg/s)'],
             'experimental_pressure_MPa':row['P (Mpa)'],
             'assumed_water_inlet_C':row['T_H2O (C)'],
             'assumed_water_flow_kg_s':row['m_H2O (kg/s)'],
             'experimental_htc_uncertainty_percent':row.experimental_htc_uncertainty_percent,
             'nominal_heat_flux_W_m2':row.nominal_heat_flux_W_m2,
             'segments':segments,
             'baseline_model':row.get('CO2_baseline_model','none'),
             'baseline_model_extrapolated':row.get('CO2_baseline_model_extrapolated',False)}
    try:
        lo=evaluate(0.)
        physical_cap=evaluate(min(.05,max_ep))
        summary.update(smooth_model_htc_W_m2K=lo+target,
                       tested_low_cap_eD=min(.05,max_ep),
                       htc_at_tested_low_cap_W_m2K=physical_cap+target)
        if max_ep>=.05:
            summary['htc_at_eD_005_W_m2K']=physical_cap+target
        if lo*physical_cap<=0:
            bounds=(0.,min(.05,max_ep))
        else:
            high=evaluate(max_ep)
            summary['htc_at_max_eD_W_m2K']=high+target
            if lo*high>0:
                summary['status']='not_bracketed'
                summary['reason']='target outside endpoint predictions; no roughness fit reported'
                return summary,None
            bounds=(min(.05,max_ep),max_ep)
        ep=brentq(evaluate,*bounds,xtol=1e-7)
        evaluate(ep)
        profile=cache[ep]
        rel=float(profile.attrs['summary']['htc_spencer_W_m2K'])/target-1
        if abs(rel)>.001:
            raise RuntimeError('Fit fails 0.1% HTC acceptance tolerance')
        summary.update(profile.attrs['summary'])
        summary.update(effective_e_over_D=ep,relative_htc_error=rel,
                       status='fit_within_0_to_005' if ep<=.05 else 'extrapolated_fit')
        summary['predicted_to_nominal_flux_ratio']=summary['heat_flux_CO2_W_m2']/row.nominal_heat_flux_W_m2
        return summary,profile
    except Exception as exc:
        summary.update(status='solver_error',reason=f'{type(exc).__name__}: {exc}')
        return summary,None


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default='Dimpled_cases.csv')
    parser.add_argument('--output',default='Dimpled-batch-output.xlsx')
    parser.add_argument('--segments',type=int,default=20,help='Screening grid; confirm important fits at 100 segments')
    parser.add_argument('--workers',type=int,default=1,help='Separate processes; property objects cannot be shared between threads')
    parser.add_argument('--ep-max',type=float,default=.05)
    parser.add_argument('--allow-extrapolation',action='store_true')
    parser.add_argument('--smooth-baseline-model', help='Optional independently smooth-trained baseline JSON')
    parser.add_argument('--allow-baseline-extrapolation',action='store_true')
    parser.add_argument('--case',help='Run only one named case')
    args=parser.parse_args()
    if not 0<args.ep_max<=1:parser.error('0 < --ep-max <= 1 required')
    if args.ep_max>.05 and not args.allow_extrapolation:parser.error('Bounds above 0.05 require --allow-extrapolation')
    if args.workers<1:parser.error('--workers must be positive')
    data=pd.read_csv(args.input)
    if args.case:
        data=data[data.Case==args.case]
    if data.empty:parser.error('No matching input cases')
    if args.smooth_baseline_model:
        from smooth_baseline import apply_model
        data=apply_model(data,args.smooth_baseline_model,args.allow_baseline_extrapolation)
    summaries=[]; profiles={}
    output=Path(args.output)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(calibrate_one,row.to_dict(),args.segments,args.ep_max):row.Case
                 for _,row in data.iterrows()}
        for future in as_completed(futures):
            summary,profile=future.result()
            summaries.append(summary)
            if profile is not None:profiles[summary['Case']]=profile
            pd.DataFrame(summaries).to_csv(output.with_suffix('.partial.csv'),index=False)
            print(f"{len(summaries)}/{len(data)} {summary['Case']}: {summary['status']}; "
                  f"e/D={summary.get('effective_e_over_D', 'unfitted')}",flush=True)
    order={case:i for i,case in enumerate(data.Case)}
    summaries.sort(key=lambda r:order[r['Case']])
    table=pd.DataFrame(summaries)
    band=table.target_htc_W_m2K*table.experimental_htc_uncertainty_percent/100
    table['target_HTC_lower_reported_band_W_m2K']=table.target_htc_W_m2K-band
    table['target_HTC_upper_reported_band_W_m2K']=table.target_htc_W_m2K+band
    for field,label in [('smooth_model_htc_W_m2K','smooth_prediction_within_reported_band'),
                        ('htc_at_eD_005_W_m2K','eD005_prediction_within_reported_band')]:
        if field in table:
            table[label]=(table[field]-table.target_htc_W_m2K).abs()<=band
    with pd.ExcelWriter(output,engine='xlsxwriter') as writer:
        table.to_excel(writer,sheet_name='Summary',index=False)
        data.to_excel(writer,sheet_name='Inputs_assumptions',index=False)
        for case in data.Case:
            if case in profiles:profiles[case].to_excel(writer,sheet_name=str(case)[:31],index=False)
    table.to_csv(output.with_suffix('.summary.csv'),index=False)
    print(table.status.value_counts().to_string(),flush=True)


if __name__=='__main__':main()
