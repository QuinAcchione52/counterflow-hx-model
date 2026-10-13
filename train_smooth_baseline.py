"""Train an optional DH correction and evaluate leave-one-case-out predictions."""
import argparse,json,os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from CounterHX_full import match_bulk_case
from smooth_baseline import FEATURES,features,fit_coefficients


def fit_factor(values,segments):
    row=pd.Series(values);cache={};target=float(row.target_htc_W_m2K)
    row['CO2_heat_friction_model']='dang_hihara'
    def evaluate(factor):
        if factor not in cache:
            trial=row.copy();trial['CO2_baseline_correction_multiplier']=factor
            cache[factor]=match_bulk_case(trial,0.,row.target_mean_CO2_C,segments).attrs['summary']
        return cache[factor]['htc_spencer_W_m2K']-target
    factor=brentq(evaluate,.4,8.,xtol=1e-7)
    return {'Case':row.Case,'factor':factor,'target':target,'baseline':evaluate(1.)+target,'status':'fit'}


def predict_case(values,factor,segments):
    row=pd.Series(values);row['CO2_heat_friction_model']='dang_hihara'
    row['CO2_baseline_correction_multiplier']=factor
    summary=match_bulk_case(row,0.,row.target_mean_CO2_C,segments).attrs['summary']
    h=summary['htc_spencer_W_m2K'];target=float(row.target_htc_W_m2K)
    return {'Case':row.Case,'target_htc_W_m2K':target,'heldout_numerical_htc_W_m2K':h,
            'heldout_factor':factor,'heldout_error_percent':100*(h/target-1),**summary}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default='Smooth_cases.csv')
    parser.add_argument('--factors',help='Reuse fitted factors from the same code, inputs and segment count')
    parser.add_argument('--segments',type=int,default=20)
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--output-dir',default='results/smooth_baseline')
    args=parser.parse_args()
    if os.environ.get('HX_PROPERTY_BACKEND','HEOS')!='HEOS':parser.error('This calibration uses HEOS; remove HX_PROPERTY_BACKEND override')
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    data=pd.read_csv(args.input)
    if args.factors:
        fitted=pd.read_csv(args.factors)
    else:
        records=[]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            tasks=[pool.submit(fit_factor,row.to_dict(),args.segments) for _,row in data.iterrows()]
            for task in as_completed(tasks):
                records.append(task.result());print('Fitted',len(records),flush=True)
        fitted=pd.DataFrame(records)
    fitted=fitted.set_index('Case').loc[data.Case].reset_index()
    if len(fitted)!=len(data) or not fitted.status.eq('fit').all():raise ValueError('All smooth factors must be fitted')
    fitted[['Case','factor','target','baseline','status']].to_csv(out/'training_factors.csv',index=False)
    feature_rows=[features(row) for _,row in data.iterrows()]
    X=np.array([item[0] for item in feature_rows]);factors=fitted.factor.to_numpy()
    fixed_keys=['D (mm)','L (mm)','wall thickness (mm)','water gap (mm)','Kp','m_H2O (kg/s)','T_H2O (C)','P (atm)','T_amb (C)','Q_leak_slope','Q_leak_intercept','U_amb','d_outer_mm']
    if any(data[key].nunique(dropna=False)!=1 for key in fixed_keys):
        raise ValueError('This baseline requires shared fixed geometry, water and leakage assumptions')
    coefficients=fit_coefficients(X,factors)
    domain={key:[min(item[1][key] for item in feature_rows),max(item[1][key] for item in feature_rows)] for key in feature_rows[0][1]}
    model={'schema':'smooth-DH-global-factor-v1','backend':'HEOS','features':FEATURES,
           'coefficients':coefficients.tolist(),'domain':domain,'training_cases':len(data),
           'fixed_inputs':{key:float(data.iloc[0][key]) for key in fixed_keys},
           'segments':args.segments,'water_inlet_C':20.7,'water_mass_flow_kg_s':.087,
           'interpretation':'Empirical correction conditional on assumed water conditions; not physical roughness or independent correlation validation.'}
    records=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        tasks=[]
        for i,(_,row) in enumerate(data.iterrows()):
            keep=np.arange(len(data))!=i
            factor=float(np.exp(X[i]@fit_coefficients(X[keep],factors[keep])))
            tasks.append(pool.submit(predict_case,row.to_dict(),factor,args.segments))
        for task in as_completed(tasks):
            records.append(task.result());print('Held-out',len(records),flush=True)
    result=pd.DataFrame(records).set_index('Case').loc[data.Case].reset_index()
    result['uncorrected_htc_W_m2K']=fitted.baseline.to_numpy()
    result['uncorrected_error_percent']=100*(result.uncorrected_htc_W_m2K/result.target_htc_W_m2K-1)
    result.to_csv(out/'heldout_predictions.csv',index=False)
    with pd.ExcelWriter(out/'heldout_predictions.xlsx',engine='xlsxwriter') as writer:result.to_excel(writer,index=False,sheet_name='Held-out smooth cases')
    model['validation']={'method':'Leave-one-case-out; fitted factor for held-out case excluded from regression',
        'uncorrected_MAE_percent':float(result.uncorrected_error_percent.abs().mean()),
        'heldout_MAE_percent':float(result.heldout_error_percent.abs().mean()),
        'heldout_within_5_percent':int(result.heldout_error_percent.abs().le(5).sum()),
        'heldout_max_abs_error_percent':float(result.heldout_error_percent.abs().max())}
    (out/'model.json').write_text(json.dumps(model,indent=2)+'\n')
    print(json.dumps(model['validation'],indent=2),flush=True)

if __name__=='__main__':main()
