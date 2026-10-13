"""Optional empirical smooth-tube correction; separate from DH and roughness."""
import json,math,os
from functools import lru_cache
from pathlib import Path
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar

FEATURES=['intercept','delta_Tpc_over_5K','delta_Tpc_squared','log_G_over_30','log_nominal_flux_over_2500']


@lru_cache(maxsize=256)
def pseudocritical_temperature(P,backend='HEOS'):
    state=CP.AbstractState(backend,'CO2')
    def negative_cp(T):state.update(CP.PT_INPUTS,P,float(T));return -state.cpmass()
    result=minimize_scalar(negative_cp,bounds=(298.,320.),method='bounded',options={'xatol':1e-8})
    if not result.success:raise RuntimeError('Pseudocritical temperature search failed')
    return float(result.x)


def features(row,backend='HEOS'):
    P=float(row['P (Mpa)'])*1e6
    T=float(row['target_mean_CO2_C'])
    G=float(row['m_CO2 (kg/s)'])/(math.pi*(float(row['D (mm)'])*1e-3)**2/4)
    q=float(row['nominal_heat_flux_W_m2'])
    if not all(math.isfinite(v) for v in [P,T,G,q]) or min(P,G,q)<=0:raise ValueError('Finite measured mean, pressure, flow, diameter and nominal flux required')
    delta=T-(pseudocritical_temperature(P,backend)-273.15)
    z=delta/5.
    X=np.array([1.,z,z*z,math.log(G/30.),math.log(q/2500.)])
    domain={'delta_Tpc_K':delta,'mass_flux_kg_m2s':G,'nominal_flux_W_m2':q,'pressure_MPa':P/1e6}
    return X,domain


def fit_coefficients(X,factors):
    factors=np.asarray(factors,dtype=float)
    if np.any(factors<=0) or not np.all(np.isfinite(factors)):raise ValueError('Positive finite calibration factors required')
    coefficients,_,rank,_=np.linalg.lstsq(np.asarray(X),np.log(factors),rcond=None)
    if rank!=len(FEATURES):raise ValueError('Calibration features are rank deficient')
    return coefficients


def correction(row,model,allow_extrapolation=False):
    backend=os.environ.get('HX_PROPERTY_BACKEND','HEOS')
    if backend!=model['backend']:raise ValueError('Baseline model backend differs from the running property backend')
    for key, expected in model.get('fixed_inputs', {}).items():
        if not math.isclose(float(row[key]), float(expected), rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError('Baseline model requires training assumption '+key)
    X,point=features(row,backend)
    outside=[key for key,value in point.items() if not model['domain'][key][0]-1e-8<=value<=model['domain'][key][1]+1e-8]
    if outside and not allow_extrapolation:raise ValueError('Outside smooth baseline calibration domain: '+', '.join(outside))
    factor=float(np.exp(X@np.asarray(model['coefficients'])))
    if not math.isfinite(factor) or not .2<=factor<=10:
        raise ValueError('Baseline factor outside numerical guard 0.2..10; investigate extrapolation')
    return factor,bool(outside)


def load_model(path):
    model=json.loads(Path(path).read_text())
    if model.get('schema')!='smooth-DH-global-factor-v1' or model.get('features')!=FEATURES:
        raise ValueError('Unsupported smooth baseline model schema')
    return model


def apply_model(data,path,allow_extrapolation=False):
    model=load_model(path);data=data.copy()
    if 'nominal_heat_flux_W_m2' not in data:
        raise ValueError('Baseline calibration requires nominal_heat_flux_W_m2; it does not impose this flux')
    for index,row in data.iterrows():
        factor,outside=correction(row,model,allow_extrapolation)
        data.loc[index,'CO2_baseline_correction_multiplier']=factor
        data.loc[index,'CO2_baseline_model_extrapolated']=outside
        data.loc[index,'CO2_baseline_model']=str(path)
    return data
