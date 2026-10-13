import math
import numpy as np
import pandas as pd
import pytest
from Correlations import DHC, laminar_annulus
from CounterHX_full import run_case, fit_case
from hx_metrics import logarithmic_mean, summarize


def test_laminar_annulus_parallel_plate_limit():
    # Inner heated, outer insulated parallel plates: Nu_Dh=70/13, f*Re=96.
    nu, poiseuille = laminar_annulus(.99)
    assert nu == pytest.approx(70/13, rel=.003)
    assert poiseuille == pytest.approx(96, rel=.001)


def test_endpoint_metric_is_not_arithmetic_htc():
    assert logarithmic_mean(2,2) == 2
    with pytest.raises(ValueError): logarithmic_mean(0,2)
    p=pd.DataFrame({'pos':[0.,1.], 'Q':[-100.,np.nan],
                    'CO2b':[310.,305.], 'CO2w':[308.,301.],
                    'H2Ob':[295.,294.], 'htc_CO2':[1000.,2000.],
                    'CO2_dp_cumulative_psi':[0.,.01], 'Re_H2O':[1000.,1000.]})
    result=summarize(p,10.,1.)
    assert result['htc_spencer_W_m2K']==pytest.approx(100/(math.pi*.01*logarithmic_mean(2,4)))
    assert result['htc_local_length_mean_W_m2K']==1500


def test_dhc_inlet_uses_same_physics_as_other_nodes():
    kwargs=dict(Tb=309.65,Tw=306.,P=8.094e6,D=.0133858,G=30.55,q=-3000,
                fluid='CO2',ep=.001,f=.02)
    a=DHC(count=0,**kwargs); b=DHC(count=1,**kwargs)
    assert a.htc.iloc[0] == pytest.approx(b.htc.iloc[0])
    assert a.Tw.iloc[0] < kwargs['Tb']
    assert a.Tw.iloc[0] == pytest.approx(kwargs['Tb']+kwargs['q']/a.htc.iloc[0],abs=2e-7)


@pytest.fixture(scope='module')
def solved():
    row=pd.read_csv('Test.csv').iloc[0]
    return row, run_case(row,.001,10)


def test_counterflow_boundaries_energy_and_walls(solved):
    row,p=solved; summary=p.attrs['summary']
    assert summary['converged']
    assert p.CO2b.iloc[0] == pytest.approx(309.65)
    assert p.H2Ob.iloc[-1] == pytest.approx(293.85,abs=2e-6)
    assert p.CO2w.iloc[0] < p.CO2b.iloc[0]
    assert np.all(p.CO2w > p.H2Ow)
    assert abs(summary['CO2_energy_residual_W'])<1e-5
    assert abs(summary['water_energy_residual_W'])<1e-5
    assert summary['water_Re_inlet'] < 2300
    assert summary['outer_diameter_mm']==pytest.approx(15.875)
    assert np.all(p.CO2b > p.H2Ob)


def test_fit_rejects_unreachable_target(monkeypatch,solved):
    row,p=solved
    monkeypatch.setattr('CounterHX_full.run_case',lambda *a,**k:p)
    with pytest.raises(ValueError,match='not bracketed'):fit_case(row,2642,segments=10)


def test_fit_recovers_known_roughness(monkeypatch,solved):
    row,p=solved
    def fake(row,ep,segments):
        result=p.copy(); result.attrs={'summary':{'htc_spencer_W_m2K':1000+10000*ep}}
        return result
    monkeypatch.setattr('CounterHX_full.run_case',fake)
    ep,result=fit_case(row,1200,segments=10)
    assert ep==pytest.approx(.02)


def test_inlet_inverse_matches_bulk_and_labels_assumption(monkeypatch,solved):
    from CounterHX_full import match_bulk_case
    row,p=solved
    def fake(row,ep,segments):
        result=p.copy()
        result.attrs={'summary':{'CO2_mean_endpoint_C':row['T_CO2 (C)']-.5}}
        return result
    monkeypatch.setattr('CounterHX_full.run_case',fake)
    result=match_bulk_case(row,.001,35.23,segments=10)
    assert result.attrs['summary']['inferred_CO2_inlet_C']==pytest.approx(35.73)
    assert result.attrs['summary']['CO2_mean_endpoint_C']==pytest.approx(35.23)
    assert result.attrs['summary']['CO2_inlet_is_inferred']


def test_large_roughness_fit_requires_explicit_extrapolation(solved):
    row,p=solved
    with pytest.raises(ValueError,match='allow_extrapolation'):
        fit_case(row,2642,bounds=(0,1),segments=10)
