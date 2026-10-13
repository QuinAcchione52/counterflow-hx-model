import math
import numpy as np
import pytest
from roughness import build_roughness
from Correlations import DHC,colebrook,heat_friction,hydraulic_friction


def test_profiles_in_physical_coordinates_and_one_step():
    x=np.linspace(0,1,11)
    np.testing.assert_allclose(build_roughness(x,.02,'Step',1),.02)
    np.testing.assert_allclose(build_roughness(x,.02,'Linear'),.02*x)
    coarse=build_roughness(x,.02,'Step',3)
    fine=build_roughness(np.linspace(0,1,101),.02,'Step',3)
    np.testing.assert_allclose(coarse,fine[::10])
    custom=np.linspace(.001,.01,11)
    np.testing.assert_array_equal(build_roughness(x,custom),custom)
    assert np.all(build_roughness(x,.02,'Sine',3)>=0)
    for ep in [float('nan'),-.001,1.1]:
        with pytest.raises(ValueError):build_roughness(x,ep)
    with pytest.raises(ValueError):build_roughness(x,.02,'Step',0)
    with pytest.raises(ValueError):build_roughness(x,.02,'Step',1.2)
    with pytest.raises(ValueError):build_roughness(x,[.01,.02])


def test_colebrook_reference_point_and_implicit_residual():
    # Standard smooth-pipe reference at Re=1e5: Darcy f approximately .01799.
    assert colebrook(1e5,0)==pytest.approx(.0179898,rel=1e-5)
    for ep in [0,.001,.01,.05]:
        f=colebrook(1e5,ep)
        assert 1/math.sqrt(f)+2*math.log10(ep/3.7+2.51/(1e5*math.sqrt(f)))==pytest.approx(0,abs=1e-10)
    assert hydraulic_friction(1000,.01)==pytest.approx(.064)


def test_rough_extension_recovers_original_dh_and_increases_friction():
    for re in [4000,1e4,1e5]:
        original=(1.82*math.log10(re)-1.64)**-2
        assert heat_friction(re,0,'dh_roughness')==pytest.approx(original)
        assert heat_friction(re,.01,'dh_roughness')>original
        with pytest.raises(ValueError):heat_friction(re,.01,'dang_hihara')


def test_original_dh_zero_flux_matches_independent_equation_evaluation():
    import CoolProp.CoolProp as CP
    # Equal wall/bulk temperature removes all property-averaging choices.
    T=309.15;P=8e6;D=.006;G=200.
    state=CP.AbstractState('HEOS','CO2');state.update(CP.PT_INPUTS,P,T)
    re=G*D/state.viscosity();pr=state.cpmass()*state.viscosity()/state.conductivity()
    f=(1.82*math.log10(re)-1.64)**-2
    nu=(f/8*(re-1000)*pr)/(1.07+12.7*math.sqrt(f/8)*(pr**(2/3)-1))
    p=DHC(T,T,P,D,G,0,'CO2',0,0,.02,friction_model='dang_hihara').iloc[0]
    assert p.htc==pytest.approx(nu*state.conductivity()/D,rel=1e-10)
    assert p.Tw==T
    assert p.Pr_branch=='bulk_cp'


def test_invalid_temperature_is_not_silently_clipped():
    with pytest.raises(ValueError):DHC(-10,300,8e6,.006,200,-1000,'CO2',0,0,.02)
    with pytest.raises(ValueError,match='turbulent'):DHC(300,299,8e6,.006,1,-1000,'CO2',0,0,.02)


def test_step_profile_satisfies_conjugate_wall_balance():
    import pandas as pd
    from CounterHX_full import match_bulk_case
    row=pd.read_csv('Test.csv').iloc[0]
    row['Profile (Constant,Step,Sine,Linear)']='Step'
    row['Steps/Harmonic (n)']=3
    p=match_bulk_case(row,.05,35.23,40)
    scale=np.maximum(1.,np.abs(p.H2O_node_heat_flux_W_m2))
    assert np.max(np.abs(p.H2O_final_wall_flux_residual_W_m2)/scale)<1e-5
    assert np.max(np.abs(p.CO2_final_wall_flux_residual_W_m2)/np.maximum(1.,np.abs(p.CO2_node_heat_flux_W_m2)))<1e-5
    assert np.max(np.abs(p.wall_conduction_residual_K))<1e-4
    assert abs(p.attrs['summary']['CO2_energy_residual_W'])<1e-5
    assert abs(p.attrs['summary']['water_energy_residual_W'])<1e-5


def test_cold_case_initialization_does_not_create_false_freezing_state():
    import pandas as pd
    from CounterHX_full import match_bulk_case
    row=pd.read_csv('Smooth_cases.csv').set_index('Case').loc['March_27_1_b'].copy()
    row['CO2_heat_friction_model']='dang_hihara'
    row['CO2_baseline_correction_multiplier']=.4
    p=match_bulk_case(row,0.,row.target_mean_CO2_C,10)
    assert p.H2Ow.min()>273.15
    assert np.all(p.CO2b>p.H2Ob)
    assert p.attrs['summary']['CO2_baseline_correction_multiplier']==.4
    assert abs(p.attrs['summary']['CO2_energy_residual_W'])<1e-5


def test_baseline_model_is_explicit_and_guards_domain_and_assumptions(monkeypatch):
    import pandas as pd
    from smooth_baseline import load_model,correction,apply_model
    monkeypatch.setenv('HX_PROPERTY_BACKEND','HEOS')
    data=pd.read_csv('Smooth_cases.csv');row=data.iloc[0].copy()
    path='results/smooth_baseline/model.json';model=load_model(path)
    factor,outside=correction(row,model)
    assert .2<factor<10 and not outside
    assert 'CO2_baseline_correction_multiplier' not in data
    corrected=apply_model(data.iloc[:1],path)
    assert corrected.iloc[0].CO2_baseline_correction_multiplier==pytest.approx(factor)
    trial=row.copy();trial['nominal_heat_flux_W_m2']=4000
    with pytest.raises(ValueError,match='Outside'):correction(trial,model)
    assert correction(trial,model,True)[1]
    trial=row.copy();trial['T_H2O (C)']=25
    with pytest.raises(ValueError,match='training assumption'):correction(trial,model,True)
    monkeypatch.setenv('HX_PROPERTY_BACKEND','REFPROP')
    with pytest.raises(ValueError,match='backend'):correction(row,model)
