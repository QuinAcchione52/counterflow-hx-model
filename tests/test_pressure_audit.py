import numpy as np
import pytest
from analysis.audit_dimpled_pressure import pressure_components


def test_uniform_pipe_matches_standard_darcy_reference_on_nonuniform_grid():
    # Re=1e5, standard smooth Darcy f=0.0179898; L=1m, D=.01m, G=100.
    v,a,f=pressure_components([0,.2,.7,1],[500]*4,[1e-5]*4,100,.01,0)
    assert v==pytest.approx(.0179898*(1/.01)*100**2/(2*500),rel=1e-5)
    assert a==0


def test_density_increase_gives_pressure_recovery_independent_of_grid():
    for n in [2,11,101]:
        rho=np.linspace(400,700,n)
        v,a,f=pressure_components(np.linspace(0,1,n),rho,np.full(n,1e-5),100,.01,0)
        assert a==pytest.approx(100**2*(1/700-1/400),rel=1e-12)
        assert a<0 and v>0


def test_invalid_pressure_geometry_or_profile_is_rejected():
    for G,D,ep in [(float('nan'),.01,0),(100,0,0),(100,.01,-.01)]:
        with pytest.raises(ValueError):pressure_components([0,1],[500,500],[1e-5,1e-5],G,D,ep)
    with pytest.raises(ValueError):pressure_components([0,0],[500,500],[1e-5,1e-5],100,.01,0)
