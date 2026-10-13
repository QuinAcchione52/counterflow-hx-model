# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 08:16:31 2024

@author: dhard

Property backend: HX_PROPERTY_BACKEND=HEOS (default) or REFPROP (licensed install).
CO2: Dang-Hihara Prandtl branches and film-Re friction; Colebrook roughness
is an effective extension, not a validated dimple correlation. Wall temperature
is iterated to convergence, including both endpoints.
Water: fully developed laminar inner-flux/outer-insulated annulus below Re=2300;
linear interpolation over 2300..4000; Gnielinski above 4000. The laminar
assumption excludes entrance and mixed-convection effects and is a baseline.
"""

import pandas as pd
import CoolProp.CoolProp as CP
import math
import os
import warnings
from functools import lru_cache
from scipy.integrate import quad

# One AbstractState per fluid, reused across calls
_AS = {}


def _state(fluid):
    backend = os.environ.get("HX_PROPERTY_BACKEND", "HEOS")
    key = (backend, fluid)
    AS = _AS.get(key)
    if AS is None:
        AS = CP.AbstractState(backend, fluid)
        _AS[key] = AS
    return AS


def _props_DHC(fluid, T, P):
    """rho, mu, Cp, k, i at (T, P) — same set the original PropsSI calls returned."""
    AS = _state(fluid)
    AS.update(CP.PT_INPUTS, P, T)
    return {'rho': AS.rhomass(), 'mu': AS.viscosity(), 'Cp': AS.cpmass(),
            'k': AS.conductivity(), 'i': AS.hmass()}


@lru_cache(maxsize=1)
def _warn_extrapolation():
    warnings.warn("e/D > 0.05 is an extrapolated calibration parameter, not validated physical roughness", RuntimeWarning, stacklevel=3)


def colebrook(Re, ep):
    """Darcy friction: solve the implicit turbulent Colebrook equation in 1/sqrt(f)."""
    if not math.isfinite(Re) or Re<=0 or not math.isfinite(ep) or not 0<=ep<=1:
        raise ValueError('Positive Reynolds number and finite e/D in 0..1 required')
    inv=7.
    for _ in range(100):
        updated=-2*math.log10(ep/3.7+2.51*inv/Re)
        if updated<=0:raise ValueError('Colebrook has no positive inverse-friction solution')
        if abs(updated-inv)<1e-11:return updated**-2
        inv=updated
    raise RuntimeError('Colebrook iteration did not converge')


def heat_friction(Re_f,ep,model='dh_roughness'):
    """Original DH smooth friction, or explicitly empirical roughness extensions.

    dh_roughness keeps the original smooth limit and uses the Colebrook rough /
    smooth friction ratio. This extension is not published Dang-Hihara physics.
    colebrook reproduces the previous direct replacement. dang_hihara is smooth only.
    """
    if not math.isfinite(Re_f) or Re_f<=0 or not math.isfinite(ep) or not 0<=ep<=1:
        raise ValueError('Heat friction requires positive finite film Re and e/D in 0..1')
    smooth=(1.82*math.log10(Re_f)-1.64)**-2
    if model=='dang_hihara':
        if ep!=0:raise ValueError('Original Dang-Hihara mode requires e/D=0')
        return smooth
    if model=='colebrook':return colebrook(Re_f,ep)
    if model=='dh_roughness':return smooth*colebrook(Re_f,ep)/colebrook(Re_f,0.)
    raise ValueError('CO2 heat friction must be dang_hihara, dh_roughness, or colebrook')


def hydraulic_friction(Re,ep):
    """Mechanical Darcy friction uses BULK Reynolds number, not film Reynolds."""
    if not math.isfinite(Re) or Re<=0:raise ValueError('Positive finite Reynolds number required')
    if Re<=2300:return 64/Re
    if Re>=4000:return colebrook(Re,ep)
    weight=(Re-2300)/1700
    return (1-weight)*64/Re+weight*colebrook(4000,ep)


def DHC(Tb,Tw,P,D,G,q,fluid,count,ep,f,friction_model='dh_roughness', baseline_multiplier=1.0):
    """DH equations as transcribed in Spencer's thesis Eqs.4–10.

    Temperatures K, pressure Pa, diameter m, mass flux kg/(m2 s), q W/m2
    into CO2 (negative for cooling). No inlet-node exception, no silent clamps.
    `count` and initial `f` are retained for caller compatibility.
    """
    if not all(math.isfinite(v) for v in [Tb,Tw,P,D,G,q,ep]):
        raise ValueError('DHC inputs must be finite')
    if min(Tb,Tw,P,D,G)<=0 or not 0<=ep<=1:
        raise ValueError('Positive temperatures, SI pressure, diameter and mass flux required; e/D in 0..1')
    if not math.isfinite(baseline_multiplier) or baseline_multiplier<=0:
        raise ValueError('Baseline correction multiplier must be positive and finite')
    if q>0:raise ValueError("Dang-Hihara is a cooling correlation; heating needs a different model")
    if ep>.05:_warn_extrapolation()
    bulk=_props_DHC(fluid,Tb,P)
    Reb=G*D/bulk['mu']
    if Reb<3000:
        raise ValueError(f'Dang-Hihara turbulent correlation cannot be used at bulk Re={Reb:.0f} (<3000)')
    def evaluate(wall_T):
        wall=_props_DHC(fluid,wall_T,P)
        Tf=(Tb+wall_T)/2
        film=_props_DHC(fluid,Tf,P)
        cpbar=bulk['Cp'] if abs(Tb-wall_T)<1e-6 else (bulk['i']-wall['i'])/(Tb-wall_T)
        if not math.isfinite(cpbar) or cpbar<=0:raise ValueError('Invalid integrated heat capacity')
        mkb=bulk['mu']/bulk['k'];mkf=film['mu']/film['k']
        if bulk['Cp']>=cpbar:
            Pr=bulk['Cp']*mkb;branch='bulk_cp'
        elif mkb>=mkf:
            Pr=cpbar*mkb;branch='integrated_cp_bulk_transport'
        else:
            Pr=cpbar*mkf;branch='integrated_cp_film_transport'
        Ref=G*D/film['mu'];ff=heat_friction(Ref,ep,friction_model)
        Nu=(ff/8*(Reb-1000)*Pr)/(1.07+12.7*math.sqrt(ff/8)*(Pr**(2/3)-1))
        Nu_before_correction=Nu
        Nu=baseline_multiplier*Nu
        h=Nu*film['k']/D
        if not math.isfinite(h) or h<=0:raise ValueError('Nonpositive or nonfinite DH HTC')
        return dict(Tw=wall_T,Tf=Tf,htc=h,Nu=Nu,f=ff,Re_bulk=Reb,Re_film=Ref,
                    Cp_bulk=bulk['Cp'],Cp_integrated=cpbar,Pr_used=Pr,k_film=film['k'],
                    Pr_branch=branch,heat_friction_model=friction_model,
                    baseline_multiplier=baseline_multiplier,Nu_before_baseline_correction=Nu_before_correction)
    for iterations in range(150):
        values=evaluate(Tw)
        next_T=Tb+q/values['htc']
        if not math.isfinite(next_T) or next_T<=0:raise ValueError('Unphysical wall temperature from q/h')
        if abs(next_T-Tw)<1e-7:
            # Return all properties evaluated at the actual returned wall temperature.
            values=evaluate(next_T)
            values['wall_iterations']=iterations+1
            values['wall_flux_residual_W_m2']=values['htc']*(next_T-Tb)-q
            if abs(values['wall_flux_residual_W_m2'])<=max(1e-4,abs(q)*1e-7):
                return pd.DataFrame([values])
        Tw=.5*Tw+.5*next_T
    raise RuntimeError('Dang-Hihara wall iteration did not converge')


@lru_cache(maxsize=64)
def laminar_annulus(a):
    """Fully developed laminar Nu: uniform inner-wall flux, insulated outer wall.

    Integrate the annular Poiseuille velocity and radial energy equation.
    Radius is normalized by outer radius; theta is temperature up to a scale.
    The returned Darcy Poiseuille number is f*Re based on hydraulic diameter.
    No turbulent floor or empirical entrance correction is applied.
    """
    if not 0 < a < 1:
        raise ValueError("Annulus diameter ratio must be between zero and one")
    c = (1-a*a)/math.log(a)
    def velocity(r):
        return 1-r*r-c*math.log(r)
    def primitive(r):
        return r*r/2-r**4/4-c*(r*r*math.log(r)/2-r*r/4)
    def flux(r):
        return primitive(1)-primitive(r)
    def theta(r):
        return quad(lambda t: flux(t)/t, r, 1, epsabs=1e-11)[0]
    integral = flux(a)
    mean_theta = quad(lambda r: r*velocity(r)*theta(r), a, 1, epsabs=1e-11)[0]/integral
    nu = 2*(1-a)*integral/(a*(theta(a)-mean_theta))
    mean_velocity = 2*integral/(1-a*a)
    poiseuille = 32*(1-a)**2/mean_velocity
    return nu, poiseuille

def Gniel(Tb,Tw,P,a,G,dh,L,q,fluid):
    """
    Heat transfer to cocentric annulous, this specific form is heat transfer
    to the inner wall assuming the outter wall is perfectly insulated

    Parameters
    ----------
    Tb : Bulk temperature of annular fluid
    Tw : Wall temperature of inner wall in annulus
    P : Pressure of annular fluid
    a : Diameter ratio
    G : Mass flux of cooling fluid
    dh : Hydraulic Diameter of annulus
    L : Flow length from the annulus inlet used in the entrance factor (1+(dh/L)^(2/3)).
        Use the full heat-exchanger length, NOT a segment length.
    q : Heat flux into fluid
    fluid : type of fluid being used

    Returns
    -------
    Tw : Wall temperature
    Tf : Film temperature
    htc : Heat Transfer Coefficient on inner wall
    Nu : Nusselt Number of inner wall
    f : Friction factor in the annulus

    """
    if not all(math.isfinite(v) for v in [Tb,Tw,P,a,G,dh,L,q]) or min(Tb,Tw,P,G,dh,L)<=0:
        raise ValueError('Gniel inputs must be finite, positive and in SI units (P in Pa)')
    if not 0<a<1:raise ValueError('Annulus diameter ratio must be in 0..1')
    P_Pa=P
    AS = _state(fluid)
    # Bulk properties do not change inside the loop (Tb fixed)
    AS.update(CP.PT_INPUTS, P_Pa, Tb)
    Bulk = {'rho': AS.rhomass(), 'mu': AS.viscosity(), 'K': AS.conductivity(), 'Pr': AS.Prandtl()}

    for i in range(100):
        previous_Tw=Tw
        Tf = (Tb+Tw)/2


        AS.update(CP.PT_INPUTS, P_Pa, Tw)
        wall = {'Pr': AS.Prandtl()}

        Re = G*dh/Bulk['mu']

        # Laminar below 2300; linear transition to turbulent at 4000.
        Re_eff = max(Re, 4000)
        Res = Re_eff*(((1+a**2)*math.log(a)+(1-a**2))/((1-a)**2*math.log(a)))
        fann = (1.8*math.log(Res,10)-1.5)**(-2)
        F = 0.75*a**(-0.17)
        K = (Bulk['Pr']/wall['Pr'])**(0.11)
        Nu = ((((fann/8)*(Re_eff-1000)*Bulk['Pr'])/(1+12.7*math.sqrt(fann/8)*(Bulk['Pr']**(2/3)-1)))*
              (1+(dh/L)**(2/3))*F*K)
        Nu_lam, poiseuille = laminar_annulus(float(a))
        if Re < 4000:
            weight = max(0.0, (Re-2300)/1700)
            Nu = (1-weight)*Nu_lam + weight*Nu
            fann = (1-weight)*poiseuille/Re + weight*fann
        htc = Nu*Bulk['K']/dh
        next_Tw=q/htc+Tb
        if not math.isfinite(next_Tw) or next_Tw<=0:raise ValueError("Unphysical water wall temperature")
        if abs(next_Tw-Tw)<1e-7:
            Tw=next_Tw;Tf=(Tb+Tw)/2;break
        Tw=.5*Tw+.5*next_Tw
    else:
        raise RuntimeError("Water correlation wall iteration did not converge")

    return pd.DataFrame([{'Tw':Tw,'Tf':Tf,'htc':htc,'Nu':Nu,'f':fann}])