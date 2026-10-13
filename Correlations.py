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


def DHC(Tb,Tw,P,D,G,q,fluid,count,ep,f):
    """
    Dang & Hihara corroelation

    Parameters
    ----------
    Tb : Bulk temperature of annular fluid
    Tw : Wall temperature of inner wall in annulus
    P : Pressure of annular fluid
    D : pipe Diameter
    G : Mass flux of fluid
    q : Heat flux into fluid
    fluid : type of fluid being used
    count : Position along pipe
    ep : Epsilon (relitive roughness (ep/Dh))

    Returns
    -------
    Tw : Wall temperature
    Tf : Film temperature
    htc : Heat Transfer Coefficient
    Nu : Nusselt Number
    f : Friction factor (Darcy)

    """
    if not 0 <= ep <= 1:
        raise ValueError("Effective e/D must be between 0 and 1")
    if ep > .05:
        _warn_extrapolation()
    # Clamp temperatures to reasonable bounds to prevent CoolProp errors
    # For CO2: minimum temperature around 216 K (triple point)
    T_min = 220.0  # Minimum temperature in K
    T_max = 1000.0  # Maximum temperature in K

    Tb = max(T_min, min(T_max, Tb))
    Tw = max(T_min, min(T_max, Tw))

    # Bulk Properties with error handling
    try:
        bulk = _props_DHC(fluid, Tb, P)
    except Exception as e:
        raise ValueError(f"Error calculating bulk properties at T={Tb:.2f}K, P={P:.2e}Pa: {e}")

    #itteration on inside wall temperature
    for l in range(100):
        Tw_previous = Tw
        try:
            wall = _props_DHC(fluid, Tw, P)
        except Exception as e:
            # If wall temperature is invalid, clamp it and retry
            Tw = max(T_min, min(T_max, Tw))
            wall = _props_DHC(fluid, Tw, P)

        Tf = (Tw + Tb) / 2
        Tf = max(T_min, min(T_max, Tf))  # Clamp film temperature too
        try:
            film = _props_DHC(fluid, Tf, P)
        except Exception as e:
            raise ValueError(f"Error calculating film properties at T={Tf:.2f}K, P={P:.2e}Pa: {e}")


        #Determination of CP_bar ensuring that the wall and bulk CP are not the same
        if abs(Tb - Tw) < 1e-6:
            Cpbar = bulk['Cp']
        else:
            Cpbar = (bulk['i'] - wall['i']) / (Tb - Tw)

        Ref = G*D/film['mu']
        Reb = G*D/bulk['mu']

        # determining the ratios of bulk and film viscosity and thermal conductivity
        MKb = bulk['mu']/bulk['k']
        MKf = film['mu']/film['k']

        #calculating the appropriate prandle number given by Dang & Hihara
        Pr_c = bulk['Cp']*bulk['mu']/bulk['k']     # default (also covers non-finite comparisons)
        if bulk['Cp'] >= Cpbar:
            Pr_c = bulk['Cp']*bulk['mu']/bulk['k']
        elif bulk['Cp'] < Cpbar and MKb >= MKf:
            Pr_c = Cpbar*bulk['mu']/bulk['k']
        elif bulk['Cp'] < Cpbar and MKb < MKf:
            Pr_c = Cpbar*film['mu']/film['k']


        #f = (1.82*math.log10(Ref)-1.64)**(-2)
        #incorporation of the white-Colbrooke friction factor to simulate
        #relative roughness on the CO2 side
        for i in range(100):
            previous_f = f
            f = (-2*math.log10(ep/3.7+2.51/(Ref*math.sqrt(f))))**(-2)
            if abs(f - previous_f) < 1e-12:
                break

        #determining the Nusselt #, heat transfer coefficent, and wall
        #temperature based on applied heat flux
        Nu = (f/8*(Reb-1000)*Pr_c)/(1.07+12.7*math.sqrt(f/8)*(Pr_c**(2/3)-1))
        htc = Nu*film['k']/D
        Tw = 0.5 * Tw_previous + 0.5 * (q/htc+Tb)

        # Clamp wall temperature to prevent invalid CoolProp calls
        Tw = max(T_min, min(T_max, Tw))
        if abs(Tw - Tw_previous) < 1e-7:
            break
    else:
        raise RuntimeError("Dang-Hihara wall iteration did not converge")

    return pd.DataFrame([{'Tw':Tw,'Tf':Tf,'htc':htc,'Nu':Nu,'f':f}])


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
    # Clamp Tb, Tw to valid range for REFPROP (avoid out-of-range CoolProp errors)
    T_lo, T_hi = 274.0, 2000.0  # K; water valid range
    Tb = max(T_lo, min(T_hi, float(Tb)))
    Tw = max(T_lo, min(T_hi, float(Tw)))
    # Convert pressure: check if in atm (small values) or already in Pa
    # If P < 100, assume it's in atm and convert to Pa
    if P < 100:
        P_Pa = P * 101325  # Convert atm to Pa
    else:
        P_Pa = P  # Already in Pa

    AS = _state(fluid)
    # Bulk properties do not change inside the loop (Tb fixed)
    AS.update(CP.PT_INPUTS, P_Pa, Tb)
    Bulk = {'rho': AS.rhomass(), 'mu': AS.viscosity(), 'K': AS.conductivity(), 'Pr': AS.Prandtl()}

    for i in range(0,3):
        Tf = (Tb+Tw)/2
        Tf = max(T_lo, min(T_hi, Tf))

        AS.update(CP.PT_INPUTS, P_Pa, Tw)
        wall = {'K': AS.conductivity(), 'Pr': AS.Prandtl()}
        AS.update(CP.PT_INPUTS, P_Pa, Tf)
        film = {'K': AS.conductivity()}

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
        Tw = q/htc+Tb
        Tw = max(T_lo, min(T_hi, Tw))

    return pd.DataFrame([{'Tw':Tw,'Tf':Tf,'htc':htc,'Nu':Nu,'f':fann}])