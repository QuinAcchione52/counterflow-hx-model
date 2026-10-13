# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 08:16:31 2024

@author: dhard

Revised Oct 2026 (Q):
    - Property calls use cached CoolProp AbstractState objects (REFPROP backend) instead of
      CP.PropsSI. Same REFPROP values, far less overhead (PropsSI re-resolves the fluid and
      backend on every call; DHC made ~55 and Gniel ~21 such calls per node per iteration).
    - DHC: Pr_c initialized to the bulk Prandtl number before the branch logic so it can never be
      undefined.
    - Equations are otherwise unchanged. Modeling choices left as-is pending review:
        * Colebrook friction factor evaluated with the FILM Reynolds number (Ref)
        * count == 0 forces bulk Cp at the inlet node
        * Gniel clamps Re to >= 2300 (turbulent correlation applied to laminar water flow)
    - NOTE for callers: Gniel's L argument enters Gnielinski's length-averaged entrance factor
      (1 + (dh/L)^(2/3)). It must be the flow length from the inlet (e.g. the full HX length),
      not a discretization segment length.
"""

import pandas as pd
import CoolProp.CoolProp as CP
import math

# One AbstractState per fluid, reused across calls
_AS = {}


def _state(fluid):
    AS = _AS.get(fluid)
    if AS is None:
        AS = CP.AbstractState('REFPROP', fluid)
        _AS[fluid] = AS
    return AS


def _props_DHC(fluid, T, P):
    """rho, mu, Cp, k, i at (T, P) — same set the original PropsSI calls returned."""
    AS = _state(fluid)
    AS.update(CP.PT_INPUTS, P, T)
    return {'rho': AS.rhomass(), 'mu': AS.viscosity(), 'Cp': AS.cpmass(),
            'k': AS.conductivity(), 'i': AS.hmass()}


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
    for l in range(0,5):
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
        if count == 0 or abs(Tb - Tw) < 1e-6:
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
        if bulk['Cp'] >= Cpbar or count == 0:
            Pr_c = bulk['Cp']*bulk['mu']/bulk['k']
        elif bulk['Cp'] < Cpbar and MKb >= MKf:
            Pr_c = Cpbar*bulk['mu']/bulk['k']
        elif bulk['Cp'] < Cpbar and MKb < MKf:
            Pr_c = Cpbar*film['mu']/film['k']


        #f = (1.82*math.log10(Ref)-1.64)**(-2)
        #incorporation of the white-Colbrooke friction factor to simulate
        #relative roughness on the CO2 side
        for i in range(0,5):
            f = (-2*math.log10(ep/3.7+2.51/(Ref*math.sqrt(f))))**(-2)

        #determining the Nusselt #, heat transfer coefficent, and wall
        #temperature based on applied heat flux
        Nu = (f/8*(Reb-1000)*Pr_c)/(1.07+12.7*math.sqrt(f/8)*(Pr_c**(2/3)-1))
        htc = Nu*film['k']/D
        Tw = q/htc+Tb

        # Clamp wall temperature to prevent invalid CoolProp calls
        Tw = max(T_min, min(T_max, Tw))

    return pd.DataFrame([{'Tw':Tw,'Tf':Tf,'htc':htc,'Nu':Nu,'f':f}])

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

        # H2O annulus: always use turbulent correlation (no laminar branch) to match rig behavior.
        Re_eff = max(Re, 2300)
        Res = Re_eff*(((1+a**2)*math.log(a)+(1-a**2))/((1-a)**2*math.log(a)))
        fann = (1.8*math.log(Res,10)-1.5)**(-2)
        F = 0.75*a**(-0.17)
        K = (Bulk['Pr']/wall['Pr'])**(0.11)
        Nu = ((((fann/8)*(Re_eff-1000)*Bulk['Pr'])/(1+12.7*math.sqrt(fann/8)*(Bulk['Pr']**(2/3)-1)))*
              (1+(dh/L)**(2/3))*F*K)
        htc = Nu*film['K']/dh
        Tw = q/htc+Tb
        Tw = max(T_lo, min(T_hi, Tw))

    return pd.DataFrame([{'Tw':Tw,'Tf':Tf,'htc':htc,'Nu':Nu,'f':fann}])