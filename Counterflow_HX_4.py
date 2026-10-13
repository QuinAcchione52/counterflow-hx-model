# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 2024

@author: Devon Hardy
Revised Oct 2026 (Q): true counterflow coupling.

Counter flow heat exchanger calculator for CO2 and water, exploring effects of roughness.
This is the main function to be called from CounterHX_full.py allowing iterations to be done for different inputs.
This code calls several other functions such as:
    Dang and Hihara for Nusselt number on surface of pipe for sCO2
    Gnielinski correlation for annular flows for heat transfer to inner surface of annulus for turbulent water
    REFPROP for fluid properties of both CO2 and Water through COOLPROP (AbstractState)

INDEXING CONVENTION (both fluids):
    Node k is at axial position x[k] for BOTH CO2 and water.
    CO2   enters at k = 0   and flows toward k = n-1  (+x).
    Water enters at k = n-1 and flows toward k = 0    (-x).
    Segment k lies between nodes k and k+1 (k = 0 .. n-2). Heat transfer in segment k is
    computed ONLY from CO2 and water temperatures at nodes k and k+1 (same physical location).

REVISION NOTES (Oct 2026):
    - Previous version paired CO2 node `count` with water node `n-1-count` (opposite ends of the
      exchanger). With water also marching n-1 -> 0, that coupling was parallel flow. Fixed.
    - Solution method (replaces the fixed 25 coupled sweeps):
        Outer loop (Picard): Dang-Hihara / Gnielinski at every node from the current profiles
          -> U(x), under-relaxed between iterations.
        Inner thermal solve for that U(x): forward march in +x from the CO2 inlet with a trial
          water OUTLET temperature at x = 0; each segment solved implicitly from the
          same-position LMTD; brentq on the water outlet temperature until the march hits the
          known water inlet temperature at x = L (shooting). Simple alternating CO2/water
          sweeps were tested and diverge at high NTU with balanced capacity rates.
        Iterates until bulk temperatures and total Q converge.
    - Constant = True (imposed q) still supported: CO2 marched 0 -> n-1, water n-1 -> 0.
    - LMTD fallback now uses the arithmetic mean (was half the DIFFERENCE, which gave Q = 0).
    - Water-side wall temperature includes the wall conduction drop (was set equal to the
      CO2-side wall temperature).
    - Segment quantities (Q, q, lmtd, Q_leak) exist for n-1 segments; the last output row is
      blank instead of a duplicated segment (old total Q double-counted one segment).
    - Linear heat leak is distributed over n-1 segments (was divided by n but applied to n-1).
    - Overall effectiveness = cumulative duty / Q_max (enthalpy-based). The old
      1 - prod(1 - eps_k) formula is not valid for counterflow segments in series.
    - Added CO2_dp_cumulative_psi column (output_improved.py plots it).
    - Internal storage uses numpy arrays instead of per-node pandas updates (much faster).
    - Gniel() now receives the full HX length L (was dx). The entrance factor (1+(dh/L)^(2/3))
      with L = dx multiplied the water Nu by ~4.2 at dx = 10 mm and changed with grid size.
    - Viscous pressure drop uses the segment-average BULK density (was film density).
"""
# import tool boxes
import pandas as pd
import numpy as np
import CoolProp.CoolProp as CP
from Correlations import DHC
from Correlations import Gniel, hydraulic_friction
from roughness import build_roughness
import math
import os
from scipy.optimize import brentq, root


#%% function call
def CO2_counter(m_dot_CO2, T_CO2, P_CO2, m_dot_H2O, T_H2O, P_H2O, D, wt, wg, L, Kp, ep, pro, p, Constant, q_CO2_i,
                whole_section=False, T_ambient=25., U_amb=0., d_outer_mm=0., Q_leak_slope=None, Q_leak_intercept=0.,
                max_iter=100, tol_T=1e-4, tol_Q=1e-6, relax=0.7, verbose=True, segments=None, water_htc_multiplier=1.0, target_mean_CO2_C=None, co2_friction_model='dh_roughness', co2_baseline_multiplier=1.0):
    """
    New optional arguments (defaults keep CounterHX_full.py working unchanged):
        max_iter : maximum number of outer (correlation) iterations
        tol_T    : convergence tolerance on max change in bulk temperature between iterations (K)
        tol_Q    : convergence tolerance on relative change in total Q between iterations (-)
        relax    : under-relaxation on U between iterations (1 = none; lower if it oscillates)
        target_mean_CO2_C : optional measured endpoint mean; infer CO2 inlet jointly
        water_htc_multiplier : diagnostic sensitivity only; default 1 preserves the model
    """

    #Fluid Phases
    PHASES={CP.iphase_liquid:'liquid',CP.iphase_gas:'gas',
            CP.iphase_twophase:'two_phase',CP.iphase_supercritical:'supercritical',
            CP.iphase_supercritical_liquid:'supercritical_liquid',
            CP.iphase_supercritical_gas:'supercritical_gas',
            CP.iphase_critical_point:'critical',CP.iphase_unknown:'unknown'}

    if not all(np.isfinite(v) for v in [m_dot_CO2,P_CO2,m_dot_H2O,T_H2O,P_H2O,D,wt,wg,L,Kp,q_CO2_i]):
        raise ValueError('Flows, pressure, geometry, water temperature and heat-flux input must be finite')
    if P_CO2<=0 or P_H2O<=0:
        raise ValueError('Both operating pressures must be positive')
    if not np.isfinite(co2_baseline_multiplier) or co2_baseline_multiplier<=0:
        raise ValueError('CO2 baseline multiplier must be positive and finite')
    if not np.isfinite(water_htc_multiplier) or water_htc_multiplier <= 0:
        raise ValueError("water_htc_multiplier must be positive and finite")
    if whole_section:
        raise ValueError("whole_section is an unvalidated placeholder; use the segmented solver")
    if min(m_dot_CO2, m_dot_H2O, D, wt, wg, L, Kp) <= 0:
        raise ValueError("Flows and geometry must be positive; D is the inner diameter in mm")
    if segments is not None and (segments < 1 or int(segments) != segments):
        raise ValueError("segments must be a positive integer")

    if target_mean_CO2_C is not None:
        if Constant or not np.isfinite(target_mean_CO2_C) or target_mean_CO2_C <= T_H2O:
            raise ValueError("Inverse bulk target must exceed water inlet and requires the coupled thermal solver")
        T_CO2 = float(target_mean_CO2_C) + .3  # initial guess only, not an experimental inlet

    if not np.isfinite(T_CO2):raise ValueError('CO2 inlet must be finite unless a valid inverse mean target is supplied')

    # inputs Pipe Geometry
    d_CO2 = D   # Inner diameter (mm) of CO2 pipe [e.g. 13.3858]
    t = wt       # Wall thickness (mm) of CO2 pipe [e.g. 1.2446]
    kp = Kp    # Thermal conductivity (W/m·K) of pipe [e.g. SS=14.4, Copper=386]
    wg = wg      # Water gap thickness (mm), radial [e.g. 28.4353 for gap to PVC ID]
    L = L     # Length of heat exchanger (mm) [e.g. 1000]

    # Inputs CO2
    F1 = 'CO2' #Hot fluid
    T1 = T_CO2  # Inlet temperature (°C) [e.g. 34.37]
    P1 = P_CO2   # Operating pressure (MPa) [e.g. 8.2]
    ASCO2 = CP.AbstractState(os.environ.get('HX_PROPERTY_BACKEND', 'HEOS'),F1)

    #Inputs Water
    F2 = "H2O" #Cooling fluid
    T2 = T_H2O  # Inlet temperature (°C) [e.g. 20.7]
    P2 = P_H2O   # Operating pressure (atm); ~1.05 atm for ambient cooling water
    ASH2O = CP.AbstractState(os.environ.get('HX_PROPERTY_BACKEND', 'HEOS'),F2)

    # other inputs

    fg = 0.02   #inital guess at friction factor
    pro_key = str(pro).strip().replace(' ', '')
    n = int(segments) + 1 if segments is not None else max(2, int(round(L/10)) + 1)      #number of points along pipe

    # conversions: all internal calculations in SI (m, m², Pa, K)
    L = L*1E-3                  # length (mm → m)
    t = t*1E-3                  # wall thickness (mm → m)
    wg = wg*1E-3                # water gap radial (mm → m)
    d_CO2 = d_CO2*1E-3          # CO2 pipe ID (mm → m)
    do = d_CO2+2*t              # CO2 pipe OD (m)
    d_H2O = do+2*wg             # water channel OD = ID of jacket (m)
    da = d_H2O-do               # hydraulic diameter of annulus (m)
    a = do/d_H2O                # diameter ratio of annulus (-)
    # Heat leakage: either linear in ΔT (real-time fit) or U_amb*A_outer*ΔT
    use_linear_Q_leak = Q_leak_slope is not None
    d_outer = d_outer_mm*1e-3 if d_outer_mm and (U_amb or use_linear_Q_leak) else 0.
    T_amb_K = (T_ambient + 273.15) if (use_linear_Q_leak or (U_amb and d_outer_mm)) else None
    A_CO2 = math.pi/4*d_CO2**2       # CO2 flow area (m²)
    A_H2O = math.pi/4*(d_H2O**2-do**2) # H2O flow area (m²)
    mdot_1 = m_dot_CO2               # mass flow rate CO2 (kg/s)
    mdot_2 = m_dot_H2O               # mass flow rate H2O (kg/s)
    G1 = mdot_1 / A_CO2              # mass flux CO2 (kg/m²·s)
    G2 = mdot_2 / A_H2O              # mass flux H2O (kg/m²·s)
    T1 = T1+273.15        # CO2 inlet (°C → K)
    T2 = T2+273.15        # H2O inlet (°C → K)
    P2 = P2*101325        # H2O pressure (atm → Pa)
    P1 = P1*1E6           # CO2 pressure (MPa → Pa)

    # ------------------------------------------------------------------
    # Discretization
    # ------------------------------------------------------------------
    x = np.linspace(0, L, n)
    dx = x[1] - x[0]
    nseg = n - 1                         # segment k spans nodes k, k+1
    As_CO2 = math.pi*d_CO2*dx            # surface area CO2 side of pipe per segment (inside, m^2)
    As_H2O = math.pi*do*dx               # surface area H2O side of pipe per segment (outside, m^2)
    A_outer = math.pi*d_outer*dx if d_outer and not use_linear_Q_leak else 0.   # for U_amb model only
    R2 = math.log(do / d_CO2) / (2 * math.pi * kp * dx)   # wall conduction resistance per segment (K/W)

    # ------------------------------------------------------------------
    # Node state storage (numpy). Suffix b = bulk, w = wall.
    # ------------------------------------------------------------------
    PROPS = [('rho', 'rhomass'), ('i', 'hmass'), ('mu', 'viscosity'), ('K', 'conductivity'),
             ('Cp', 'cpmass'), ('Pr', 'Prandtl'), ('a', 'speed_sound')]

    def new_state(P):
        S = {f'{nm}{s}': np.zeros(n) for nm, _ in PROPS for s in 'bw'}
        S['Tb'] = np.zeros(n)
        S['Tw'] = np.zeros(n)
        S['P'] = np.full(n, float(P))
        S['Phaseb'] = [''] * n
        S['Phasew'] = [''] * n
        return S

    def store(S, k, AS, s):
        """Copy the AbstractState's current state into node k (s = 'b' bulk or 'w' wall)."""
        S[f'T{s}'][k] = AS.T()
        for nm, meth in PROPS:
            S[f'{nm}{s}'][k] = getattr(AS, meth)()
        S[f'Phase{s}'][k] = PHASES.get(AS.phase(), 'unknown')

    CO2 = new_state(P1)
    H2O = new_state(P2)
    ASCO2.update(CP.PT_INPUTS, P1, T1)
    for k in range(n):
        store(CO2, k, ASCO2, 'b'); store(CO2, k, ASCO2, 'w')
    ASH2O.update(CP.PT_INPUTS, P2, T2)
    for k in range(n):
        store(H2O, k, ASH2O, 'b'); store(H2O, k, ASH2O, 'w')

    rough=build_roughness(x,ep,pro,p)
    if co2_friction_model=='dang_hihara' and np.any(rough!=0):
        raise ValueError('Original smooth Dang-Hihara mode requires an all-zero profile')

    # ------------------------------------------------------------------
    # Heat transfer storage
    #   Segment arrays (length n, entry n-1 unused): Q, Q_leak, lmtd, dT
    #   Node arrays: Nu, htc, R1, R3, U, f, Tf
    #   Sign: Q < 0 means heat leaves CO2 (CO2 is cooled).
    # ------------------------------------------------------------------
    Q = np.zeros(n)
    if Constant:
        Q[:nseg] = q_CO2_i*As_CO2  # prescribed physical flux
    # Coupled runs start at zero flux: the arbitrary seed must not create a
    # subfreezing water-wall state before the first heat-balance solve.
    Q_leak = np.zeros(n)
    lmtd = np.zeros(n);   dT = np.zeros(n)
    Nu_CO2 = np.zeros(n); htc_CO2 = np.zeros(n)
    Nu_H2O = np.zeros(n); htc_H2O = np.zeros(n)
    R1 = np.zeros(n);     R3 = np.zeros(n);  U = np.zeros(n)
    f_CO2 = np.zeros(n);  f_H2O = np.zeros(n)
    Tf_CO2 = np.zeros(n); Tf_H2O = np.zeros(n)
    diagnostic={key:np.zeros(n) for key in ["Re_film","Cp_bulk","Cp_integrated","Pr_used","k_film","wall_flux_residual_W_m2"]}
    pr_branch=[""]*n

    # Linear heat leakage: one total from your fit (based on inlet water T), spread over the segments
    if use_linear_Q_leak:
        dT_C = T_H2O - T_ambient                            # inlet H2O (°C) - T_amb (°C)
        y_leak = Q_leak_slope * dT_C + Q_leak_intercept     # fit: heat into water (W)
        total_heat_to_ambient = -y_leak                     # positive = heat leaving water to ambient
        Q_leak_per_seg = total_heat_to_ambient / nseg

    def node_avg(seg):
        """Segment values (length nseg) -> node values (length n): mean of adjacent segments."""
        out = np.empty(n)
        out[0] = seg[0]
        out[-1] = seg[-1]
        out[1:-1] = 0.5 * (seg[:-1] + seg[1:])
        return out

    def lmtd_fn(dTa, dTb):
        """Segment LMTD from the two end temperature differences (hot - cold at the same x)."""
        if dTa > 0 and dTb > 0:
            if abs(dTa - dTb) <= 1e-6 * max(dTa, dTb):
                return 0.5 * (dTa + dTb)
            return (dTa - dTb) / math.log(dTa / dTb)
        # temperature cross inside the segment (should not persist at convergence)
        return 0.5 * (dTa + dTb)

    h_H2O_in = H2O['ib'][n - 1]     # water inlet enthalpy (fixed boundary condition at x = L)
    T_lo, T_hi = min(T1, T2) - 60.0, max(T1, T2) + 60.0   # runaway guard for trial marches

    def leak_seg(k, Tw_a, Tw_b):
        """Heat leaving the water to ambient in segment k (W)."""
        if use_linear_Q_leak:
            return Q_leak_per_seg
        if A_outer and T_amb_K is not None:
            return U_amb * A_outer * (0.5 * (Tw_a + Tw_b) - T_amb_K)
        return 0.0

    Q_seed = np.zeros(n)            # last-march segment Q, used to start each local solve

    trial_CO2_out = np.nan

    def march(h_w_out, store_props):
        """
        Forward march in +x for a trial water OUTLET enthalpy at x = 0 (frozen U profile).
        CO2 enthalpy known at k=0 (inlet); water enthalpy at k=0 is the trial outlet value.
        Segment k is solved implicitly: Q_k = -U_seg*As*LMTD(dT_k, dT_k+1), with both
        end temperatures taken at the same axial positions (counterflow).
        Water energy balance: h_w[k] = h_w[k+1] + (-Q_k - Q_leak_k)/mdot_2
                          ->  h_w[k+1] = h_w[k] + (Q_k + Q_leak_k)/mdot_2
        Returns residual h_w[n-1] - h_w_in (zero when the water inlet condition is met).
        """
        nonlocal trial_CO2_out
        trial_CO2_out = np.nan
        hc = CO2['ib'][0]
        hw = h_w_out
        Tc = CO2['Tb'][0]
        try:
            ASH2O.update(CP.HmassP_INPUTS, hw, H2O['P'][0])
        except Exception:
            return -1e12 if h_w_out < h_H2O_in else 1e12
        Tw = ASH2O.T()
        if store_props:
            store(H2O, 0, ASH2O, 'b')
        for k in range(nseg):
            U_seg = 0.5 * (U[k] + U[k + 1])
            # start from this segment's Q in the previous march (explicit estimate the first time)
            Qk = Q_seed[k] if Q_seed[k] != 0.0 else -U_seg * As_CO2 * (Tc - Tw)
            Ql = leak_seg(k, Tw, Tw)
            Tc1, Tw1 = Tc, Tw
            try:
                for _ in range(30):                   # local fixed point (segment NTU << 1)
                    hc1 = hc + Qk / mdot_1
                    ASCO2.update(CP.HmassP_INPUTS, hc1, CO2['P'][k + 1]); Tc1 = ASCO2.T()
                    Ql = leak_seg(k, Tw, Tw1)
                    hw1 = hw + (Qk + Ql) / mdot_2
                    ASH2O.update(CP.HmassP_INPUTS, hw1, H2O['P'][k + 1]); Tw1 = ASH2O.T()
                    lm = lmtd_fn(Tc - Tw, Tc1 - Tw1)
                    Q_new = -U_seg * As_CO2 * lm
                    if abs(Q_new - Qk) <= 1e-8 * max(1.0, abs(Q_new)):
                        break                         # property states already match Qk
                    Qk = Q_new
            except Exception:
                if store_props:
                    raise
                return hw - h_H2O_in                  # trial guess ran out of property range
            if not (T_lo < Tc1 < T_hi and T_lo < Tw1 < T_hi):
                if store_props:
                    raise RuntimeError("Final march left the valid temperature bracket")
                return hw1 - h_H2O_in                 # runaway trial; sign still points the right way
            if store_props:
                Q[k] = Qk
                Q_leak[k] = Ql
                lmtd[k] = lmtd_fn(Tc - Tw, Tc1 - Tw1)
                dT[k] = Tc - Tw
                store(CO2, k + 1, ASCO2, 'b')
                store(H2O, k + 1, ASH2O, 'b')
            Q_seed[k] = Qk
            hc, hw, Tc, Tw = hc1, hw1, Tc1, Tw1
        trial_CO2_out = Tc
        return hw - h_H2O_in

    #% Outer Picard loop: correlations <-> counterflow thermal solve
    converged = False
    dT_change = np.nan
    dQ_rel = np.nan
    U_prev = None
    Tw_out_prev = None
    for j in range(max_iter):
        wall_prev = CO2['Tw'].copy()
        Tc_prev = CO2['Tb'].copy()
        Th2o_prev = H2O['Tb'].copy()
        Qtot_prev = Q[:nseg].sum()

        # node heat rate per segment-length, from current segment Q
        Qn = (node_avg(Q[:nseg]) if Constant or U_prev is None
              else -U_prev*As_CO2*(CO2['Tb']-H2O['Tb']))

        # ---- 1. Correlations at every node (CO2 and water at the SAME position k) ----
        for k in range(n):
            ## CO2 side: Dang & Hihara
            DH = DHC(CO2['Tb'][k], CO2['Tw'][k], CO2['P'][k], d_CO2, G1,
                     Qn[k] / As_CO2, F1, k, float(rough[k]), fg, friction_model=co2_friction_model, baseline_multiplier=co2_baseline_multiplier)
            Nu_CO2[k], htc_CO2[k] = DH.Nu[0], DH.htc[0]
            f_CO2[k], Tf_CO2[k] = DH.f[0], DH.Tf[0]
            Tw_inner = DH.Tw[0]
            for key in diagnostic:diagnostic[key][k]=DH[key].iloc[0]
            pr_branch[k]=DH.Pr_branch.iloc[0]
            ASCO2.update(CP.PT_INPUTS, CO2['P'][k], Tw_inner)
            store(CO2, k, ASCO2, 'w')

            ## Water side: Gnielinski (annulus); heat flux into water is -Q/As_H2O
            # Pass the FULL exchanger length: Gnielinski's (1 + (dh/L)^(2/3)) entrance factor is a
            # length-averaged correction from the flow inlet. Passing dx (old code) inflated the
            # water htc ~4x at dx = 10 mm and made it grid-dependent.
            Gne = Gniel(H2O['Tb'][k], H2O['Tw'][k], H2O['P'][k], a, G2, da, L,
                        -Qn[k] / As_H2O / water_htc_multiplier, F2)
            Nu_H2O[k], htc_H2O[k] = water_htc_multiplier*Gne.Nu[0], water_htc_multiplier*Gne.htc[0]
            f_H2O[k], Tf_H2O[k] = Gne.f[0], Gne.Tf[0]
            # Outer wall = inner wall minus conduction drop across the tube wall.
            # Heat leaving CO2 through the wall = -Qn  ->  Tw_outer = Tw_inner + Qn*R2
            Tw_outer = Tw_inner + Qn[k] * R2
            ASH2O.update(CP.PT_INPUTS, H2O['P'][k], Tw_outer)
            store(H2O, k, ASH2O, 'w')
            # thermal resistances and U at node k (both sides at the same position)
            h1 = max(htc_CO2[k], 1.0)
            h3 = max(htc_H2O[k], 1.0)
            R1[k] = 1 / (h1 * As_CO2)
            R3[k] = 1 / (h3 * As_H2O)
            U[k] = 1 / As_CO2 * (1 / (R1[k] + R2 + R3[k]))

        # under-relax U between Picard iterations (htc depends strongly on q and Tw near the
        # pseudocritical point)
        if U_prev is not None:
            U[:] = (1 - relax) * U_prev + relax * U
        U_prev = U.copy()

        # ---- 2. Thermal solve for this U profile ----
        if Constant:
            # Imposed heat flux: CO2 marched 0 -> n-1, water marched n-1 -> 0 (inlet known)
            for k in range(nseg):
                Q_leak[k] = leak_seg(k, H2O['Tb'][k], H2O['Tb'][k + 1])
                ASCO2.update(CP.HmassP_INPUTS, CO2['ib'][k] + Q[k] / mdot_1, CO2['P'][k + 1])
                store(CO2, k + 1, ASCO2, 'b')
            for k in range(nseg - 1, -1, -1):
                ASH2O.update(CP.HmassP_INPUTS, H2O['ib'][k + 1] + (-Q[k] - Q_leak[k]) / mdot_2,
                             H2O['P'][k])
                store(H2O, k, ASH2O, 'b')
            for k in range(nseg):
                lmtd[k] = lmtd_fn(CO2['Tb'][k] - H2O['Tb'][k], CO2['Tb'][k + 1] - H2O['Tb'][k + 1])
                dT[k] = CO2['Tb'][k] - H2O['Tb'][k]
        elif target_mean_CO2_C is not None:
            # Solve the two inlet/outlet constraints together for frozen U:
            # unknown CO2 inlet, unknown water outlet; mean CO2 and water inlet fixed.
            target_K = float(target_mean_CO2_C)+273.15
            def inverse_residual(temperatures):
                hot_in, water_out = map(float, temperatures)
                if not (target_K-10 < hot_in < target_K+30 and 274. < water_out < 449.):
                    return [1e3+abs(hot_in-target_K), 1e3+abs(water_out-T2)]
                ASCO2.update(CP.PT_INPUTS, P1, hot_in)
                store(CO2, 0, ASCO2, 'b')
                ASH2O.update(CP.PT_INPUTS, P2, water_out)
                water_residual = march(ASH2O.hmass(), False)
                if not np.isfinite(trial_CO2_out):
                    raise RuntimeError("Inverse trial march did not complete")
                return [water_residual/4200., (hot_in+trial_CO2_out)/2-target_K]
            guess=[T1, Tw_out_prev if Tw_out_prev is not None else T2+.2]
            solution=root(inverse_residual,guess,tol=1e-9)
            final_residual=np.asarray(inverse_residual(solution.x))
            if not np.all(np.isfinite(final_residual)) or np.max(np.abs(final_residual))>1e-6:
                raise RuntimeError(f"Joint inlet/bulk solve did not converge: {final_residual}")
            T1,Tw_out_prev=map(float,solution.x)
            ASCO2.update(CP.PT_INPUTS,P1,T1)
            store(CO2,0,ASCO2,'b')
            ASH2O.update(CP.PT_INPUTS,P2,Tw_out_prev)
            inlet_residual=march(ASH2O.hmass(),True)
            if abs(inlet_residual)>1e-3:
                raise RuntimeError("Inverse solution missed water inlet boundary")
        else:
            # Shooting on the water outlet temperature (x = 0) so the forward march hits the
            # known water inlet temperature at x = L.
            def resid_T(Tw_out):
                ASH2O.update(CP.PT_INPUTS, H2O['P'][0], Tw_out)
                return march(ASH2O.hmass(), False)
            T_min = max(T2 - 10.0, 274.5)
            T_max = min(max(T1, T2) + 10.0, 449.0)
            if Tw_out_prev is None:
                T_a, T_b = T_min, T_max
                f_a, f_b = resid_T(T_a), resid_T(T_b)
            else:
                # warm start: the root moves very little between outer iterations
                half = 0.05
                while True:
                    T_a, T_b = max(Tw_out_prev - half, T_min), min(Tw_out_prev + half, T_max)
                    f_a, f_b = resid_T(T_a), resid_T(T_b)
                    if f_a * f_b <= 0 or (T_a == T_min and T_b == T_max):
                        break
                    half *= 5
            if f_a * f_b > 0:
                raise RuntimeError(
                    f"Could not bracket water outlet temperature between {T_a-273.15:.1f} and "
                    f"{T_b-273.15:.1f} °C (residuals {f_a:.3g}, {f_b:.3g} J/kg).")
            Tw_out = brentq(resid_T, T_a, T_b, xtol=1e-6, rtol=1e-12, maxiter=200)
            Tw_out_prev = Tw_out
            ASH2O.update(CP.PT_INPUTS, H2O['P'][0], Tw_out)
            inlet_residual = march(ASH2O.hmass(), True)  # final march stores the solution profiles
            if abs(inlet_residual) > 1e-3:
                raise RuntimeError("Final water inlet enthalpy residual exceeds 0.001 J/kg")

        # ---- 3. Convergence of the Picard loop ----
        dT_change = max(np.max(np.abs(CO2['Tb'] - Tc_prev)), np.max(np.abs(H2O['Tb'] - Th2o_prev)))
        Qtot = Q[:nseg].sum()
        dQ_rel = abs(Qtot - Qtot_prev) / max(abs(Qtot), 1e-12)
        local_node_Q=node_avg(Q[:nseg]) if Constant else -U*As_CO2*(CO2['Tb']-H2O['Tb'])
        water_wall_residual=htc_H2O*(H2O['Tw']-H2O['Tb'])+local_node_Q/As_H2O
        wall_balance_relative=float(np.max(np.abs(water_wall_residual)/np.maximum(1.,np.abs(local_node_Q/As_H2O))))
        if j >= 1 and dT_change < tol_T and dQ_rel < tol_Q and np.max(np.abs(CO2['Tw']-wall_prev)) < tol_T and (Constant or wall_balance_relative<1e-5):
            converged = True
            break

    n_sweeps = j + 1
    if verbose:
        status = "converged" if converged else "NOT converged"
        print(f"Counterflow solver: {n_sweeps} outer iterations ({status}; max dT change {dT_change:.2e} K, "
              f"rel dQ {dQ_rel:.2e})")
    if not converged:
        raise RuntimeError("Solver did not converge; lower relax or increase max_iter")

    # ------------------------------------------------------------------
    # Pressure drop and pumping power (on the converged state)
    # CO2 moves k-1 -> k ; water moves k+1 -> k
    # ------------------------------------------------------------------
    u_CO2 = G1 / CO2['rhob']
    u_H2O = G2 / H2O['rhob']
    Re_hydraulic=G1*d_CO2/CO2['mub']
    f_hydraulic=np.array([hydraulic_friction(float(re),float(e)) for re,e in zip(Re_hydraulic,rough)])
    CO2_Visc = np.zeros(n); CO2_acc = np.zeros(n)
    H2O_Visc = np.zeros(n); H2O_acc = np.zeros(n)
    for k in range(n):
        if k >= 1:
            rho_avg = 0.5 * (CO2['rhob'][k - 1] + CO2['rhob'][k])     # bulk density over segment
            CO2_Visc[k] = dx/d_CO2*G1**2/4*(f_hydraulic[k-1]/CO2['rhob'][k-1]+f_hydraulic[k]/CO2['rhob'][k])
            CO2_acc[k] = G1**2 * (1 / CO2['rhob'][k] - 1 / CO2['rhob'][k - 1])
        if k <= n - 2:
            rho_avg = 0.5 * (H2O['rhob'][k] + H2O['rhob'][k + 1])
            H2O_Visc[k] = f_H2O[k] * dx / da * G2**2 / (2 * rho_avg)
            H2O_acc[k] = G2**2 * (1 / H2O['rhob'][k] - 1 / H2O['rhob'][k + 1])
    CO2_Pp = (CO2_Visc + CO2_acc) * u_CO2 * A_CO2
    H2O_Pp = (H2O_Visc + H2O_acc) * u_H2O * A_H2O

    Re_CO2b = G1 * d_CO2 / CO2['mub']
    Re_H2Ob = G2 * da / H2O['mub']

    #%% plots and outputs
    def seg_col(arr):
        """Segment array for output: n-1 segment values, blank last row."""
        out = arr.astype(float).copy()
        out[-1] = np.nan
        return out

    node_Q_final=node_avg(Q[:nseg]) if Constant else -U*As_CO2*(CO2['Tb']-H2O['Tb'])
    q_CO2_seg = Q / As_CO2
    q_H2O_seg = -Q / As_H2O

    plot = pd.DataFrame({
        'pos'       : x,
        'roughness' : np.ravel(rough),
        'Re_CO2'    : Re_CO2b,
        'Re_H2O'    : Re_H2Ob,
        'CO2b'      : CO2['Tb'],
        'CO2w'      : CO2['Tw'],
        'H2Ob'      : H2O['Tb'],
        'H2Ow'      : H2O['Tw'],
        'H2Ob_at_pos' : H2O['Tb'].copy(),
        'dT_bulk'   : CO2['Tb'] - H2O['Tb'],
        'Q'         : seg_col(Q),
        'q_CO2'     : seg_col(q_CO2_seg),
        'q_H2O'     : seg_col(q_H2O_seg),
        'Q_leak'    : seg_col(Q_leak),
        'lmtd'      : seg_col(lmtd),
        'CO2_heat_friction_factor':f_CO2,
        'CO2_hydraulic_friction_factor':f_hydraulic,
        'CO2_phase_bulk':CO2['Phaseb'],
        'CO2_Pr_branch':pr_branch,
        'CO2_Re_film':diagnostic['Re_film'],
        'CO2_Cp_bulk_J_kgK':diagnostic['Cp_bulk'],
        'CO2_Cp_integrated_J_kgK':diagnostic['Cp_integrated'],
        'CO2_Pr_used':diagnostic['Pr_used'],
        'CO2_k_film_W_mK':diagnostic['k_film'],
        'CO2_node_heat_flux_W_m2':node_Q_final/As_CO2,
        'H2O_node_heat_flux_W_m2':-node_Q_final/As_H2O,
        'CO2_correlation_wall_flux_residual_W_m2':diagnostic['wall_flux_residual_W_m2'],
        'CO2_final_wall_flux_residual_W_m2':htc_CO2*(CO2['Tw']-CO2['Tb'])-node_Q_final/As_CO2,
        'H2O_final_wall_flux_residual_W_m2':htc_H2O*(H2O['Tw']-H2O['Tb'])+node_Q_final/As_H2O,
        'wall_conduction_residual_K':CO2['Tw']-H2O['Tw']+node_Q_final*R2,
        'Nu_CO2'    : Nu_CO2,
        'htc_CO2'   : htc_CO2,
        'Nu_H2O'    : Nu_H2O,
        'htc_H2O'   : htc_H2O,
        'CO2_mCp'   : mdot_1 * CO2['Cpb'],
        'H2O_mCp'   : mdot_2 * H2O['Cpb'],
        'r_CO2'     : R1,
        'r_cond'    : np.full(n, R2),
        'r_H2O'     : R3,
        'U'         : U,
        'CO2_visc'  : CO2_Visc,
        'CO2_acc'   : CO2_acc,
        'CO2_tot'   : CO2_Visc + CO2_acc,
        'Ve_CO2'    : u_CO2,
        'H2O_visc'  : H2O_Visc,
        'H2O_acc'   : H2O_acc,
        'H2O_tot'   : H2O_Visc + H2O_acc,
        'Ve_H2O'    : u_H2O,
        'CO2_dp_cumulative_psi': np.cumsum(CO2_Visc + CO2_acc) / 6894.757,
    })

    # Local (segment) effectiveness from node NTU, counterflow formula (C_r = 1 handled)
    Cmin = plot[['CO2_mCp', 'H2O_mCp']].min(axis=1)
    Cmax = plot[['CO2_mCp', 'H2O_mCp']].max(axis=1)
    plot['C'] = Cmin / Cmax
    plot['NTU'] = (plot['U'] * As_CO2) / Cmin
    Cr = plot['C'].values
    NTU = plot['NTU'].values
    with np.errstate(divide='ignore', invalid='ignore'):
        eps_seg = np.where(np.abs(1 - Cr) < 1e-9,
                           NTU / (1 + NTU),
                           (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(-NTU * (1 - Cr))))
    plot['effectiveness'] = eps_seg * 100

    # Overall effectiveness: cumulative duty / Q_max (enthalpy-based, valid with variable cp)
    ASCO2.update(CP.PT_INPUTS, P1, T1); h_CO2_in = ASCO2.hmass()
    ASCO2.update(CP.PT_INPUTS, P1, T2); h_CO2_at_Tc = ASCO2.hmass()
    ASH2O.update(CP.PT_INPUTS, P2, T2); h_H2O_in = ASH2O.hmass()
    ASH2O.update(CP.PT_INPUTS, P2, max(274.0, min(450.0, T1))); h_H2O_at_Th = ASH2O.hmass()
    Q_max = min(mdot_1 * (h_CO2_in - h_CO2_at_Tc), mdot_2 * (h_H2O_at_Th - h_H2O_in))
    cum_duty = np.concatenate([[0.0], np.cumsum(-Q[:nseg])])
    plot['overall effectiveness'] = cum_duty / Q_max * 100

    # ------------------------------------------------------------------
    # Printed summary
    # ------------------------------------------------------------------
    T_CO2_out_K = CO2['Tb'][n - 1]
    T_H2O_out_K = H2O['Tb'][0]
    total_Q = Q[:nseg].sum()                         # W, negative = heat out of CO2
    total_Q_leak = Q_leak[:nseg].sum()               # W, positive = water -> ambient
    duty_CO2 = mdot_1 * (CO2['ib'][0] - CO2['ib'][n - 1])
    duty_H2O = mdot_2 * (H2O['ib'][0] - H2O['ib'][n - 1])
    avg_q_CO2 = np.mean(q_CO2_seg[:nseg])
    avg_q_H2O = np.mean(q_H2O_seg[:nseg])
    lmtd_valid = lmtd[:nseg][lmtd[:nseg] != 0]
    avg_lmtd = float(lmtd_valid.mean()) if len(lmtd_valid) else 0.0

    if verbose:
        print(f"Inlet temperature: CO2 {T1-273.15:.2f} °C, H2O {T2-273.15:.2f} °C")
        print(f"Outlet temperature: CO2 {T_CO2_out_K-273.15:.2f} °C, H2O {T_H2O_out_K-273.15:.2f} °C")
        print(f"Delta T: CO2 {T1 - T_CO2_out_K:.2f} K, H2O {T_H2O_out_K - T2:.2f} K")
        print(f"Inlet Reynolds: CO2 {Re_CO2b[0]:.0f}, H2O {Re_H2Ob[n-1]:.0f}")
        print(f"Average heat flux: CO2 {avg_q_CO2:.0f} W/m², H2O {avg_q_H2O:.0f} W/m²")
        print(f"Average HTC: CO2 {np.mean(htc_CO2):.0f} W/(m²·K), H2O {np.mean(htc_H2O):.0f} W/(m²·K)")
        print(f"Total heat transfer Q: {total_Q:.2f} W  |  Avg LMTD: {avg_lmtd:.2f} K  |  "
              f"Avg U: {np.mean(U):.1f} W/(m²·K)  |  Effectiveness: {cum_duty[-1]/Q_max*100:.1f} %")
        # Energy balance: CO2 duty = |sum Q| ; water gain = |sum Q| - leakage
        imb_CO2 = duty_CO2 - (-total_Q)
        imb_H2O = duty_H2O - (-total_Q - total_Q_leak)
        print(f"Energy balance: CO2 duty {duty_CO2:.3f} W, water gain {duty_H2O:.3f} W, "
              f"leakage to ambient {total_Q_leak:.3f} W  "
              f"(residuals {imb_CO2:.2e} W, {imb_H2O:.2e} W)")
        if total_Q_leak != 0:
            if use_linear_Q_leak:
                print(f"  (linear leak: slope={Q_leak_slope} W/°C, intercept={Q_leak_intercept} W, "
                      f"heat into water from ambient = {-total_Q_leak:.2f} W)")
            elif A_outer and U_amb:
                print(f"  (T_amb={T_ambient} °C, U_amb={U_amb} W/m²·K)")
        if np.any(CO2['Tb'] - H2O['Tb'] <= 0):
            print("NOTE: CO2 bulk temperature falls to/below water temperature somewhere along the HX.")

    from hx_metrics import summarize
    plot.attrs['summary'] = summarize(plot, D, wt)
    plot.attrs['summary'].update({
        'property_backend': os.environ.get('HX_PROPERTY_BACKEND', 'HEOS'),
        'iterations': n_sweeps, 'converged': converged,
        'water_htc_multiplier': float(water_htc_multiplier),
        'CO2_heat_friction_model':co2_friction_model,
        'CO2_baseline_correction_multiplier':float(co2_baseline_multiplier),
        'CO2_baseline_is_calibrated':bool(co2_baseline_multiplier!=1),
        'CO2_pressure_friction_model':'bulk_Re_Colebrook',
        'CO2_mass_flux_kg_m2s':G1,
        'water_wall_balance_max_relative_residual':wall_balance_relative,
        'CO2_diameter_outside_cited_DH_database':not (0.001<=d_CO2<=0.006),
        'CO2_mass_flux_below_cited_200_kg_m2s':G1<200,
        'water_htc_mean_W_m2K': float(np.mean(htc_H2O)),
        'CO2_energy_residual_W': float(duty_CO2 + total_Q),
        'water_energy_residual_W': float(duty_H2O + total_Q + total_Q_leak),
    })
    if target_mean_CO2_C is not None:
        if abs(plot.attrs['summary']['CO2_mean_endpoint_C']-target_mean_CO2_C)>1e-4:
            raise RuntimeError("Inverse solution missed mean bulk temperature")
        plot.attrs['summary'].update(inferred_CO2_inlet_C=T1-273.15,
                                     target_bulk_C=float(target_mean_CO2_C),
                                     CO2_inlet_is_inferred=True)
    if verbose:
        print(f"Spencer-definition HTC: {plot.attrs['summary']['htc_spencer_W_m2K']:.1f} W/(m²·K)")
    #%% return results
    return plot