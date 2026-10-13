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
from Correlations import Gniel
import math
from scipy.optimize import brentq


#%% function call
def CO2_counter(m_dot_CO2, T_CO2, P_CO2, m_dot_H2O, T_H2O, P_H2O, D, wt, wg, L, Kp, ep, pro, p, Constant, q_CO2_i,
                whole_section=False, T_ambient=25., U_amb=0., d_outer_mm=0., Q_leak_slope=None, Q_leak_intercept=0.,
                max_iter=100, tol_T=1e-4, tol_Q=1e-6, relax=0.7, verbose=True):
    """
    New optional arguments (defaults keep CounterHX_full.py working unchanged):
        max_iter : maximum number of outer (correlation) iterations
        tol_T    : convergence tolerance on max change in bulk temperature between iterations (K)
        tol_Q    : convergence tolerance on relative change in total Q between iterations (-)
        relax    : under-relaxation on U between iterations (1 = none; lower if it oscillates)
    """

    #Fluid Phases
    PHASES = {
        0: "invalid",
        1: "liquid",
        2: "two-phase",
        3: "gas",
        4: "supercritical",
        5: "supercritical_liquid",
        6: "supercritical_gas",
        7: "critical"
    }

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
    ASCO2 = CP.AbstractState('REFPROP',F1)

    #Inputs Water
    F2 = "H2O" #Cooling fluid
    T2 = T_H2O  # Inlet temperature (°C) [e.g. 20.7]
    P2 = P_H2O   # Operating pressure (atm); ~1.05 atm for ambient cooling water
    ASH2O = CP.AbstractState('REFPROP',F2)

    # other inputs
    ep_max_input = float(ep) if np.isscalar(ep) else float(np.ravel(ep)[0])
    fg = 0.02   #inital guess at friction factor
    pro_key = str(pro).strip().replace(' ', '')
    n = round(int(L/10+1))      #number of points along pipe

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

    # --- Whole test section (one LMTD, one Q for full length) ---
    # NOTE: lumped placeholder; effectiveness and wall temperatures here are not physical.
    if whole_section:
        As_CO2_total = math.pi * d_CO2 * L
        As_H2O_total = math.pi * do * L
        R2_total = math.log(do / d_CO2) / (2 * math.pi * kp * L)
        ep_val = float(ep) if np.isscalar(ep) else float(ep[0])
        T_CO2_out_K = T1 - 5.0
        T_H2O_out_K = T2 + 5.0
        Q_total = 0.0
        LMTD = 10.0
        for _ in range(60):
            dT1 = T1 - T_H2O_out_K
            dT2 = T_CO2_out_K - T2
            if dT1 <= 0 or dT2 <= 0:
                break
            if abs(dT1 - dT2) < 1e-6:
                LMTD = (dT1 + dT2) / 2
            else:
                LMTD = (dT1 - dT2) / math.log(dT1 / dT2)
            T_CO2_avg = (T1 + T_CO2_out_K) / 2
            T_H2O_avg = (T2 + T_H2O_out_K) / 2
            q_guess = -abs(Q_total / As_CO2_total) if Q_total != 0 else -3000.0
            DH = DHC(T_CO2_avg, T_CO2_avg - 2, P1, d_CO2, G1, q_guess, F1, 0, ep_val, fg)
            q_h2o_guess = abs(Q_total / As_H2O_total) if Q_total != 0 else 3000.0
            Gne = Gniel(T_H2O_avg, T_H2O_avg + 2, P2, a, G2, da, L, q_h2o_guess, F2)
            htc_co2 = max(DH.htc[0], 1.0)
            htc_h2o = max(Gne.htc[0], 1.0)
            R1 = 1 / (htc_co2 * As_CO2_total)
            R3 = 1 / (htc_h2o * As_H2O_total)
            U = 1 / As_CO2_total * (1 / (R1 + R2_total + R3))
            Q_total = U * As_CO2_total * LMTD
            h_CO2_in = CP.PropsSI('H', 'T', T1, 'P', P1, 'REFPROP::CO2')
            h_CO2_out = h_CO2_in - Q_total / mdot_1
            h_H2O_in = CP.PropsSI('H', 'T', T2, 'P', P2, 'REFPROP::H2O')
            h_H2O_out = h_H2O_in + Q_total / mdot_2
            try:
                T_CO2_out_K = CP.PropsSI('T', 'H', h_CO2_out, 'P', P1, 'REFPROP::CO2')
                T_H2O_out_K = CP.PropsSI('T', 'H', h_H2O_out, 'P', P2, 'REFPROP::H2O')
            except Exception:
                break
        q_CO2_whole = Q_total / As_CO2_total
        q_H2O_whole = -Q_total / As_H2O_total
        n_wh = 2
        x_wh = np.array([0.0, L])
        plot_whole = pd.DataFrame({
            'pos': x_wh,
            'roughness': [ep_val] * n_wh,
            'CO2b': [T1, T_CO2_out_K],
            'CO2w': [T1 - 1, T_CO2_out_K - 1],
            'H2Ob': [T_H2O_out_K, T2],
            'H2Ow': [T_H2O_out_K + 1, T2 + 1],
            'Q': [Q_total / 2, Q_total / 2],
            'q_CO2': [q_CO2_whole, q_CO2_whole],
            'q_H2O': [q_H2O_whole, q_H2O_whole],
            'lmtd': [LMTD, LMTD],
            'Nu_CO2': [DH.Nu[0], DH.Nu[0]],
            'htc_CO2': [htc_co2, htc_co2],
            'Nu_H2O': [Gne.Nu[0], Gne.Nu[0]],
            'htc_H2O': [htc_h2o, htc_h2o],
            'CO2_mCp': [mdot_1 * CP.PropsSI('C', 'T', T1, 'P', P1, 'REFPROP::CO2'), mdot_1 * CP.PropsSI('C', 'T', T_CO2_out_K, 'P', P1, 'REFPROP::CO2')],
            'H2O_mCp': [mdot_2 * CP.PropsSI('C', 'T', T_H2O_out_K, 'P', P2, 'REFPROP::H2O'), mdot_2 * CP.PropsSI('C', 'T', T2, 'P', P2, 'REFPROP::H2O')],
            'r_CO2': [R1, R1],
            'r_cond': [R2_total, R2_total],
            'r_H2O': [R3, R3],
            'U': [U, U],
            'CO2_visc': [0, 0], 'CO2_acc': [0, 0], 'CO2_tot': [0, 0], 'Ve_CO2': [G1 / CP.PropsSI('D', 'T', T1, 'P', P1, 'REFPROP::CO2'), 0],
            'H2O_visc': [0, 0], 'H2O_acc': [0, 0], 'H2O_tot': [0, 0], 'Ve_H2O': [G2 / CP.PropsSI('D', 'T', T2, 'P', P2, 'REFPROP::H2O'), 0],
            'CO2_dp_cumulative_psi': [0, 0],
        })
        plot_whole['C'] = plot_whole[['CO2_mCp', 'H2O_mCp']].min(axis=1) / plot_whole[['CO2_mCp', 'H2O_mCp']].max(axis=1)
        plot_whole['NTU'] = (plot_whole['U'] * As_CO2_total) / plot_whole[['CO2_mCp', 'H2O_mCp']].min(axis=1)
        plot_whole['effectiveness'] = 50.0
        plot_whole['overall effectiveness'] = 50.0
        print(f"Whole section: Q_total = {Q_total:.1f} W, q_CO2 = {q_CO2_whole:.0f} W/m², q_H2O = {abs(q_H2O_whole):.0f} W/m²")
        print(f"  T_CO2 out = {T_CO2_out_K - 273.15:.2f} °C, T_H2O out = {T_H2O_out_K - 273.15:.2f} °C")
        return plot_whole

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

    #% roughness profile (Constant, Step, Linear, Sine)
    rough = np.zeros(n)
    if pro_key == 'Constant':
        rough = ep_max_input * np.ones(n, dtype=np.float64)
    elif pro_key == 'Step':
        num = np.linspace(0, ep_max_input, int(round(p)))
        step = np.linspace(0, n, int(round(p)) + 1, dtype=int)
        for i in range(int(round(p))):
            rough[step[i]:step[i + 1]] = num[i]
    elif pro_key == 'Linear':
        rough = ep_max_input * (x) / L
    elif pro_key == 'Sine':
        rough = ep_max_input * np.sin(p * np.pi * x / L)
    else:
        raise ValueError(
            f"Unknown roughness profile '{pro}'. "
            "Use Constant, Step, Linear, or Sine."
        )

    # ------------------------------------------------------------------
    # Heat transfer storage
    #   Segment arrays (length n, entry n-1 unused): Q, Q_leak, lmtd, dT
    #   Node arrays: Nu, htc, R1, R3, U, f, Tf
    #   Sign: Q < 0 means heat leaves CO2 (CO2 is cooled).
    # ------------------------------------------------------------------
    Q = np.zeros(n);      Q[:nseg] = q_CO2_i * As_CO2      # initial guess
    Q_leak = np.zeros(n)
    lmtd = np.zeros(n);   dT = np.zeros(n)
    Nu_CO2 = np.zeros(n); htc_CO2 = np.zeros(n)
    Nu_H2O = np.zeros(n); htc_H2O = np.zeros(n)
    R1 = np.zeros(n);     R3 = np.zeros(n);  U = np.zeros(n)
    f_CO2 = np.zeros(n);  f_H2O = np.zeros(n)
    Tf_CO2 = np.zeros(n); Tf_H2O = np.zeros(n)

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
                return hw - h_H2O_in                  # trial guess ran out of property range
            if not (T_lo < Tc1 < T_hi and T_lo < Tw1 < T_hi):
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
        return hw - h_H2O_in

    #% Outer Picard loop: correlations <-> counterflow thermal solve
    converged = False
    dT_change = np.nan
    dQ_rel = np.nan
    U_prev = None
    Tw_out_prev = None
    for j in range(max_iter):
        Tc_prev = CO2['Tb'].copy()
        Th2o_prev = H2O['Tb'].copy()
        Qtot_prev = Q[:nseg].sum()

        # node heat rate per segment-length, from current segment Q
        Qn = node_avg(Q[:nseg])

        # ---- 1. Correlations at every node (CO2 and water at the SAME position k) ----
        for k in range(n):
            ## CO2 side: Dang & Hihara
            DH = DHC(CO2['Tb'][k], CO2['Tw'][k], CO2['P'][k], d_CO2, G1,
                     Qn[k] / As_CO2, F1, k, float(rough[k]), fg)
            Nu_CO2[k], htc_CO2[k] = DH.Nu[0], DH.htc[0]
            f_CO2[k], Tf_CO2[k] = DH.f[0], DH.Tf[0]
            Tw_inner = DH.Tw[0]
            if k != 0:   # inlet node keeps wall = bulk, as in the original version
                ASCO2.update(CP.PT_INPUTS, CO2['P'][k], Tw_inner)
                store(CO2, k, ASCO2, 'w')

            ## Water side: Gnielinski (annulus); heat flux into water is -Q/As_H2O
            # Pass the FULL exchanger length: Gnielinski's (1 + (dh/L)^(2/3)) entrance factor is a
            # length-averaged correction from the flow inlet. Passing dx (old code) inflated the
            # water htc ~4x at dx = 10 mm and made it grid-dependent.
            Gne = Gniel(H2O['Tb'][k], H2O['Tw'][k], H2O['P'][k], a, G2, da, L,
                        -Qn[k] / As_H2O, F2)
            Nu_H2O[k], htc_H2O[k] = Gne.Nu[0], Gne.htc[0]
            f_H2O[k], Tf_H2O[k] = Gne.f[0], Gne.Tf[0]
            if k != n - 1:   # water inlet node keeps wall = bulk, as in the original version
                # Outer wall = inner wall minus conduction drop across the tube wall.
                # Heat leaving CO2 through the wall = -Qn  ->  Tw_outer = Tw_inner + Qn*R2
                Tw_outer = Tw_inner + Qn[k] * R2
                Tw_clamp = max(274.0, min(450.0, Tw_outer))   # keep water liquid in REFPROP
                ASH2O.update(CP.PT_INPUTS, H2O['P'][k], Tw_clamp)
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
            march(ASH2O.hmass(), True)       # final march stores the solution profiles

        # ---- 3. Convergence of the Picard loop ----
        dT_change = max(np.max(np.abs(CO2['Tb'] - Tc_prev)), np.max(np.abs(H2O['Tb'] - Th2o_prev)))
        Qtot = Q[:nseg].sum()
        dQ_rel = abs(Qtot - Qtot_prev) / max(abs(Qtot), 1e-12)
        if j >= 1 and dT_change < tol_T and dQ_rel < tol_Q:
            converged = True
            break

    n_sweeps = j + 1
    if verbose:
        status = "converged" if converged else "NOT converged"
        print(f"Counterflow solver: {n_sweeps} outer iterations ({status}; max dT change {dT_change:.2e} K, "
              f"rel dQ {dQ_rel:.2e})")
    if not converged:
        print("WARNING: solver hit max_iter without converging. Try lowering relax (e.g. 0.4) "
              "or raising max_iter.")

    # ------------------------------------------------------------------
    # Pressure drop and pumping power (on the converged state)
    # CO2 moves k-1 -> k ; water moves k+1 -> k
    # ------------------------------------------------------------------
    u_CO2 = G1 / CO2['rhob']
    u_H2O = G2 / H2O['rhob']
    CO2_Visc = np.zeros(n); CO2_acc = np.zeros(n)
    H2O_Visc = np.zeros(n); H2O_acc = np.zeros(n)
    for k in range(n):
        if k >= 1:
            rho_avg = 0.5 * (CO2['rhob'][k - 1] + CO2['rhob'][k])     # bulk density over segment
            CO2_Visc[k] = f_CO2[k] * dx / d_CO2 * G1**2 / (2 * rho_avg)
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

    q_CO2_seg = Q / As_CO2
    q_H2O_seg = -Q / As_H2O

    plot = pd.DataFrame({
        'pos'       : x,
        'roughness' : np.ravel(rough),
        'Re_CO2'    : Re_CO2b,
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

    #%% return results
    return plot