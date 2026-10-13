# -*- coding: utf-8 -*-
"""
Compiler and main code for Counter flow Heat exchanger prediction ensure that input parameter CSV file is correct with
conditions that are to be ran. Outputs are generated and stored in an Output.CSV file with each case having a 
different sheet within the Output file.

Inputs include: 
    Mass flow rate (kg/s): m_CO2, m_H2O
    Inlet Temperature (Ti): for CO2 and H2O
    Operating Pressure (P): for CO2 and H2O
    Heat Exchanger geometry: 
            Diameter (mm)
            wall thickness (mm)
            thermal conductivity (W/m-K)
            cooling channel thickness (mm)
            length of Heat Exchanger (mm)
            
Outputs include: (at each segment along the pipe)
    axial position along pipe
    Bulk and Wall temperatures (K): for CO2 and H2O
    Total Heat Transfered (W)
    Heat FLux (W/m^2): for CO2 and H2O
    Temperature difference between CO2 adn H2O (lmtd)
    Overall heat transfer coefficeint U
    Thermal resistances (K/W): R sCO2, R Wall, R H2O
    Heat compacity (W/K): for CO2 and H2O
    Heat Transfer Coeffcient (W/m^2-K): for CO2 and H2O
    Fluid velocity (m/s): for CO2 and H2O
    Pressure drop (Pa):
        due to Viscous effects
        due to acceleration effects
"""
import pandas as pd
import time
from Counterflow_HX_4 import CO2_counter

#df = pd.read_csv('Parameters.csv')
#reading the input csv file and creating a dataframe with the operating conditions desired for each case
name= 'Test'

# Constant=False: compute Q from LMTD (required to match experiment). Use True only to impose fixed q.
Constant = False
# Whole section: one LMTD and one Q for full length (no segments). Set True to try lumped model.
whole_section = False
# Initial guess for q into CO2 (W/m²). Only used when Constant=False; should be near expected flux (e.g. ~2800 for Test 7a).
q_CO2_i = -3000.0

try:
    df = pd.read_csv(f"{name}.csv", encoding="utf-8")
except UnicodeDecodeError:
    df = pd.read_csv(f"{name}.csv", encoding="cp1252")
df.columns = [c.replace("\xa0", " ").strip() for c in df.columns]
df.rename(columns={'P (Mpa)': 'P_CO2', 'P (atm)': 'P_H2O', 'T_CO2 (C)': 'T_CO2', 'T_H2O (C)': 'T_H2O',
                      'D (mm)': 'D', 'wall thickness (mm)': 'wt','water gap (mm)': 'wg', 'L (mm)': 'L',
                      'ep_max (e/Dh)': 'ep', 'Profile (Constant,Step,Sine,Linear)': 'pro','Steps/Harmonic (n)': 'p',
                      'm_CO2 (kg/s)': 'm_CO2', 'm_H2O (kg/s)': 'm_H2O'}, inplace=True)
# Accept m_CO2/m_H2O if already named that (no units in header)
if 'm_CO2' not in df.columns:
    for c in df.columns:
        if 'm_CO2' in c or (c.strip().lower().startswith('m_') and 'co2' in c.lower()):
            df.rename(columns={c: 'm_CO2'}, inplace=True)
            break
if 'm_H2O' not in df.columns:
    for c in df.columns:
        if 'm_H2O' in c or (c.strip().lower().startswith('m_') and 'h2o' in c.lower()):
            df.rename(columns={c: 'm_H2O'}, inplace=True)
            break
if 'm_CO2' not in df.columns or 'm_H2O' not in df.columns:
    raise ValueError("CSV must have mass flow rate columns 'm_CO2' and 'm_H2O' (kg/s). Current columns: " + str(list(df.columns)))
# Optional heat leakage: linear model Q_leak = slope*(T_amb - T_water) + intercept (from your fit),
# or U_amb model: CSV columns T_amb (C), U_amb (W/m2.K), d_outer_mm. Prefer linear if Q_leak_slope present.
fullrun = []
cases=len(df)
Ts = time.time()


with pd.ExcelWriter(f"{name}-output.xlsx",engine="xlsxwriter") as writer:
    for i in range(cases):
        t0 = time.time()
        profile = df.pro[i]
        #Running Counter flow Heat Exchanger code for each case in the provided 
        #data frame
        print(
            f"\nCase {i + 1}/{cases}: "
            f"Profile = {df.pro[i]}, "
            f"e/D = {df.ep[i]:.6f}",
            flush=True
        )
        T_amb = df['T_amb (C)'].iloc[i] if 'T_amb (C)' in df.columns else 25.0
        U_amb = df['U_amb'].iloc[i] if 'U_amb' in df.columns else 0.0
        d_outer = df['d_outer_mm'].iloc[i] if 'd_outer_mm' in df.columns else 0.0
        Q_leak_slope = df['Q_leak_slope'].iloc[i] if 'Q_leak_slope' in df.columns and pd.notna(df['Q_leak_slope'].iloc[i]) else None
        Q_leak_intercept = float(df['Q_leak_intercept'].iloc[i]) if 'Q_leak_intercept' in df.columns and pd.notna(df['Q_leak_intercept'].iloc[i]) else 0.0
        results = CO2_counter(df['m_CO2'].iloc[i],df.T_CO2[i],df.P_CO2[i],df['m_H2O'].iloc[i],df.T_H2O[i],df.P_H2O[i],
                          df.D[i],df.wt[i],df.wg[i],df.L[i],df.Kp[i],df.ep[i],str(df.pro[i]),df.p[i],
                          Constant,q_CO2_i,whole_section,T_ambient=T_amb,U_amb=U_amb,d_outer_mm=d_outer,
                          Q_leak_slope=Q_leak_slope,Q_leak_intercept=Q_leak_intercept)
        #save results to post process immediately
        fullrun.append(results)
         
        #writing results to Excel file to be stored and saved for late
        sheet = df.Case[i]
        results.to_excel(writer, sheet_name=sheet, index=False)
        #timing each iteration
        t1 = time.time()
        elapsed = t1-t0
        print(f"Case {i+1} done in {elapsed:.1f} s", flush=True)

Tf = time.time()
seconds = Tf-Ts
minutes= seconds/60
print(f"\nAll done: {cases} case(s) in {seconds:.1f} s ({minutes:.2f} min)", flush=True)
print(f"Output written to {name}-output.xlsx", flush=True)
    
