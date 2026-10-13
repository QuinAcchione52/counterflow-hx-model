"""Separate pressure-drop calibration on an unchanged thermal solution.

The hydraulic multiplier scales viscous losses only; acceleration retains its
physical sign. It is empirical, and does not change the HTC roughness or inlet.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
import os
from scipy.optimize import brentq
from CounterHX_full import run_case, match_bulk_case

PSI_TO_PA = 6894.757


def fit_hydraulic_multiplier(profile, target_psi):
    if not np.isfinite(target_psi) or target_psi <= 0:
        raise ValueError('Target pressure drop must be positive and finite (psi)')
    viscous = float(profile.CO2_visc.sum())
    acceleration = float(profile.CO2_acc.sum())
    if viscous <= 0:
        raise ValueError('Positive baseline viscous loss required')
    multiplier = (target_psi * PSI_TO_PA - acceleration) / viscous
    if multiplier <= 0:
        raise ValueError('Target would require negative viscous loss')
    fitted = profile.copy(deep=True)
    fitted.attrs = {'summary': dict(profile.attrs['summary'])}
    fitted['CO2_visc'] *= multiplier
    fitted['CO2_tot'] = fitted.CO2_visc + fitted.CO2_acc
    fitted['CO2_dp_cumulative_psi'] = fitted.CO2_tot.cumsum() / PSI_TO_PA
    # Inlet dynamic pressure from G = rho*u and the known hydraulic diameter.
    summary = fitted.attrs['summary']
    summary.update(
        pressure_fit_is_empirical=True,
        hydraulic_viscous_multiplier=multiplier,
        baseline_viscous_drop_psi=viscous/PSI_TO_PA,
        acceleration_drop_psi=acceleration/PSI_TO_PA,
        baseline_total_drop_psi=(viscous+acceleration)/PSI_TO_PA,
        target_pressure_drop_psi=target_psi,
        CO2_pressure_drop_psi=float(fitted.CO2_dp_cumulative_psi.iloc[-1]),
        pressure_fit_relative_error=float(fitted.CO2_dp_cumulative_psi.iloc[-1]/target_psi-1),
    )
    return multiplier, fitted



def fit_pressure_only_roughness(profile, row, target_psi):
    """Invert Colebrook on a frozen thermal solution; diagnostic, not geometry.

    This independent hydraulic e/D is allowed up to 3 only to reveal how
    unphysical the required equivalent roughness is. Never feed it to DHC.
    """
    if not np.isfinite(target_psi) or target_psi <= 0:
        raise ValueError('Target pressure drop must be positive and finite (psi)')
    diameter = float(row['D (mm)']) * 1e-3
    flux = float(row['m_CO2 (kg/s)']) / (np.pi*diameter**2/4)
    density = flux / profile.Ve_CO2.to_numpy()
    film_temperature = (profile.CO2b.to_numpy()+profile.CO2w.to_numpy())/2
    state = CP.AbstractState(os.environ.get('HX_PROPERTY_BACKEND','HEOS'),'CO2')
    viscosity = []
    for temperature in film_temperature:
        state.update(CP.PT_INPUTS, float(row['P (Mpa)'])*1e6, float(temperature))
        viscosity.append(state.viscosity())
    reynolds = flux*diameter/np.asarray(viscosity)
    dx = np.diff(profile.pos.to_numpy())
    rho_avg = (density[:-1]+density[1:])/2
    acceleration = float(profile.CO2_acc.sum())
    def viscous(ep):
        f = np.full(len(reynolds), .02)
        for _ in range(200):
            updated = (-2*np.log10(ep/3.7+2.51/(reynolds*np.sqrt(f))))**-2
            if np.max(np.abs(updated-f)) < 1e-12:
                f=updated; break
            f=updated
        else:
            raise RuntimeError('Hydraulic Colebrook iteration did not converge')
        return np.r_[0., f[1:]*dx/diameter*flux**2/(2*rho_avg)]
    def residual(ep):
        return (float(viscous(ep).sum())+acceleration)/PSI_TO_PA-target_psi
    low, high=residual(0.), residual(3.)
    if low*high>0:
        raise ValueError('Pressure target not bracketed even by diagnostic hydraulic e/D=0..3')
    ep=brentq(residual,0.,3.,xtol=1e-10)
    return ep, (float(viscous(ep).sum())+acceleration)/PSI_TO_PA

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='Test.csv')
    parser.add_argument('--reference', default='experimental_targets.csv')
    parser.add_argument('--output', default='Pressure-fit.xlsx')
    parser.add_argument('--roughness', type=float, help='Thermal effective e/D; held fixed during pressure fit')
    parser.add_argument('--segments', type=int, default=100)
    parser.add_argument('--target-pressure-psi', type=float, help='Override target for a single case')
    args = parser.parse_args()
    cases = pd.read_csv(args.input)
    references = pd.read_csv(args.reference).set_index('Case')
    if not references.index.is_unique:
        parser.error('Reference Case names must be unique')
    if args.target_pressure_psi is not None and len(cases) != 1:
        parser.error('--target-pressure-psi requires a single-case input')
    summaries, profiles = [], []
    for _, row in cases.iterrows():
        reference = references.loc[row.Case]
        ep = args.roughness if args.roughness is not None else row['ep_max (e/Dh)']
        base = match_bulk_case(row, ep, reference.mean_CO2_C, args.segments)
        target = args.target_pressure_psi if args.target_pressure_psi is not None else reference.pressure_drop_psi
        multiplier, fitted = fit_hydraulic_multiplier(base, target)
        hydraulic_ep, hydraulic_pressure = fit_pressure_only_roughness(base,row,target)
        summary = {'Case':row.Case,'thermal_effective_e_over_D':ep,**fitted.attrs['summary'],
                   'pressure_only_equivalent_e_over_D':hydraulic_ep,
                   'pressure_only_equivalent_roughness_mm':hydraulic_ep*row['D (mm)'],
                   'pressure_from_equivalent_roughness_psi':hydraulic_pressure}
        fitted.attrs['summary'] = summary
        profiles.append((row.Case,fitted)); summaries.append(summary)
        print(f'{row.Case}: hydraulic multiplier={multiplier:.6g}; pressure={summary["CO2_pressure_drop_psi"]:.6g} psi; '
              f'unchanged HTC={summary["htc_spencer_W_m2K"]:.2f} W/(m² K)',flush=True)
    with pd.ExcelWriter(args.output,engine='xlsxwriter') as writer:
        for case,profile in profiles:
            profile.to_excel(writer,sheet_name=str(case)[:31],index=False)
        pd.DataFrame(summaries).to_excel(writer,sheet_name='Summary',index=False)
    pd.DataFrame(summaries).to_csv(Path(args.output).with_suffix('.summary.csv'),index=False)


if __name__ == '__main__':
    main()
