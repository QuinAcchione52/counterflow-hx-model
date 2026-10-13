"""Run cases or fit effective e/D using Spencer's bulk-to-wall LMTD HTC."""
import argparse
from pathlib import Path
import pandas as pd
from scipy.optimize import brentq
from Counterflow_HX_4 import CO2_counter


def run_case(row, roughness=None, segments=100, verbose=False):
    def optional(name, default=None):
        value = row.get(name, default)
        return default if pd.isna(value) else value
    return CO2_counter(
        row['m_CO2 (kg/s)'], row['T_CO2 (C)'], row['P (Mpa)'],
        row['m_H2O (kg/s)'], row['T_H2O (C)'], row['P (atm)'],
        row['D (mm)'], row['wall thickness (mm)'], row['water gap (mm)'],
        row['L (mm)'], row['Kp'],
        row['ep_max (e/Dh)'] if roughness is None else roughness,
        row['Profile (Constant,Step,Sine,Linear)'], row['Steps/Harmonic (n)'],
        False, -3000., T_ambient=optional('T_amb (C)', 25.),
        U_amb=optional('U_amb', 0.), d_outer_mm=optional('d_outer_mm', 0.),
        Q_leak_slope=optional('Q_leak_slope'),
        Q_leak_intercept=optional('Q_leak_intercept', 0.),
        segments=segments, verbose=verbose, max_iter=200)



def match_bulk_case(row, roughness, target_bulk, segments=100):
    """Infer the unknown CO2 inlet for a measured endpoint-mean bulk temperature.

    Water inlet remains the supplied condition. This is an inverse model,
    not a recovered experimental measurement.
    """
    if not target_bulk > row['T_H2O (C)']:
        raise ValueError('Target bulk temperature must exceed water inlet temperature')
    cache = {}
    def residual(inlet):
        trial = row.copy()
        trial['T_CO2 (C)'] = inlet
        profile = run_case(trial, roughness, segments)
        cache[inlet] = profile
        return profile.attrs['summary']['CO2_mean_endpoint_C'] - target_bulk
    lower = float(target_bulk)
    upper = max(lower + 1., float(row['T_CO2 (C)']))
    low = residual(lower)
    high = residual(upper)
    while low*high > 0 and upper < lower+20:
        upper = min(lower+20, upper+2.)
        high = residual(upper)
    if low*high > 0:
        raise ValueError('Could not bracket an inlet for the target bulk temperature')
    inlet = brentq(residual, lower, upper, xtol=1e-5)
    if inlet not in cache:
        residual(inlet)
    profile = cache[inlet]
    profile.attrs['summary'].update(inferred_CO2_inlet_C=inlet,
                                    target_bulk_C=target_bulk,
                                    CO2_inlet_is_inferred=True)
    if abs(profile.attrs['summary']['CO2_mean_endpoint_C']-target_bulk)>1e-4:
        raise RuntimeError('Inferred inlet failed bulk-temperature acceptance tolerance')
    return profile

def fit_case(row, target, bounds=(0., .05), segments=100, target_bulk=None, allow_extrapolation=False):
    """Bracket first; refuse out-of-range targets. Never fit unconverged solutions."""
    if target <= 0 or not 0 <= bounds[0] < bounds[1] <= 1:
        raise ValueError('Positive HTC target and 0 <= lower < upper <= 1 required')
    if bounds[1] > .05 and not allow_extrapolation:
        raise ValueError('Bounds above e/D=0.05 require explicit allow_extrapolation; these are empirical, not physical fits')
    cache = {}
    def evaluate(ep):
        if ep not in cache:
            cache[ep] = (run_case(row, ep, segments) if target_bulk is None
                         else match_bulk_case(row, ep, target_bulk, segments))
        return cache[ep].attrs['summary']['htc_spencer_W_m2K']-target
    low, high = [evaluate(ep) for ep in bounds]
    if low*high > 0:
        raise ValueError(f'Target {target:g} not bracketed: HTC at e/D={bounds[0]:g} is '
                         f'{low+target:.1f}; at {bounds[1]:g} is {high+target:.1f} W/(m² K). '
                         'Check inputs and correlation applicability; do not force a match.')
    ep = brentq(evaluate, *bounds, xtol=1e-7)
    evaluate(ep)
    result = cache[ep]
    if abs(result.attrs['summary']['htc_spencer_W_m2K']-target)/target > .001:
        raise RuntimeError('Fitted HTC failed the 0.1% acceptance tolerance')
    return ep, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='Test.csv')
    parser.add_argument('--output', default='Test-output.xlsx')
    parser.add_argument('--segments', type=int, default=100)
    inlet_options = parser.add_mutually_exclusive_group()
    inlet_options.add_argument('--match-bulk', dest='match_bulk', action='store_true', default=None,
                               help='Infer CO2 inlet to match target_mean_CO2_C or reference mean_CO2_C')
    inlet_options.add_argument('--fixed-inlet', dest='match_bulk', action='store_false',
                               help='Use the CSV inlet as a fixed measured condition instead of inferring it')
    parser.add_argument('--bulk-temp', type=float, help='Single-case measured endpoint-mean CO2 temperature (C)')
    parser.add_argument('--roughness', type=float, help='Override CSV effective e/D')
    parser.add_argument('--reference', help='Reference CSV keyed by Case, e.g. experimental_targets.csv')
    parser.add_argument('--fit', action='store_true', help='Fit each target_htc_W_m2K in CSV')
    parser.add_argument('--target-htc', type=float, help='Target HTC for a single-case CSV')
    parser.add_argument('--allow-extrapolation', action='store_true', help='Allow empirical e/D fits beyond 0.05; not a physical roughness estimate')
    parser.add_argument('--ep-min', type=float, default=0.)
    parser.add_argument('--ep-max', type=float, default=.05)
    args = parser.parse_args()
    if args.fit and args.roughness is not None:
        parser.error('--fit and --roughness cannot be combined')
    df = pd.read_csv(args.input, encoding='utf-8-sig')
    df.columns = df.columns.str.replace('\xa0',' ').str.strip()
    if args.target_htc is not None and (not args.fit or len(df) != 1):
        parser.error('--target-htc requires --fit and a single-case CSV')
    if args.reference:
        reference = pd.read_csv(args.reference)
        if reference.Case.duplicated().any():
            parser.error('Reference Case names must be unique')
        # Keep experimental observations separate from solver inputs.
        for _, ref in reference.iterrows():
            for key in ['target_htc_W_m2K', 'mean_CO2_C', 'pressure_drop_psi']:
                if key in ref:
                    df.loc[df.Case == ref.Case, 'experimental_' + key] = ref[key]
        if 'target_htc_W_m2K' not in df:
            df['target_htc_W_m2K'] = df['experimental_target_htc_W_m2K']
    if args.bulk_temp is not None and (not args.match_bulk or len(df) != 1):
        parser.error('--bulk-temp requires --match-bulk and a single-case CSV')
    summaries=[]
    # Finish calculations before opening an output workbook, so a failed fit does not overwrite it.
    results=[]
    for _,row in df.iterrows():
        target=None
        target_bulk=None
        match_bulk = (args.match_bulk if args.match_bulk is not None
                      else pd.notna(row.get('target_mean_CO2_C', float('nan'))))
        if match_bulk:
            target_bulk = (args.bulk_temp if args.bulk_temp is not None
                           else row.get('target_mean_CO2_C', row.get('experimental_mean_CO2_C')))
            if target_bulk is None or pd.isna(target_bulk):
                parser.error('--match-bulk requires --bulk-temp or a matched reference mean_CO2_C')
        if args.fit:
            target=args.target_htc if args.target_htc is not None else row['target_htc_W_m2K']
            ep,result=fit_case(row,target,(args.ep_min,args.ep_max),args.segments,target_bulk,args.allow_extrapolation)
        else:
            ep=row['ep_max (e/Dh)'] if args.roughness is None else args.roughness
            result=(run_case(row,ep,args.segments,verbose=True) if target_bulk is None
                    else match_bulk_case(row,ep,target_bulk,args.segments))
        summary={'Case':row['Case'],'effective_e_over_D':ep,
                 'roughness_extrapolated':bool(ep > .05),**result.attrs['summary']}
        for key in ['experimental_target_htc_W_m2K', 'experimental_mean_CO2_C', 'experimental_pressure_drop_psi']:
            if key in row and pd.notna(row[key]):
                summary[key] = row[key]
        if target is not None:
            summary.update(target_htc_W_m2K=target,relative_htc_error=summary['htc_spencer_W_m2K']/target-1)
        summaries.append(summary); results.append((row['Case'],result))
        print(f"{row['Case']}: e/D={ep:.7g}, Spencer-definition HTC={summary['htc_spencer_W_m2K']:.1f}",flush=True)
    with pd.ExcelWriter(args.output,engine='xlsxwriter') as writer:
        for case,result in results:
            result.to_excel(writer,sheet_name=str(case)[:31],index=False)
        pd.DataFrame(summaries).to_excel(writer,sheet_name='Summary',index=False)
    pd.DataFrame(summaries).to_csv(Path(args.output).with_suffix('.summary.csv'),index=False)


if __name__=='__main__':
    main()
