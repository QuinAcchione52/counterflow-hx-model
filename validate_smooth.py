"""Compare smooth-tube HTC with e/D=0; no HTC fitting or hydraulic correction."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
from CounterHX_full import match_bulk_case


def validate_one(values,segments):
    row=pd.Series(values)
    summary={'Case':row.Case,'Source test label':row.source_test_label,'Test group':row.test_group,
             'Spencer target HTC (W/m2 K)':row.target_htc_W_m2K,
             'Spencer mean CO2 temperature (C)':row.target_mean_CO2_C,
             'Reported HTC uncertainty (%)':row.experimental_htc_uncertainty_percent,
             'Nominal heat flux group (W/m2)':row.nominal_heat_flux_W_m2,
             'Assumed water inlet (C)':row['T_H2O (C)'],'Assumed water flow (kg/s)':row['m_H2O (kg/s)'],
             'CO2 mass flow (kg/s)':row['m_CO2 (kg/s)'],'Inlet pressure (MPa)':row['P (Mpa)'],
             'Fixed e/D':0.,'Segments':segments}
    try:
        p=match_bulk_case(row,0.,float(row.target_mean_CO2_C),segments)
        s=p.attrs['summary'];h=float(s['htc_spencer_W_m2K']);error=100*(h/row.target_htc_W_m2K-1)
        summary.update({'Numerical HTC (W/m2 K)':h,'Signed HTC error (%)':error,
                        'Absolute HTC error (%)':abs(error),'HTC within 5 percent':abs(error)<=5,
                        'Within reported HTC uncertainty':abs(error)<=row.experimental_htc_uncertainty_percent,
                        'Numerical mean CO2 temperature (C)':s['CO2_mean_endpoint_C'],
                        'Inferred CO2 inlet (C)':s['inferred_CO2_inlet_C'],
                        'Numerical CO2 outlet (C)':s['CO2_out_C'],
                        'Numerical heat flux (W/m2)':s['heat_flux_CO2_W_m2'],
                        'Difference from nominal flux (%)':100*(s['heat_flux_CO2_W_m2']/row.nominal_heat_flux_W_m2-1),
                        'Numerical pressure drop (psi)':s['CO2_pressure_drop_psi'],
                        'CO2 energy residual (W)':s['CO2_energy_residual_W'],
                        'Water energy residual (W)':s['water_energy_residual_W'],
                        'Property backend':s['property_backend'],'Status':'completed'})
        return summary,p
    except Exception as e:
        summary.update(Status='solver_error',Error=str(e))
        return summary,None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',default='Smooth_cases.csv');ap.add_argument('--output',default='Smooth-comparison.xlsx')
    ap.add_argument('--segments',type=int,default=100);ap.add_argument('--workers',type=int,default=1)
    args=ap.parse_args();data=pd.read_csv(args.input);records=[];profiles={};output=Path(args.output)
    if args.workers<1:ap.error('workers must be positive')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(validate_one,row.to_dict(),args.segments) for _,row in data.iterrows()]
        for future in as_completed(futures):
            record,profile=future.result();records.append(record)
            if profile is not None:profiles[record['Case']]=profile
            pd.DataFrame(records).to_csv(output.with_suffix('.partial.csv'),index=False)
            print(f"{len(records)}/{len(data)} {record['Case']}: {record['Status']}; "
                  f"HTC error={record.get('Signed HTC error (%)','unavailable')}",flush=True)
    order={case:i for i,case in enumerate(data.Case)};records.sort(key=lambda r:order[r['Case']]);table=pd.DataFrame(records)
    table.to_csv(output.with_suffix('.csv'),index=False,encoding='utf-8-sig')
    with pd.ExcelWriter(output,engine='xlsxwriter') as writer:
        table.to_excel(writer,sheet_name='Comparison',index=False);data.to_excel(writer,sheet_name='Inputs_assumptions',index=False)
        sheet=writer.sheets['Comparison'];sheet.freeze_panes(1,1);sheet.set_column(0,len(table.columns)-1,24)
        sheet.add_table(0,0,len(table),len(table.columns)-1,{'columns':[{'header':c} for c in table.columns]})
        for case in data.Case:
            if case in profiles:profiles[case].to_excel(writer,sheet_name=str(case)[:31],index=False)
        notes=pd.DataFrame({'Item':['Purpose','HTC','Boundary conditions','Limitations','Friction'],
                            'Description':['38 smooth runs from thesis Table 12; e/D=0; no HTC fitting',
                            'Experimental-definition endpoint bulk-to-inner-wall LMTD HTC',
                            'Published mean bulk held fixed by inferring CO2 inlet; water inlet20.7 C and flow0.087 kg/s assumed',
                            'Conditional base-model comparison, not validation with independently measured inlet conditions; nominal flux is not actual per-case measured flux',
                            'Same Colebrook zero-roughness extension and fully developed laminar water baseline as dimpled workflow']})
        notes.to_excel(writer,sheet_name='Read_me',index=False)
    ok=table[table.Status=='completed'];print(f"Completed {len(ok)}/{len(data)}; within 5%: {int(ok['HTC within 5 percent'].sum())}; "
          f"mean absolute HTC error: {ok['Absolute HTC error (%)'].mean():.2f}%",flush=True)


if __name__=='__main__':main()
