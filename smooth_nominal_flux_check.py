"""CO2-only smooth diagnostic at nominal flux; water-side assumptions bypassed."""
import os,math
import argparse
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP
from scipy.optimize import brentq
from Correlations import DHC
from hx_metrics import logarithmic_mean

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--input',default='Smooth_cases.csv')
 ap.add_argument('--output',default='Smooth-CO2-only-nominal-flux.csv')
 ap.add_argument('--segments',type=int,default=100)
 args=ap.parse_args()
 if args.segments<1:ap.error('segments must be positive')
 records=[]
 for _,row in pd.read_csv(args.input).iterrows():
  result={'Case':row.Case,'Spencer HTC (W/m2 K)':row.target_htc_W_m2K}
  try:
   state=CP.AbstractState(os.environ.get('HX_PROPERTY_BACKEND','HEOS'),'CO2')
   P=row['P (Mpa)']*1e6;D=row['D (mm)']*1e-3;L=row['L (mm)']*1e-3;m=row['m_CO2 (kg/s)'];Tmean=row.target_mean_CO2_C+273.15
   q=row.nominal_heat_flux_W_m2;Q=q*math.pi*D*L
   def enthalpy(t):state.update(CP.PT_INPUTS,P,float(t));return state.hmass()
   def balance(t):return enthalpy(t)-enthalpy(2*Tmean-t)-Q/m
   inlet=brentq(balance,Tmean,Tmean+20,xtol=1e-9)
   x=np.linspace(0,L,args.segments+1);h=enthalpy(inlet)-Q/m*x/L;bulk=[];wall=[];htc=[]
   G=m/(math.pi*D**2/4)
   for k,energy in enumerate(h):
    state.update(CP.HmassP_INPUTS,float(energy),P);t=state.T()
    s=DHC(t,t-2,P,D,G,-q,'CO2',k,0.,.02)
    bulk.append(t);wall.append(s.Tw.iloc[0]);htc.append(s.htc.iloc[0])
   lm=logarithmic_mean(bulk[0]-wall[0],bulk[-1]-wall[-1]);prediction=q/lm
   result.update({'CO2-only HTC at nominal flux (W/m2 K)':prediction,
    'CO2-only HTC error (%)':100*(prediction/row.target_htc_W_m2K-1),
    'Nominal-flux-implied inlet (C)':inlet-273.15,'Nominal-flux-implied outlet (C)':bulk[-1]-273.15,
    'Nominal heat flux (W/m2)':q,'Status':'completed'})
  except Exception as e:result.update(Status='error',Error=str(e))
  records.append(result)
 pd.DataFrame(records).to_csv(args.output,index=False)
 print(pd.DataFrame(records)[['Case','CO2-only HTC error (%)']].to_string(index=False))

if __name__=='__main__':main()
