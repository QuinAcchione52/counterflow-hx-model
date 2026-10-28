"""One-axis red/blue comparison; unavailable fits marked in a separate status row."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/hx-redblue-mpl')
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROOT/'results/dimpled_trends/All_40_old_vs_new_eD.csv')
x=np.arange(len(d));matched=d.new_fitted_eD.notna()
fig,ax=plt.subplots(figsize=(15,6.7))
ax.scatter(x,d.old_fitted_eD,color='red',s=42,label='Original DH-based fitted e/D (40)',zorder=3)
ax.scatter(x[matched],d.loc[matched,'new_fitted_eD'],color='blue',s=42,label='Revised + smooth-trained multiplier: fitted e/D (24)',zorder=4)
ax.axhspan(-.042,-.01,color='#f0f0f0',zorder=0)
ax.scatter(x[~matched],np.full((~matched).sum(),-.026),marker='x',color='blue',s=45,label='Revised unfitted cases (16): status symbols, NOT e/D values',zorder=4)
ax.axhline(.05,color='gray',ls='--',lw=1,label='Tested revised e/D cap = 0.05')
ax.set_xticks(x,d.Case,rotation=90,fontsize=8)
ax.set_yticks(np.arange(0,.71,.1));ax.set_ylim(-.045,.70)
ax.set(xlabel='Dimpled case',ylabel='Fitted effective e/D',title='All 40 dimpled cases: original DH-based fits vs revised fits')
ax.grid(axis='y',alpha=.2);ax.legend(loc='upper right',fontsize=9)
fig.text(.075,.02,'Fitted e/D values are model-derived, not directly reported by Spencer. Blue × marks indicate no revised fit within 0–0.05.',fontsize=10)
fig.tight_layout(rect=(0,.045,1,1))
for ext in ['png','pdf']:fig.savefig(ROOT/f'results/dimpled_trends/All_40_eD_Red_vs_Blue.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
