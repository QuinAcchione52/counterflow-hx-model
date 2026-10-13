"""Dimensionless local e/D profiles defined in physical axial coordinates."""
import numpy as np


def build_roughness(x, ep, profile='Constant', p=1):
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.all(np.isfinite(x)) or np.any(np.diff(x)<=0):
        raise ValueError('x must contain at least two finite, increasing positions')
    z=(x-x[0])/(x[-1]-x[0])
    values=np.asarray(ep,dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values<0) or np.any(values>1):
        raise ValueError('Every local e/D must be finite and between 0 and 1')
    if values.ndim:
        if values.ndim!=1 or len(values)!=len(x):
            raise ValueError('A custom e/D array must have one value per axial node')
        return values.copy()
    maximum=float(values)
    name=str(profile).strip().lower()
    if name=='constant':return np.full(len(x),maximum)
    if name=='linear':return maximum*z
    if name=='step':
        if not np.isfinite(p) or p<1 or int(p)!=p:
            raise ValueError('Step count must be a positive integer')
        count=int(p)
        levels=np.array([maximum]) if count==1 else np.linspace(0,maximum,count)
        index=np.minimum(np.floor(z*count+1e-12).astype(int),count-1)
        return levels[index]
    if name=='sine':
        if not np.isfinite(p) or p<=0:raise ValueError('Sine harmonic must be positive and finite')
        return maximum*np.abs(np.sin(float(p)*np.pi*z))
    raise ValueError('Profile must be Constant, Step, Linear, Sine, or a custom node array')
