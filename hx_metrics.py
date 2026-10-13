"""Experimental-comparison metrics (thesis equations 16, 18, 19)."""
import math
import numpy as np


def logarithmic_mean(a, b):
    if not (a > 0 and b > 0):
        raise ValueError('Bulk-to-wall endpoint differences must both be positive')
    if math.isclose(a, b, rel_tol=1e-8):
        return (a+b)/2
    return (a-b)/math.log(a/b)


def summarize(profile, inner_diameter_mm, wall_thickness_mm):
    length = float(profile.pos.iloc[-1]-profile.pos.iloc[0])
    area = math.pi*inner_diameter_mm*1e-3*length
    duty = -float(profile.Q.sum())
    delta = profile.CO2b-profile.CO2w
    lmtd = logarithmic_mean(float(delta.iloc[0]), float(delta.iloc[-1]))
    return {
        'htc_spencer_W_m2K': duty/area/lmtd,
        'htc_local_length_mean_W_m2K': float(np.trapezoid(profile.htc_CO2, profile.pos)/length),
        'duty_CO2_W': duty,
        'heat_flux_CO2_W_m2': duty/area,
        'bulk_wall_lmtd_K': lmtd,
        'CO2_out_C': float(profile.CO2b.iloc[-1]-273.15),
        'water_out_C': float(profile.H2Ob.iloc[0]-273.15),
        'CO2_mean_endpoint_C': float((profile.CO2b.iloc[0]+profile.CO2b.iloc[-1])/2-273.15),
        'CO2_pressure_drop_psi': float(profile.CO2_dp_cumulative_psi.iloc[-1]),
        'water_Re_inlet': float(profile.Re_H2O.iloc[-1]),
        'inner_diameter_mm': float(inner_diameter_mm),
        'outer_diameter_mm': float(inner_diameter_mm+2*wall_thickness_mm),
    }
