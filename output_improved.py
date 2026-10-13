import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D




#importing Excel file
name = 'Test'
input_file = pd.read_csv(f"{name}.csv")
output_file = pd.read_excel(f"{name}-output.xlsx", sheet_name=None)

input_file.columns = [c.replace("\xa0", " ").strip() for c in input_file.columns]
# Filename prefix from first case (e.g. T_CO2 and P for figures)
row0 = input_file.iloc[0]
t_co2 = row0.get('T_CO2 (C)', row0.get('T_CO2', 0))
p_co2 = row0.get('P (Mpa)', row0.get('P_CO2', 0))
fname_prefix = f"{name}_T{t_co2}_P{p_co2}"

line_styles = ['-', ':', '-.', '--']
profile_names = ['Case 1', 'Case 2', 'Case 3', 'Case 4']
legend_handles = [
    Line2D([0], [0], color='black', linestyle=line_styles[i % len(line_styles)], linewidth=2, label=profile_names[i])
    for i in range(min(4, len(output_file)))
]
#%% Roughness Profile
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos']
    CO2 = df['roughness']

    ax.plot(x_data, CO2, color='black', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Roughness Value ($\epsilon/d$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_roughness.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%%Temperature profile plot
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos']
    CO2 = df['CO2b']
    H2O = df['H2Ob']

    ax.plot(x_data, CO2, color='green', linestyle=line_styles[idx % len(line_styles)])
    ax.plot(x_data, H2O, color='blue', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Temperature ($K$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_temperature.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Log Mean Temperature Difference (LMTD)
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][0:-1]
    CO2 = df['lmtd'][0:-1]

    ax.plot(x_data, CO2, color='r', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Log Mean Temperature Difference ($K$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_LMTD.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Heat Transfer Coefficient (HTC)
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][1:-1]
    CO2 = df['htc_CO2'][1:-1]
    H2O = df['htc_H2O'][1:-1]

    ax.plot(x_data, CO2, color='green', linestyle=line_styles[idx % len(line_styles)])
    ax.plot(x_data, H2O, color='blue', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Heat transfer coefficient (W/m$^2$ K)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_HTC.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Local pressure drop  (dP Total )
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][1:]
    CO2 = df['CO2_dp_cumulative_psi'][1:]
    H2O = None

    ax.plot(x_data, CO2, color='green', linestyle=line_styles[idx % len(line_styles)])
    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Cumulative CO2 pressure drop (psi)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_dP.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Segment Effectiveness
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][1:-1]
    eff = df['effectiveness'][1:-1]

    ax.plot(x_data, eff, color='red', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Segment Effectiveness ($\%$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_effectiveness.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

# %% overall Effectiveness
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos']
    eff = df['overall effectiveness']

    ax.plot(x_data, eff, color='red', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Overall Effectiveness ($\%$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_overall_effectiveness.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Heat Capacitance  (mdotCp)
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][1:-1]
    CO2 = df['CO2_mCp'][1:-1]
    H2O = df['H2O_mCp'][1:-1]

    ax.plot(x_data, CO2, color='green', linestyle=line_styles[idx % len(line_styles)])
    ax.plot(x_data, H2O, color='blue', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Heat Capacitance ($W/K$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_heat_capacitance.png", dpi=150, bbox_inches='tight')
plt.show(block=False)

#%% Thermal Resistance
fig, ax = plt.subplots()

for idx, df in enumerate(output_file.values()):
    x_data = df['pos'][1:-1]
    CO2 = df['r_CO2'][1:-1]
    H2O = df['r_H2O'][1:-1]
    conduction = df['r_cond'][1:-1]

    ax.plot(x_data, CO2, color='green', linestyle=line_styles[idx % len(line_styles)])
    ax.plot(x_data, H2O, color='blue', linestyle=line_styles[idx % len(line_styles)])
    ax.plot(x_data, conduction, color='#800000', linestyle=line_styles[idx % len(line_styles)])

    ax.set_xlabel('Axial position (m)',fontsize=12)
    ax.set_ylabel(r'Thermal Resistance ($K/W$)',fontsize=12)
    legend = ax.legend(handles=legend_handles, loc="best", fontsize=10, frameon=True)
    ax.add_artist(legend)
plt.savefig(f"{fname_prefix}_thermal_resistance.png", dpi=150, bbox_inches='tight')
plt.show()