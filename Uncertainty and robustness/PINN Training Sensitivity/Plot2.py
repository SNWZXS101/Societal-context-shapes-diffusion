
import os
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# --- Configuration ---
RESULTS_DIR   = 'results'
TS_FILE       = 'time_series.xlsx'  # must contain Day, Se, S, I, R, f
RAW_FILE      = 'daily_cases.xlsx'
FIG_FILE      = 'Figs_single_instance.png'
SMOOTH_WINDOW = 5
M_SAMPLES     = 200   # bootstrap samples
NOISE_SCALE   = 0.1   # relative noise on f
START_DATE    = '2022-2-8'  # replace with actual start date

# Ensure directory
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- Load & smooth observed data ---
df_raw    = pd.read_excel(RAW_FILE)
daily_raw = df_raw.iloc[:,1].to_numpy().astype(float)
w = SMOOTH_WINDOW; pad = w//2
padded = np.pad(daily_raw, (pad,pad), 'constant')
kernel = np.ones(w)
smoothed = np.convolve(padded, kernel, mode='valid')
counter = np.convolve(np.pad(np.ones_like(daily_raw),(pad,pad)), kernel, mode='valid')
daily_obs = smoothed / counter

# --- Load PINN results & params ---
ts     = pd.read_excel(os.path.join(RESULTS_DIR, TS_FILE))
params = pd.read_excel(os.path.join(RESULTS_DIR, 'parameters.xlsx'))

t    = ts['Day'].to_numpy()
Se   = ts['Se'].to_numpy()
S    = ts['S'].to_numpy()
I0   = ts['I'].to_numpy()
R    = ts['R'].to_numpy()
f    = ts['f'].to_numpy()
beta  = float(params.loc[params['parameter']=='beta','value'].iloc[0])
gamma = float(params.loc[params['parameter']=='gamma','value'].iloc[0])
y0    = [Se[0], S[0], I0[0], R[0]]

# --- Bootstrap ODE ensemble ---
I_ode_samples = []
f_samples     = []
for _ in range(M_SAMPLES):
    noise = np.random.normal(0, NOISE_SCALE, size=f.shape)
    f_s   = f * (1 + noise)
    f_samples.append(f_s)
    def seir(y, ti):
        Se0,S0,I0_v,R0 = y
        ft = np.interp(ti, t, f_s)
        return [
            -ft,
            ft - beta*S0*I0_v/(S0+I0_v+R0),
            beta*S0*I0_v/(S0+I0_v+R0)-gamma*I0_v,
            gamma*I0_v
        ]
    sol = odeint(seir, y0, t)
    I_ode = beta * sol[:,1] * sol[:,2] / (sol[:,1] + sol[:,2] + sol[:,3])
    I_ode_samples.append(I_ode)

I_arr     = np.stack(I_ode_samples)
median_I  = np.percentile(I_arr, 50, axis=0)
lower_I   = np.percentile(I_arr, 2.5, axis=0)
upper_I   = np.percentile(I_arr, 97.5, axis=0)
f_arr     = np.stack(f_samples)
median_f  = np.percentile(f_arr, 50, axis=0)
lower_f   = np.percentile(f_arr, 2.5, axis=0)
upper_f   = np.percentile(f_arr, 97.5, axis=0)

# --- Plotting (2 rows) ---
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(12,8), gridspec_kw={'hspace':0.15})

# Common style: larger tick labels
for ax in (ax1, ax2):
    ax.tick_params(axis='x', labelsize=18, width=2)

# Top: Daily infections
ax1.scatter(t, daily_obs, s=30, color='black', alpha=0.6)
ax1.plot(t, median_I, color='#1f77b4', lw=4)
ax1.fill_between(t, lower_I, upper_I, color='#1f77b4', alpha=0.3)
ax1.set_ylabel('Daily New Infections', color='#1f77b4', fontsize=20)
ax1.tick_params(axis='y', labelsize=18, colors='#1f77b4', width=2)
ax1.set_title('Daily New Infections with 95% CI', loc='left', fontsize=22)
# Minimalist axes
for spine in ['top','right','bottom','left']:
    ax1.spines[spine].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.grid(alpha=0.1)

# Bottom: Learned flow f(t) with updated color
flow_color = '#2ca02c'
ax2.plot(t, median_f, color=flow_color, lw=4)
ax2.fill_between(t, lower_f, upper_f, color=flow_color, alpha=0.3)
ax2.set_ylabel('Learned Flow f(t)', color=flow_color, fontsize=20)
ax2.tick_params(axis='y', labelsize=18, colors=flow_color, width=2)
ax2.set_xlabel('Day', fontsize=20)
ax2.set_title('Learned Flow f(t) with 95% CI', loc='left', fontsize=22)
for spine in ['top','right','bottom','left']:
    ax2.spines[spine].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.grid(alpha=0.1)

# --- Add start date textbox under x=0 ---
x0 = t[0]
y_min = ax2.get_ylim()[0]
ax2.text(x0, y_min - 0.05*(ax2.get_ylim()[1]-y_min), START_DATE,
         ha='center', va='top', fontsize=16,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=1.5))

# Adjust layout manually
fig.subplots_adjust(top=0.95, bottom=0.1, left=0.12, right=0.95)

# Save figure
fig_path = os.path.join(RESULTS_DIR, FIG_FILE)
fig.savefig(fig_path, dpi=300)
print(f"Saved figure to {fig_path}")
