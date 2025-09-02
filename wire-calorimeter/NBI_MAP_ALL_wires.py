#nbi map
#hussain gajani

#imports
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import plotly.io as pio


#FILE SELECTION
file_list = [Path("/Users/husky/Desktop/wcal/111525.txt")]

#FINDING WHERE TO START DATA
def find_data_start(path: Path) -> int:
    with open(path, "r", errors="ignore") as f:
        end = -1
        for i, line in enumerate(f):
            if line.strip() == "***End_of_Header***":
                end = i
    return end + 1

#LOADING LVM TO TEXT FILE IN A DATA FRAME
def load_lvm(path: Path):
    start_row = find_data_start(path)
    df = pd.read_csv(path, sep="\t", skiprows=start_row, engine="python")
    df = df.dropna(axis=1, how="all").reset_index(drop=True)
    df["Time (s)"] = df.index * 1e-3  # adjust Δt if needed OR 0.002
    for c in df.columns:
        if c not in ["Time (s)", "Comment"]:
            df[c] = pd.to_numeric(df[c]) # CONVERT data to numeric
    return df

#PULSE AND PLATEAU DETECTION PARAMTERS
TIME_WINDOW = (1.25, 1.35)  # seconds
EPS_DERIV = 0.01            # how flat slope must be to call plateau
MIN_DUR = 0.015             # plateau duration
DERIV_FLOOR = 0.05          # minimum negative slope to count as pulse

#slope helper calcualting slope using moving window size (5), giving rate of change at each point in time series
def rolling_slope(t, y, win=5):
    t = np.asarray(t, float); y = np.asarray(y, float)
    if len(t) <= win:
        return np.full_like(t, np.nan)
    dt = t[win:] - t[:-win]
    m = (y[win:] - y[:-win]) / np.where(dt == 0, np.nan, dt)
    return np.r_[m, np.full(win, np.nan)]

# FINDING PULSE AND PLTEAU
def pick_pulse_plateau_gate(t_ser, y_ser, win=5):
# Pulse = last near-zero slope before strong negative drop inside TIME_WINDOW.
# Plateau = after dip, first place slope>0 then flattens APPROX close to 0 for MIN_DUR.
# Plateau marker = 5th sample after that flattening.
    t = np.asarray(t_ser, float)
    y = np.asarray(y_ser, float)

    # restrict to time window
    mask = (t >= TIME_WINDOW[0]) & (t <= TIME_WINDOW[1])
    if not np.any(mask):
        return 0, 0
    idx0, idx1 = np.where(mask)[0][[0, -1]]

    m = rolling_slope(t, y, win=win)

    #finding pulse_idx (intital pulse point)
    #pick the first strong negative slope as the pulse and then search 
    #for a positive slope and flattening for the plateau regardless of the order of slopes in the raw data
    pulse_idx = None
    for i in range(idx0, idx1):
        if m[i] < -DERIV_FLOOR:
            # look backwards to last near-zero slope
            back = np.where(np.abs(m[idx0:i]) < EPS_DERIV)[0]
            if len(back):
                pulse_idx = back[-1]
            else:
                pulse_idx = i
            break

    if pulse_idx is None:
        # fallback = strongest negative slope inside the window
        rel_idx = np.nanargmin(m[idx0:idx1+1])
        pulse_idx = idx0 + int(rel_idx)

    #making sure pulse idx is in the window we want
    pulse_idx = max(idx0, min(pulse_idx, idx1))

    #finding plateau_idx (point where slope flattens after pulse)
    plateau_idx = idx1
    found = False
    for i in range(pulse_idx+1, idx1-win):
        # slope > 0 and flattens for MIN_DUR
        if m[i] > 0:
            dur = int(MIN_DUR / (t[1]-t[0]))
            if i+dur < len(m) and np.all(np.abs(m[i:i+dur]) < EPS_DERIV):
                plateau_idx = min(i+5, len(t)-1)
                found = True
                break
    if not found:
        plateau_idx = idx1

    return int(pulse_idx), int(plateau_idx)



#load data set
datasets = {f.name: load_lvm(f) for f in file_list}

#plot points from first file for a pulse idx and plateau idx plot to visualize the points
first_file = list(datasets.keys())[0]
df = datasets[first_file]
time_col = "Time (s)"
channel_cols = [c for c in df.columns if c not in [time_col, "Comment"]]

fig = go.Figure()

for col in channel_cols:
    y = df[col]
    pulse_idx, plateau_idx = pick_pulse_plateau_gate(df[time_col], y)
    fig.add_trace(go.Scatter(x=df[time_col], y=y, mode="lines", name=col))
    fig.add_trace(go.Scatter(x=[df[time_col].iloc[pulse_idx]], y=[y.iloc[pulse_idx]],
                             mode="markers", marker=dict(color="green", size=9, symbol="circle"),
                             name=f"{col} pulse"))
    fig.add_trace(go.Scatter(x=[df[time_col].iloc[plateau_idx]], y=[y.iloc[plateau_idx]],
                             mode="markers", marker=dict(color="red", size=9, symbol="diamond"),
                             name=f"{col} plateau"))

#drop down gui menu to pick different files (have to reopen the code and assign the file to see the map)
buttons = []
for fname, df in datasets.items():
    traces_x, traces_y, traces_names = [], [], []
    for col in channel_cols:
        y = df[col]
        pulse_idx, plateau_idx = pick_pulse_plateau_gate(df[time_col], y)

        traces_x.extend([
            df[time_col],
            [df[time_col].iloc[pulse_idx]],
            [df[time_col].iloc[plateau_idx]]
        ])
        traces_y.extend([
            y,
            [y.iloc[pulse_idx]],
            [y.iloc[plateau_idx]]
        ])
        traces_names.extend([
            f"{col} trace",
            f"{col} pulse",
            f"{col} plateau"
        ])
    buttons.append(dict(
        label=fname,
        method="update",
        args=[{"x": traces_x, "y": traces_y, "name": traces_names},
              {"title": f"Wire Calorimeter – {fname}"}]
    ))

fig.update_layout(
    title=f"Wire Calorimeter – {first_file}",
    xaxis_title="Time (s)",
    yaxis_title="Voltage (V)",
    updatemenus=[dict(
        buttons=buttons,
        direction="down",
        x=1.05,
        y=1.15
    )]
)

fig.update_yaxes(range=[-1.75, 1.75])
fig.show()

# tungsten wire properties
WIRE_CURRENT = 0.6
TCR = 0.0045

#building a table to manually confirm data
def build_power_table(df):
    channel_cols = [c for c in df.columns if c not in ["Time (s)", "Comment"]]
    rows = []

    for col in channel_cols:
        y = df[col].values
        pulse_idx, plateau_idx = pick_pulse_plateau_gate(df[time_col], y)

        resistance = y / WIRE_CURRENT  # Ohm’s law

        # baseline before pulse
        R_1 = np.mean(resistance[max(0, pulse_idx-5):pulse_idx])
        # plateau after pulse
        R_2 = np.mean(resistance[plateau_idx:plateau_idx+5])

        # ΔR and ΔT
        delta_R = R_2 - R_1
        delta_T = delta_R / (R_1 * TCR) if R_1 != 0 else np.nan
        
        
        #FLIPPING SIGNS UCOMMENT BELOW
        # if delta_R < 0:
        #     delta_R = abs(delta_R)
        #     delta_T = abs(delta_T)
        #     note = "Sign flipped"
        # else:
        #     note = "OK"

        rows.append({
            "Channel": col,
            "Pulse_idx": pulse_idx,
            "Plateau_idx": plateau_idx,
            "R_1 (Ω)": R_1,
            "R_2 (Ω)": R_2,
            "ΔR (Ω)": delta_R,
            "ΔT (°C)": delta_T
        })

    return pd.DataFrame(rows)

#print table
table = build_power_table(df)
print("\n=== Power Table ===")
print(table.to_string(index=False))

# split vertical (0–7) vs horizontal (8–15)
vert  = table.iloc[0:8].reset_index(drop=True)
horiz = table.iloc[8:16].reset_index(drop=True)



#plotting the ΔR and ΔT profiles NO GAUSSIAN FIT

fig = plt.plot(figsize=(12, 6))


# ΔT profile
plt.plot(horiz.index, horiz["ΔT (°C)"], "s-b", label="Horizontal ΔT")
plt.plot(vert.index,  vert["ΔT (°C)"], "s-r", label="Vertical ΔT")
plt.title("ΔT profile")
plt.xlabel("Wire index")
plt.ylabel("ΔT (°C)")
plt.legend()
plt.grid(True)

plt.tight_layout()





#GASSIAN FITTING
def gaussian(x, A, mu, sigma, offset):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2)) + offset

# Fit and plot helper
def fit_and_plot(ax, x, y, label, color, marker):
    # Scatter the data
    ax.plot(x, y, marker, label=f"{label} data", color=color)

    # Initial guess: amplitude, mean, sigma, offset
    p0 = [max(y) - min(y), np.mean(x), len(x)/4, np.min(y)]
    try:
        popt, _ = curve_fit(gaussian, x, y, p0=p0)
        xfit = np.linspace(min(x), max(x), 200)
        ax.plot(xfit, gaussian(xfit, *popt), "-", color=color,
                label=f"{label} fit (μ={popt[1]:.2f}, σ={popt[2]:.2f})")
    except RuntimeError:
        print(f"Fit failed for {label}")
    ax.legend()


fig = plt.plot(figsize=(12, 6))


# ΔT profile
plt.plot(horiz.index, horiz["ΔT (°C)"], "s-b", label="Horizontal ΔT")
plt.plot(vert.index,  vert["ΔT (°C)"], "s-r", label="Vertical ΔT")
plt.title("ΔT profile")
plt.xlabel("Wire index")
plt.ylabel("ΔT (°C)")
plt.legend()
plt.grid(True)

plt.tight_layout()

# Plotting ΔR and ΔT profiles with Gaussian fits 
fig = plt.plot(figsize=(12, 6))

# ΔT overlay
fit_and_plot(horiz.index, horiz["ΔT (°C)"], "Horizontal ΔT", "blue", "s")
fit_and_plot(vert.index,  vert["ΔT (°C)"],  "Vertical ΔT",   "red",  "s")
plt.set_title("ΔT profile")
plt.set_xlabel("Wire index")
plt.set_ylabel("ΔT (°C)")
plt.grid(True)

plt.tight_layout()
plt.show()

