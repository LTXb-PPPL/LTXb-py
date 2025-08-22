"""
This code processes data from a tungsten wire calorimeter experiment to see a map of where the beam hits.
You can use the slider at the bottom to select which wire to visualize (0 = wire 1, 1 = wire 2, etc.).
"""

#import necessary libraries
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.optimize import curve_fit

file_path = Path("/Users/husky/Desktop/wcal/111215 0.6A.txt")  

#Settings
WIRE_CURRENT = 0.6     # Current
WIN_R1 = 200           # samples before pulse for baseline
WIN_R2 = 200           # samples for plateau averaging

#TUNGSTEN WIRE PROPERTIES & FORMULAS
TCR = 0.0045           # tungsten temperature coefficient of resistance
TUNGSTEN_DENSITY = 19300       # kg/m³
SPECIFIC_HEAT    = 134         # J/(kg·K)
WIRE_LENGTH      = 0.2286
WIRE_DIAMETER    = 0.000254
WIRE_RADIUS      = WIRE_DIAMETER / 2.0
WIRE_AREA        = np.pi * WIRE_RADIUS**2 
WIRE_VOLUME      = WIRE_AREA * WIRE_LENGTH
WIRE_MASS       = TUNGSTEN_DENSITY * WIRE_VOLUME    # kg

#Interpreting The Files 
def find_data_start(path: Path) -> int:
    #finding where data starts in the text file
    with open(path, "r", errors="ignore") as f:
        end = -1
        for i, line in enumerate(f):
            if line.strip() == "***End_of_Header***":
                end = i
    if end < 0:
        raise RuntimeError(f"No ***End_of_Header*** found in {path}")
    return end + 1

def load_lvm(path: Path) -> pd.DataFrame:
    start_row = find_data_start(path)
    data= pd.read_csv(path, sep="\t", skiprows=start_row, engine="python")
    data= data.dropna(axis=1, how='all')
    data["Time (s)"] = data.index * 2e-1  # assume 0.2s per sample
    return data

#Pulse Detection and Plateau Detection
def find_pulse_and_plateau(signal, win=10): #win 10 is how many data points to take
    dV = np.diff(signal) # slope at each point
    pulse_idx = np.argmax(np.abs(dV)) # finds where largest change occurs
    search = dV[pulse_idx+win:]
    plateau_rel = np.argmin(np.abs(search))
    plateau_idx = pulse_idx + win + plateau_rel
    return pulse_idx, plateau_idx

#Table building
def build_power_table(data: pd.DataFrame) -> pd.DataFrame:
    channel_cols = [c for c in data.columns if c not in ["Time (s)", "Comment"]]
    rows = []
    pulse_idx, plateau_idx = find_pulse_and_plateau(data[channel_cols[0]].values)
    dt_sample = data["Time (s)"].iloc[1] - data["Time (s)"].iloc[0]
    delta_t = 0.005
    print(f"Pulse at ~{data['Time (s)'].iloc[pulse_idx]:.4f}s, Δt={delta_t:.4f}s")

    for col in channel_cols:
        voltage = data[col].values
        resistance = voltage / WIRE_CURRENT  # Ohm's law

        # baseline before pulse
        R_1 = np.mean(resistance[max(0, pulse_idx - WIN_R1):pulse_idx])
        # plateau after pulse
        R_2 = np.mean(resistance[plateau_idx:plateau_idx + WIN_R2])

        # change in resistance → temperature rise
        delta_R = R_2 - R_1
        delta_T = delta_R / (R_1 * TCR) if R_1 != 0 else float("something wrong in line 87")

        # beam power per unit cross-sectional area
        P_per_area = TUNGSTEN_DENSITY * WIRE_LENGTH * SPECIFIC_HEAT * (delta_T / delta_t) if delta_t > 0 else float("line 90 issue")

        #adding rows to the table
        rows.append({
            "Channel": col,
            "R_1 (Ω)": R_1,
            "R_2 (Ω)": R_2,
            "ΔR (Ω)": delta_R,
            "ΔT (°C)": delta_T,
            "P_per_area (W/m²)": P_per_area
        })
    return pd.DataFrame(rows)

#GAUSSIAN MODEL FOR FITTING
def gaussian(x, A, mu, sigma, offset):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2)) + offset

#creating plot of data and overlaying the gaussian fit
def fit_and_plot(ax, x, y, label, color="blue"):
    ax.plot(x, y, "o", label=f"{label} data", color=color)
    p0 = [max(y)-min(y), np.mean(x), len(x)/4, min(y)]
    popt, _ = curve_fit(gaussian, x, y, p0=p0)
    xfit = np.linspace(min(x), max(x), 200)
    ax.plot(xfit, gaussian(xfit, *popt), "-", color=color, alpha=0.7,
            label=f"{label} fit (μ={popt[1]:.2f}, σ={popt[2]:.2f})")
    ax.legend()

# ===== DATA PROCESSING =====
data = load_lvm(file_path)
shot_name = file_path.stem
print(f"\n===== SHOT NAME: {shot_name} =====\n")

channel_cols = [c for c in data.columns if c not in ["Time (s)", "Comment"]]
table = build_power_table(data)

# remove channel 16 (if present)
table = table[~table["Channel"].str.contains("16", case=False, na=False)].reset_index(drop=True)

# ===== PLOT WITH SLIDER =====
time = data["Time (s)"].values

fig, ax = plt.subplots(2,1, figsize=(10,6))
plt.subplots_adjust(bottom=0.25)

line_deltaT, = ax[0].plot([], [], "b-", label="ΔT(t)")
line_power, = ax[1].plot([], [], "orange", label="P_per_area(t)")

ax[0].set_ylabel("ΔT (°C)")
ax[0].set_xlabel("Time (s)")
ax[0].grid(True)
ax[0].set_title(f"{shot_name} - Channel ΔT(t)")

ax[1].set_ylabel("P_per_area (W/m²)")
ax[1].set_xlabel("Time (s)")
ax[1].grid(True)
ax[1].set_title(f"{shot_name} - Channel Power per Unit Area")

# Slider axis
ax_slider = plt.axes([0.25, 0.1, 0.5, 0.03])
slider = Slider(ax_slider, 'Channel', 0, len(channel_cols)-1, valinit=0, valstep=1)

def update(val):
    idx = int(slider.val)
    col = channel_cols[idx]
    voltage = data[col].values
    resistance = voltage / WIRE_CURRENT
    pulse_idx, plateau_idx = find_pulse_and_plateau(voltage)
    R_1 = np.mean(resistance[max(0, pulse_idx - WIN_R1):pulse_idx])
    delta_T = (resistance - R_1) / (R_1 * TCR)
    dt_sample = data["Time (s)"].iloc[1] - data["Time (s)"].iloc[0]
    P_per_area = TUNGSTEN_DENSITY * WIRE_LENGTH * SPECIFIC_HEAT * (delta_T / dt_sample)
    
    line_deltaT.set_data(time, delta_T)
    line_power.set_data(time, P_per_area)
    
    ax[0].relim()
    ax[0].autoscale_view()
    ax[1].relim()
    ax[1].autoscale_view()
    
    ax[0].set_title(f"{shot_name} - {col} ΔT(t)")
    ax[1].set_title(f"{shot_name} - {col} Power per Unit Area")
    fig.canvas.draw_idle()

slider.on_changed(update)
update(0)  # initialize with first channel

plt.show()
