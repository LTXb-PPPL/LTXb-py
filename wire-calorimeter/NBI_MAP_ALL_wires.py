"""This code processes data from a tungsten wire calorimeter experiment to see a map of where the beam hits.
we are looking at the heat distribution across the wires.
It calculates the temperature rise in each wire based on resistance changes, and computes the power per cross sectional unit area.
It also detects pulses in the data and averages the resistance before and after the pulse to determine the temperature change.
The results are stored in a structured table format, and Gaussian fits are applied to the data for visualization.
The script assumes a specific tungsten wire setup and uses constants for density and specific heat.
It also plots the temperature profiles and beam intensity profiles for both horizontal and vertical wires.
The Gaussian fitting is used to model the data and provide a better understanding of the distribution of heat
across the wires."""

#import necessary libraries
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
   
   #Detect pulse (max slope) and steady plateau after it
    dV = np.diff(signal) #slope at each point
    pulse_idx = np.argmax(np.abs(dV)) #finds where largest change occurs
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

# data processing
data= load_lvm(file_path)
table = build_power_table(data)
shot_name = file_path.stem
print(f"\n===== SHOT NAME: {shot_name} =====\n")
# remove channel 16 (if present)
table = table[~table["Channel"].str.contains("16", case=False, na=False)].reset_index(drop=True)

print(table.to_string(index=False))

# rows 0–7 = vertical, rows 8–15 = horizontal
vert  = table.iloc[0:8].reset_index(drop=True)
horiz = table.iloc[8:16].reset_index(drop=True)

def plot_channel_with_pulse(data: pd.DataFrame, channel: str, table: pd.DataFrame):
    """Plot one channel's time series with markers at pulse and plateau indices."""

    time = data["Time (s)"].values
    voltage = data[channel].values / WIRE_CURRENT  # resistance (Ω)

    # grab row for this channel from the table
    row = table[table["Channel"] == channel].iloc[0]
    pulse_idx = int(row["Pulse index"])
    plateau_idx = int(row["Plateau index"])

    plt.figure(figsize=(10, 5))
    plt.plot(time, voltage, label=f"{channel} resistance", color="blue")

    # mark pulse and plateau points
    plt.axvline(time[pulse_idx], color="red", linestyle="--", label=f"Pulse @ {time[pulse_idx]:.4f}s")
    plt.axvline(time[plateau_idx], color="green", linestyle="--", label=f"Plateau @ {time[plateau_idx]:.4f}s")

    plt.xlabel("Time (s)")
    plt.ylabel("Resistance (Ω)")
    plt.title(f"Channel {channel}: pulse/plateau detection")
    plt.legend()
    plt.grid(True)
    plt.show()
#PLOTTING
fig, axes = plt.subplots(2, 2, figsize=(12,8))

# Horizontal ΔT
axes[0,0].set_title(f"{shot_name} - Horizontal Wires ΔT Profile")
fit_and_plot(axes[0,0], horiz.index, horiz["ΔT (°C)"], "Horizontal Wires ΔT", "blue")
axes[0,0].set_xlabel("Wire index (0–7)")
axes[0,0].set_ylabel("ΔT (°C)")
axes[0,0].grid(True)

# Vertical ΔT
fit_and_plot(axes[0,1], vert.index, vert["ΔT (°C)"], "Vertical Wires ΔT", "red")
axes[0,1].set_title("Vertical ΔT Wires Profile")
axes[0,1].set_xlabel("Wire index (8–15)")
axes[0,1].grid(True)

# Horizontal beam intensity (P_per_area)
fit_and_plot(axes[1,0], horiz.index, horiz["P_per_area (W/m²)"], "Horizontal Wires Beam Intensity", "blue")
axes[1,0].set_title("Horizontal Wires Beam Intensity Profile")
axes[1,0].set_xlabel("Wire index (0–7)")
axes[1,0].set_ylabel("Power per unit area (W/m²)")
axes[1,0].grid(True)

# Vertical beam intensity (P_per_area)
fit_and_plot(axes[1,1], vert.index, vert["P_per_area (W/m²)"], "Vertical Wires Beam Intensity", "red")
axes[1,1].set_title("Vertical Wires Beam Intensity Profile")
axes[1,1].set_xlabel("Wire index (8–15)")
axes[1,1].set_ylabel("Power per unit area (W/m²)")
axes[1,1].grid(True)

plt.tight_layout()
plt.show()
