#library imports
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go

folder_path = Path("/Users/husky/Desktop/wcal/") 
file_list = sorted(folder_path.glob("*.txt"))

#data sorting
def find_data_start(path: Path) -> int:
    """Find the line after the last ***End_of_Header*** marker."""
    with open(path, "r", errors="ignore") as f:
        end = -1
        for i, line in enumerate(f):
            if line.strip() == "***End_of_Header***":
                end = i
    if end < 0:
        raise RuntimeError(f"No ***End_of_Header*** found in {path}")
    return end + 1

def load_lvm(path: Path):
    """Load LVM file and create time axis."""
    start_row = find_data_start(path)
    df = pd.read_csv(path, sep="\t", skiprows=start_row, engine="python")
    df = df.dropna(axis=1, how='all')
    df["Time (s)"] = df.index * 1e-3  # adjust if Delta_X is different
    return df

# Preload data
datasets = {f.name: load_lvm(f) for f in file_list}

# Build figure for the first file only
first_file = list(datasets.keys())[0]
df = datasets[first_file]
time_col = "Time (s)"
channel_cols = [c for c in df.columns if c not in [time_col, "Comment"]]

fig = go.Figure()

for col in channel_cols:
    fig.add_trace(go.Scatter(
        x=df[time_col],
        y=df[col],
        mode='lines',
        name=col
    ))

# Dropdown menu to reload traces
buttons = []
for fname in datasets.keys():
    buttons.append(dict(
        label=fname,
        method="update",
        args=[{"x": [datasets[fname][time_col]] * len(channel_cols),
               "y": [datasets[fname][col] for col in channel_cols],
               "name": channel_cols},
              {"title": f"Wire Calorimeter - {fname}"}]
    ))

 
fig.update_layout(
    title=f"Wire Calorimeter - {first_file}",
    xaxis_title="Time (s)",
    yaxis_title="Voltage (V)",
    updatemenus=[dict(
        buttons=buttons,
        direction="down",
        x=1.05,
        y=1.15
    )]
)

#UPDATE THE Y AXES INCASE U WANT TO SEE CLOSER DATA BUT CAN USE THE INTEARCTIVE GUI AS WELL
fig.update_yaxes(range=[-1.75, 1.75])

fig.show()
