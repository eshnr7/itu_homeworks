import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.ndimage import generic_filter


file_path = "/postProcessing/probes/0/U.air"
with open(file_path, 'r') as f:
    lines = f.readlines()
coords = []
for line in lines:
    if line.startswith("# Probe"):
        match = re.search(r'\(([^)]+)\)', line)
        if match:
            x, y, z = map(float, match.group(1).split())
            coords.append((x, y))

x_coords = [c[0] for c in coords]
y_coords = [c[1] for c in coords]
x_unique = sorted(set(x_coords))
y_unique = sorted(set(y_coords))
nx = len(x_unique)
ny = len(y_unique)

data_lines = [line.strip() for line in lines if not line.startswith("#") and line.strip()]

time_values = []
for line in data_lines:
    match = re.match(r'^\s*([\d.eE+-]+)', line)
    if match:
        time_values.append(float(match.group()))
    else:
        time_values.append(None)

target_time = 0.1
time_array = np.array([t for t in time_values if t is not None])
valid_indices = [i for i, t in enumerate(time_values) if t is not None]
index_t = valid_indices[np.argmin(np.abs(time_array - target_time))]
actual_time = time_values[index_t]


line_t = data_lines[index_t]
vectors = re.findall(r'\(([^)]+)\)', line_t)
velocity_vectors = [list(map(float, v.split())) for v in vectors]
U_magnitude = [np.linalg.norm(v) for v in velocity_vectors]


U_grid = np.full((ny, nx), np.nan)
for (x, y), u in zip(coords, U_magnitude):
    ix = x_unique.index(x)
    iy = y_unique.index(y)
    U_grid[iy, ix] = u


def nanmean_filter(values):
    valid = values[~np.isnan(values)]
    return np.mean(valid) if len(valid) > 0 else np.nan

U_grid_filled = U_grid.copy()
U_grid_filled = generic_filter(U_grid_filled, nanmean_filter, size=3, mode='constant', cval=np.nan)


plt.figure(figsize=(8, 6))
img = plt.imshow(U_grid_filled, origin='lower',
                 extent=[min(x_unique), max(x_unique), min(y_unique), max(y_unique)],
                 cmap='coolwarm', aspect='equal')
plt.colorbar(img, label='|U| [m/s] at t ≈ {:.3f} s'.format(actual_time))
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title(f"Velocity Contour at t ≈ {actual_time:.3f} s")
plt.grid(False)
plt.tight_layout()
plt.show()
