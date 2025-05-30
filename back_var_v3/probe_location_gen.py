import numpy as np
import matplotlib.pyplot as plt

x_min, x_max, x_step = -0.35, 1.35, 0.05
z_val = 0.0
y_max = 0.035
dy = 0.003

x_vals = np.arange(x_min, x_max + x_step, x_step)

points = []

for x in x_vals:
    y_min = 0.0 if x < 0 else -0.025
    y_vals = np.arange(y_min, y_max + dy, dy)

    for y in y_vals:
        points.append((x, y, z_val))

x_coords = [p[0] for p in points]
y_coords = [p[1] for p in points]

with open("system/probeLocations", "w") as f:
    f.write("probeLocations\n(\n")
    for x, y, z in points:
        f.write(f"    ({x:.6f} {y:.6f} {z:.6f})\n")
    f.write(");\n")

plt.figure(figsize=(12, 5))
plt.scatter(x_coords, y_coords, s=12, color='black', marker='o')
plt.title("Probe Locations", fontsize=14)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid(True, linestyle='--', alpha=0.5)
plt.axis('equal')

info_text = f"dx = {x_step} m\n" + f"dy = {dy} m"
plt.text(x_max - 0.45, y_max + 0.005, info_text, fontsize=12, color='blue')


plt.tight_layout()
plt.show()
