import numpy as np
import os

nu = 1.5e-5                 # m^2/s, kin. viscosity, air
U_max_target = 10.5         # m/s, U_max
utau_guess = U_max_target / 17  # friction velocity, estimated
kappa = 0.41
B = 5.0


y_vals = np.linspace(0.0005, 0.0395, 9)     # y [0.5 mm – 39.5 mm]
z_vals = np.linspace(-0.2285, 0.2285, 11)   # z [–228.5 mm – 228.5 mm]
x_val = -0.4                                # x = -400 mm


os.makedirs("constant/boundaryData/inlet/0", exist_ok=True)


u_y = []
for y in y_vals:
    yPlus = (y * utau_guess) / nu
    if yPlus < 1:
        u = 0
    else:
        u = utau_guess * ((1 / kappa) * np.log(yPlus) + B)
    u_y.append(u)


u_max_actual = max(u_y)
scaling_factor = U_max_target / u_max_actual
print(f"Real max velocity: {u_max_actual:.4f}, Scaled: {scaling_factor:.4f}")


with open("constant/boundaryData/inlet/points", "w") as f_points:
    f_points.write("(\n")
    for z in z_vals:
        for y in y_vals:
            f_points.write(f"({x_val:.6f} {y:.6f} {z:.6f})\n")
    f_points.write(")\n")


with open("constant/boundaryData/inlet/0/U.air", "w") as f_U:
    f_U.write("(\n")
    for z in z_vals:
        for i, y in enumerate(y_vals):
            u_scaled = u_y[i] * scaling_factor
            f_U.write(f"({u_scaled:.6f} 0 0)\n")
    f_U.write(")\n")

