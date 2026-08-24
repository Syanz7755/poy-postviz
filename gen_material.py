"""
Generate 27-point silica material file from 53-point data.
Format: omega (rad/s)  Re(epsilon)  Im(epsilon)
"""
import numpy as np

# Read 53-point data
freqs = []  # Hz
eps_real = []
eps_img = []

with open("aligned-SiO2_Franta-300C-sparse_53 - real.txt") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            freqs.append(float(parts[0]))
            eps_real.append(float(parts[1]))

with open("aligned-SiO2_Franta-300C-sparse_53 - img.txt") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            eps_img.append(float(parts[1]))

freqs = np.array(freqs)
eps_real = np.array(eps_real)
eps_img = np.array(eps_img)

print(f"Read {len(freqs)} points")
print(f"Freq range: {freqs[0]:.3e} - {freqs[-1]:.3e} Hz")

# Select 27 points: capture important features
# Strategy: take every 2nd point, plus key peaks/valleys
# Total 53 points -> 27 points: roughly every 2nd point

indices = list(range(0, 53, 2))  # 0, 2, 4, ..., 52 -> 27 points
print(f"Selected {len(indices)} points at indices: {indices}")

# Convert frequency to angular frequency (rad/s)
omega = 2 * np.pi * freqs[indices]

# Write material file
# scuff-em format: omega  Re(epsilon)  Im(epsilon)
with open("materials/silica-27p.txt", "w") as f:
    f.write("# Silica (Franta 300C) - 27 points\n")
    f.write("# omega(rad/s)  Re(eps)  Im(eps)\n")
    for i, idx in enumerate(indices):
        f.write(f"{omega[i]:.6e}  {eps_real[idx]:.6e}  {eps_img[idx]:.6e}\n")

print(f"\nWrote materials/silica-27p.txt with 27 points")
print(f"Omega range: {omega[0]:.3e} - {omega[-1]:.3e} rad/s")
