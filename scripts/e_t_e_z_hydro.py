import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Estilo global
# -----------------------------
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 18,
    "axes.titlesize": 16,
    "legend.fontsize": 13,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "lines.linewidth": 2,
    "lines.markersize": 3,
    "figure.figsize": (9,10)
})

filename = "fort.10"

nx, ny, nz = 150, 150, 452
ix0, iy0, kz0 = nx//2, ny//2, nz//2

dt = 0.006671282

regions = {
    "center": (ix0, iy0, kz0),
    "mid":    (ix0 + 10, iy0, kz0),
    "edge":   (ix0 + 20, iy0, kz0),
}

data = {
    name: {"t": [], "e": [], "vz": [], "T0z": []}
    for name in regions
}

current_cycle = None
reading = False

# -----------------------------
# Lectura
# -----------------------------
with open(filename, "r") as f:
    for line in f:

        if "# Time(ncyc)=" in line:
            parts = line.split()
            try:
                current_cycle = int(parts[2])
            except:
                continue
            reading = True
            continue

        if line.startswith("#"):
            continue

        if not reading:
            continue

        parts = line.split()
        if len(parts) < 10:
            continue

        i = int(parts[0])

        if i == 9999:
            reading = False
            continue

        j = int(parts[1])
        k = int(parts[2])

        T00 = float(parts[3])
        vz  = float(parts[6])
        T0z = float(parts[9])

        for name, (ix, iy, kz) in regions.items():
            if i == ix and j == iy and k == kz:
                time = current_cycle * dt

                data[name]["t"].append(time)
                data[name]["e"].append(T00)
                data[name]["vz"].append(vz)
                data[name]["T0z"].append(T0z)

# -----------------------------
# Plot
# -----------------------------
fig, axes = plt.subplots(3, 1, sharex=True)

colors = {
    "center": "tab:blue",
    "mid": "tab:orange",
    "edge": "tab:green"
}

for name in regions:
    t = np.array(data[name]["t"])

    idx = np.argsort(t)
    t = t[idx]

    e   = np.array(data[name]["e"])[idx]
    vz  = np.array(data[name]["vz"])[idx]
    T0z = np.array(data[name]["T0z"])[idx]

    axes[0].plot(t, e,  '-o', color=colors[name], label=name)
    axes[1].plot(t, vz, '-o', color=colors[name])
    axes[2].plot(t, T0z,'-o', color=colors[name])

# -----------------------------
# Etiquetas
# -----------------------------
axes[0].set_ylabel(r"$\varepsilon\ \mathrm{[MeV/fm^3]}$")
axes[1].set_ylabel(r"$v_z$")
axes[2].set_ylabel(r"$T^{0z}\ \mathrm{[MeV/fm^3]}$")
axes[2].set_xlabel(r"$t\ \mathrm{[fm}/c]$")

# Grid
for ax in axes:
    ax.grid(True, linestyle="--", alpha=0.5)

# Leyenda solo arriba
axes[0].legend(frameon=False)

plt.tight_layout()
plt.savefig("hydro_fields_vs_time.pdf", dpi=300)
plt.show()
