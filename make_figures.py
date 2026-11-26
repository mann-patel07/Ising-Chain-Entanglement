"""Make the figures from results/. Reads only from results/."""
import csv
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from central_charge import chord_periodic, chord_open, fit_c

plt.rcParams.update({"font.size": 11})
R = "results"
os.makedirs("figures", exist_ok=True)


def read(name):
    with open(f"{R}/{name}") as f:
        return list(csv.DictReader(f))


def save(fig, name):
    fig.savefig(f"figures/{name}", dpi=150, bbox_inches="tight")
    plt.close(fig)


with open(f"{R}/central_charge.json") as f:
    cc = json.load(f)

prof = read("ed_profiles.csv")
Lp = int(prof[0]["L"])


def ed_S(boundary, h):
    rows = [r for r in prof if r["boundary"] == boundary and float(r["h"]) == h]
    rows.sort(key=lambda r: int(r["l"]))
    return np.array([int(r["l"]) for r in rows]), np.array([float(r["S"]) for r in rows])


# Fig 1
fig, ax = plt.subplots(figsize=(6, 4))
for h in (0.5, 1.0, 2.0):
    l, S = ed_S("periodic", h)
    ax.plot(l, S, "o-", ms=4, label="h = %.1f" % h)
ax.set_xlabel("block size l (sites)")
ax.set_ylabel("entanglement entropy S (dimensionless)")
ax.set_title("Periodic chain, L = %d, exact diagonalization" % Lp)
ax.legend()
save(fig, "fig1_entropy_profiles.png")

# Fig 2
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
l, S = ed_S("periodic", 1.0)
m = (l >= 2) & (l <= Lp - 2)
c, b = fit_c(l[m], S[m], Lp, "periodic")
x = chord_periodic(Lp, l[m])
ax = axes[0]
ax.plot(x, S[m], "o", label="exact diagonalization")
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, (c / 3) * xs + b, "-", label="fit, c = %.3f" % c)
ax.set_xlabel("ln[(L/pi) sin(pi l / L)] (dimensionless)")
ax.set_ylabel("entanglement entropy S (dimensionless)")
ax.set_title("Periodic, L = %d, h = 1" % Lp)
ax.legend()
ffp = read("ff_profiles.csv")
Lfit = int(cc["ff_open_pairs"][-2]["L2"])
pts = [(int(r["l"]), float(r["S"])) for r in ffp if int(r["L"]) == Lfit and float(r["h"]) == 1.0]
pts.sort()
ls = np.array([p[0] for p in pts])
Sf = np.array([p[1] for p in pts])
c2, b2 = fit_c(ls, Sf, Lfit, "open")
x = chord_open(Lfit, ls)
ax = axes[1]
ax.plot(x, Sf, "o", label="free fermions")
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, (c2 / 6) * xs + b2, "-", label="fit, c = %.3f" % c2)
ax.set_xlabel("ln[(2L/pi) sin(pi l / L)] (dimensionless)")
ax.set_ylabel("entanglement entropy S (dimensionless)")
ax.set_title("Open, L = %d, h = 1" % Lfit)
ax.legend()
save(fig, "fig2_calabrese_cardy.png")

# Fig 3
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot([1 / p["L2"] for p in cc["ff_open_pairs"]], [p["c"] for p in cc["ff_open_pairs"]], "o-", label="free fermion, open")
ax.plot([1 / p["L2"] for p in cc["ed_periodic_pairs"]], [p["c"] for p in cc["ed_periodic_pairs"]], "s-", label="ED, periodic")
ax.axhline(0.5, ls="--", color="k", label="c = 1/2")
ax.set_xscale("log")
ax.set_xlabel("1 / L2 (1/sites)")
ax.set_ylabel("c from two sizes (dimensionless)")
ax.legend()
save(fig, "fig3_central_charge_vs_size.png")

# Fig 4
scan = read("field_scan.csv")
fig, ax = plt.subplots(figsize=(6, 4))
for L in sorted({int(r["L"]) for r in scan}):
    rows = sorted([r for r in scan if int(r["L"]) == L], key=lambda r: float(r["h"]))
    ax.plot([float(r["h"]) for r in rows], [float(r["S_half"]) for r in rows], "o-", ms=4, label="L = %d" % L)
ax.axvline(1.0, ls=":", color="k")
ax.set_xlabel("transverse field h (units of J)")
ax.set_ylabel("half chain entropy S (dimensionless)")
ax.legend()
save(fig, "fig4_half_chain_vs_field.png")

# Fig 5
spec = read("entanglement_spectrum.csv")
fig, ax = plt.subplots(figsize=(6, 4))
for h in (0.5, 1.0, 2.0):
    rows = sorted([r for r in spec if r["boundary"] == "periodic" and float(r["h"]) == h], key=lambda r: int(r["rank"]))
    ax.plot([int(r["rank"]) for r in rows], [-np.log(float(r["lambda"])) for r in rows], "o-", ms=4, label="h = %.1f" % h)
ax.set_xlabel("rank (dimensionless)")
ax.set_ylabel("-ln(lambda) (dimensionless)")
ax.legend()
save(fig, "fig5_entanglement_spectrum.png")
print("figures written")
