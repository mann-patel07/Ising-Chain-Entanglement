"""Free fermion solution for periodic chains (even parity sector). Added after the first version
(see UPGRADE_PLAN.md). Writes results/ff_periodic_*.csv, results/ed_vs_ff_periodic.csv and
results/central_charge_periodic.json."""
import csv
import json
import time

import numpy as np

from tfim_ed import ground_state
from entropy import entropy_profile
from free_fermion import majorana_matrix_periodic, ground_state_covariance, block_entropy
from central_charge import fit_c, c_from_two_sizes

H = 1.0
LS = [16, 32, 64, 128, 256, 512, 1024, 2048]
FIT_L = 1024
CROSS_LS = [8, 12, 16, 20]
CROSS_HS = (1.0, 1.5, 2.0)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([repr(float(x)) if isinstance(x, (float, np.floating)) else x for x in r])


def ff_solve(L, h):
    return ground_state_covariance(majorana_matrix_periodic(L, h))


def main():
    t0 = time.time()
    half, rows_half = {}, []
    fit_ls, fit_S, rows_prof = [], [], []
    for L in LS:
        G, E0 = ff_solve(L, H)
        s = block_entropy(G, L // 2)
        half[L] = s
        rows_half.append([L, H, E0, s])
        if L == FIT_L:
            fit_ls = list(range(64, 513, 16))
            for l in fit_ls:
                v = block_entropy(G, l)
                fit_S.append(v)
                rows_prof.append([L, H, l, v])
        print("L = %4d  S_half = %.10f" % (L, s), flush=True)
    write_csv("results/ff_periodic_half_chain.csv", ["L", "h", "E0", "S_half"], rows_half)
    write_csv("results/ff_periodic_profiles.csv", ["L", "h", "l", "S"], rows_prof)

    cross = []
    for L in CROSS_LS:
        for h in CROSS_HS:
            E_ed, psi = ground_state(L, h, "periodic")
            S_ed = entropy_profile(psi, L)
            G, E_ff = ff_solve(L, h)
            dS = max(abs(S_ed[l] - block_entropy(G, l)) for l in range(L + 1))
            cross.append([L, h, dS, abs(E_ed - E_ff)])
    write_csv("results/ed_vs_ff_periodic.csv", ["L", "h", "max_abs_dS", "abs_dE"], cross)

    pairs = [{"L1": a, "L2": b, "c": c_from_two_sizes(half[a], a, half[b], b, "periodic")}
             for a, b in zip(LS[:-1], LS[1:])]
    c_fit, _ = fit_c(fit_ls, fit_S, FIT_L, "periodic")
    out = {
        "ff_periodic_pairs": pairs,
        "c_periodic_final": pairs[-1]["c"],
        "c_periodic_systematic": abs(pairs[-1]["c"] - pairs[-2]["c"]),
        "c_fit_ff_periodic_L1024": c_fit,
        "max_abs_dS": max(r[2] for r in cross),
        "max_abs_dE": max(r[3] for r in cross),
        "run_time_seconds": time.time() - t0,
    }
    with open("results/central_charge_periodic.json", "w") as f:
        json.dump(out, f, indent=2)
    print("c_periodic_final = %.6f  systematic = %.2e  fit = %.6f" % (out["c_periodic_final"], out["c_periodic_systematic"], c_fit))
    print("max dS = %.2e  max dE = %.2e" % (out["max_abs_dS"], out["max_abs_dE"]))
    print("run time %.1f s" % out["run_time_seconds"])


if __name__ == "__main__":
    main()
