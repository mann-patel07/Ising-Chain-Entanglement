"""Run every calculation and write the numbers to results/ (or results_quick/ with --quick)."""
import argparse
import csv
import json
import os
import platform
import time

import numpy as np
import scipy
import matplotlib

from tfim_ed import ground_state, sector_energies
from entropy import entanglement_entropy, entropy_profile, schmidt_values
from free_fermion import majorana_matrix, ground_state_covariance, block_entropy
from central_charge import fit_c, c_from_two_sizes

BOUNDARIES = ("periodic", "open")


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([repr(float(x)) if isinstance(x, (float, np.floating)) else x for x in r])


def ff_solve(L, h):
    return ground_state_covariance(majorana_matrix(L, h))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    q = args.quick
    out = "results_quick" if q else "results"
    os.makedirs(out, exist_ok=True)

    if q:
        L_prof, L_en_max, L_half_max = 12, 12, 12
        ff_Ls = [16, 32, 64, 128, 256]
        L_ffprof = 64
        cross_Ls = [8, 12]
        scan_Ls = [8, 12]
        step = 0.3
    else:
        L_prof, L_en_max, L_half_max = 20, 20, 20
        ff_Ls = [16, 32, 64, 128, 256, 512, 1024, 2048]
        L_ffprof = 256
        cross_Ls = [8, 12, 16, 20]
        scan_Ls = [8, 12, 16, 20]
        step = 0.1
    fields = [0.5, 1.0, 2.0]
    fit_L = ff_Ls[-2]  # 1024 in the full run
    times = {}

    def timed(name):
        class T:
            def __enter__(s):
                s.t = time.time()

            def __exit__(s, *a):
                times[name] = time.time() - s.t
                print("%-28s %8.1f s" % (name, times[name]), flush=True)

        return T()

    # 1. ED profiles
    prof = {}
    half_profile = {}
    with timed("1_ed_profiles"):
        rows, srows = [], []
        for b in BOUNDARIES:
            for h in fields:
                _, psi = ground_state(L_prof, h, b)
                S = entropy_profile(psi, L_prof)
                prof[(b, h)] = S
                half_profile["%s_L%d_h%s" % (b, L_prof, h)] = float(S[L_prof // 2])
                for l in range(L_prof + 1):
                    rows.append([b, L_prof, h, l, S[l]])
                lam = schmidt_values(psi, L_prof, L_prof // 2)[:16]
                for r, v in enumerate(lam):
                    srows.append([b, L_prof, h, r + 1, v])
        write_csv(f"{out}/ed_profiles.csv", ["boundary", "L", "h", "l", "S"], rows)
        write_csv(f"{out}/entanglement_spectrum.csv", ["boundary", "L", "h", "rank", "lambda"], srows)

    # 2. ED energies
    with timed("2_ed_energies"):
        rows = []
        for b in BOUNDARIES:
            for L in range(8, L_en_max + 1, 2):
                for h in fields:
                    Ee, Eo = sector_energies(L, h, b)
                    rows.append([b, L, h, Ee, Eo, Eo - Ee])
        write_csv(f"{out}/ed_energies.csv", ["boundary", "L", "h", "E_even", "E_odd", "gap"], rows)

    # 3. ED half chain entropy at h = 1, periodic
    ed_half = {}
    with timed("3_ed_half_chain"):
        rows = []
        for L in range(8, L_half_max + 1, 2):
            _, psi = ground_state(L, 1.0, "periodic")
            s = entanglement_entropy(psi, L, L // 2)
            ed_half[L] = s
            rows.append(["periodic", L, 1.0, s])
        write_csv(f"{out}/ed_half_chain.csv", ["boundary", "L", "h", "S_half"], rows)

    # 4. Free fermion
    ff_half = {}
    ff_pts = {}
    with timed("4_free_fermion"):
        hrows, prows = [], []
        for L in ff_Ls:
            G, E0 = ff_solve(L, 1.0)
            s = block_entropy(G, L // 2)
            ff_half[L] = s
            hrows.append([L, 1.0, E0, s])
            if L == fit_L:
                ls = list(range(L // 16, L // 2 + 1, L // 64))
                vals = []
                for l in ls:
                    v = block_entropy(G, l)
                    vals.append(v)
                    prows.append([L, 1.0, l, v])
                ff_pts = {"ls": ls, "S": vals}
        for h in (1.0, 1.5, 2.0):
            G, E0 = ff_solve(L_ffprof, h)
            for l in range(L_ffprof + 1):
                prows.append([L_ffprof, h, l, block_entropy(G, l)])
        write_csv(f"{out}/ff_half_chain.csv", ["L", "h", "E0", "S_half"], hrows)
        write_csv(f"{out}/ff_profiles.csv", ["L", "h", "l", "S"], prows)

    # 5. Cross check
    with timed("5_cross_check"):
        rows = []
        for L in cross_Ls:
            for h in (1.0, 1.5, 2.0):
                E_ed, psi = ground_state(L, h, "open")
                S_ed = entropy_profile(psi, L)
                G, E_ff = ff_solve(L, h)
                dS = max(abs(S_ed[l] - block_entropy(G, l)) for l in range(L + 1))
                rows.append([L, h, dS, abs(E_ed - E_ff)])
        write_csv(f"{out}/ed_vs_ff.csv", ["L", "h", "max_abs_dS", "abs_dE"], rows)
        max_dS = max(r[2] for r in rows)
        max_dE = max(r[3] for r in rows)

    # 6. Field scan
    with timed("6_field_scan"):
        hs = np.round(np.arange(0.2, 2.0001, step), 10)
        rows = []
        for L in scan_Ls:
            for h in hs:
                _, psi = ground_state(L, float(h), "periodic")
                rows.append([L, float(h), entanglement_entropy(psi, L, L // 2)])
        write_csv(f"{out}/field_scan.csv", ["L", "h", "S_half"], rows)

    # 7. Central charge
    with timed("7_central_charge"):
        ff_pairs = [
            {"L1": a, "L2": b, "c": c_from_two_sizes(ff_half[a], a, ff_half[b], b, "open")}
            for a, b in zip(ff_Ls[:-1], ff_Ls[1:])
        ]
        eds = sorted(ed_half)
        ed_pairs = [
            {"L1": a, "L2": b, "c": c_from_two_sizes(ed_half[a], a, ed_half[b], b, "periodic")}
            for a, b in zip(eds[:-1], eds[1:])
        ]
        c_final = ff_pairs[-1]["c"]
        c_sys = abs(ff_pairs[-1]["c"] - ff_pairs[-2]["c"])
        c_fit_ff, _ = fit_c(ff_pts["ls"], ff_pts["S"], fit_L, "open")
        ls_ed = list(range(2, L_prof - 1))
        c_fit_ed, _ = fit_c(ls_ed, prof[("periodic", 1.0)][ls_ed], L_prof, "periodic")
        cc = {
            "ff_open_pairs": ff_pairs,
            "ed_periodic_pairs": ed_pairs,
            "c_final": c_final,
            "c_final_systematic": c_sys,
            "c_fit_ff_open_L1024": c_fit_ff,
            "c_fit_ed_periodic_L20": c_fit_ed,
            "c_ed_periodic_pair_18_20": ed_pairs[-1]["c"],
            "consistent": bool(abs(c_final - 0.5) < 3 * max(c_sys, 1e-4)),
        }
        with open(f"{out}/central_charge.json", "w") as f:
            json.dump(cc, f, indent=2)

    # 8. Summary
    summary = {
        "python": platform.python_version(),
        "packages": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "quick": q,
        **cc,
        "max_abs_dS": max_dS,
        "max_abs_dE": max_dE,
        "half_chain_entropies_profiles": half_profile,
        "run_time_seconds": times,
    }
    with open(f"{out}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("c_final = %.4f  systematic = %.4f  consistent = %s" % (c_final, c_sys, cc["consistent"]))
    print("max dS = %.2e  max dE = %.2e" % (max_dS, max_dE))


if __name__ == "__main__":
    main()
