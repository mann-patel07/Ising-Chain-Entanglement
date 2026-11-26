"""Consistency checks between results/, README.md and the files to upload."""
import json
import os
import re
import sys

fails = []
n = 0


def check(ok, msg):
    global n
    n += 1
    if not ok:
        fails.append(msg)


with open("results/summary.json") as f:
    s = json.load(f)
with open("README.md") as f:
    readme = f.read()

for key in ["c_final", "c_final_systematic", "c_fit_ff_open_L1024", "c_fit_ed_periodic_L20", "c_ed_periodic_pair_18_20"]:
    t = f"{s[key]:.4f}"
    check(t in readme, f"README is missing {key} = {t}")
# Rounding-level differences are stated as a bound, because their exact values depend on the computer.
for key in ["max_abs_dS", "max_abs_dE"]:
    check(s[key] < 1e-11, f"{key} = {s[key]:.1e} is not below the 1e-11 stated in the README")

# periodic free fermion section (added after the first version)
with open("results/central_charge_periodic.json") as f:
    cp = json.load(f)
for q in cp["ff_periodic_pairs"]:
    row = f"| Two sizes, free fermion, periodic, h = 1 | L = {q['L1']} and {q['L2']} | {q['c']:.6f} |"
    check(row in readme, f"README is missing the periodic row {row}")
check(f"| Fit of S(l), free fermion, periodic, h = 1 | L = 1024, l = 64 to 512 in steps of 16 | {cp['c_fit_ff_periodic_L1024']:.6f} |" in readme, "README is missing the periodic fit row")
check(f"give c = {cp['c_periodic_final']:.6f}." in readme, "README first paragraph is missing the periodic c")
check(f"The last doubling gives c = {cp['c_periodic_final']:.6f}." in readme, "README is missing the last periodic doubling")
for key in ["max_abs_dS", "max_abs_dE"]:
    check(cp[key] < 1e-11, f"periodic {key} = {cp[key]:.1e} is not below the 1e-11 stated in the README")
check(cp["c_periodic_systematic"] < 1e-6, "periodic systematic is not below the 1e-6 stated in the README")
check(abs(cp["c_periodic_final"] - 0.5) < 1e-6, "periodic distance from 1/2 is not below the 1e-6 stated in the README")
check(f"the open chain is {abs(s['c_final'] - 0.5):.1e} away from 1/2 and the periodic chain is less than 1e-6 away" in readme, "README distances for L = 1024 and 2048 are wrong")
check(f"the open chain is {abs(s['ff_open_pairs'][0]['c'] - 0.5):.1e} away and the periodic chain is {abs(cp['ff_periodic_pairs'][0]['c'] - 0.5):.1e} away" in readme, "README distances for L = 16 and 32 are wrong")
check(all(q["c"] < 0.5 for q in cp["ff_periodic_pairs"]), "periodic pairs should all be below 0.5")
check(os.path.exists("UPGRADE_PLAN.md"), "UPGRADE_PLAN.md is missing")

if s["consistent"]:
    check("consistent" in readme, "README should say consistent")
else:
    check("not consistent" in readme, "README should say not consistent")

for p in re.findall(r"\]\((figures/[^)]+)\)", readme):
    check(os.path.exists(p), f"missing figure {p}")

EXTS = (".py", ".md", ".txt", ".csv", ".json")
bad = [chr(c) for c in (0x2014, 0x2013, 0x201C, 0x201D, 0x2018, 0x2019)]
bad += ["/Users/"]
# Only the files that are uploaded are scanned.
files = ["README.md", "PREREGISTRATION.md", "UPGRADE_PLAN.md", "LICENSE", ".gitignore", "requirements.txt"]
files += sorted(f for f in os.listdir(".") if f.endswith(".py"))
for d in ("tests", "results", "figures"):
    files += sorted(os.path.join(d, f) for f in os.listdir(d) if not f.startswith(".") and f != "__pycache__")
for p in files:
    if p == "check_consistency.py":
        continue  # this file names the patterns it looks for
    if p.endswith(EXTS) or os.path.basename(p) in ("LICENSE", ".gitignore"):
        with open(p, encoding="utf-8") as f:
            txt = f.read()
        for b in bad:
            check(b not in txt, f"{p} contains {b!r}")

table = [l for l in readme.splitlines() if l.startswith("|")]
if table:
    k = table[0].count("|")
    for l in table:
        check(l.count("|") == k, f"table row has {l.count('|')} pipes, header has {k}: {l}")

check(len(files) < 100, f"too many files: {len(files)}")
for p in files:
    check(os.path.getsize(p) < 25 * 1024 * 1024, f"{p} is over 25 MiB")

if fails:
    for m in fails:
        print("FAIL:", m)
    sys.exit(1)
print("all checks passed (%d checks)" % n)
