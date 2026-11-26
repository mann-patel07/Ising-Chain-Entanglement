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
for key in ["max_abs_dS", "max_abs_dE"]:
    t = f"{s[key]:.1e}"
    check(t in readme, f"README is missing {key} = {t}")

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
files = ["README.md", "PREREGISTRATION.md", "LICENSE", ".gitignore", "requirements.txt"]
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
