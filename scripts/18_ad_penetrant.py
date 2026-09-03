#!/usr/bin/env python3
"""Re-run the inherited-dominant cases with the transmitting parent coded affected.

The main run codes both parents unaffected, which for a dominant disease asserts
a non-penetrant carrier. Exomiser's inheritance filter correctly rejects a
heterozygous variant carried by an unaffected parent under AUTOSOMAL_DOMINANT,
so those cases were removed by the filter rather than mis-ranked. Coding the
father affected turns the flaw into a controlled contrast, and the difference
between the two is what the trio inheritance filter costs when penetrance is
incomplete.
"""
import json, os, sys, math

REPO = os.path.expanduser("~/rdacmg/repo")
JOBS = os.path.expanduser("~/rdacmg/jobs_pen")
OUT = os.path.expanduser("~/rdacmg/arma")
BDIR = os.path.expanduser("~/rdacmg/batches_pen")
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 6

sys.path.insert(0, f"{REPO}/scripts")
from make_job import job

for d in (JOBS, BDIR):
    os.makedirs(d, exist_ok=True)
for f in os.listdir(BDIR):
    os.remove(os.path.join(BDIR, f))

PED = f"{REPO}/config/trio_penetrant.ped"
manifest = [c for c in json.load(open(f"{REPO}/data/manifest.json"))
            if c["model"] == "ad_inherited"]
print(f"{len(manifest)} inherited-dominant cases")

lines = []
for c in manifest:
    name = f'{c["case_id"]}__ad_penetrant__as_reported__pos'
    y = job(c["vcf_pos"], c["hpo_sets"]["as_reported"], name, OUT, True)
    y = y.replace(f'ped: {REPO}/config/trio.ped', f'ped: {PED}')
    path = f"{JOBS}/{name}.yml"
    open(path, "w").write(y)
    lines.append(f"--analysis {path}")

per = math.ceil(len(lines) / NW)
n = 0
for w in range(NW):
    chunk = lines[w * per:(w + 1) * per]
    if not chunk:
        continue
    open(f"{BDIR}/informed_p{w}.batch", "w").write("\n".join(chunk) + "\n")
    n += 1
print(f"wrote {n} batch files to {BDIR}")
