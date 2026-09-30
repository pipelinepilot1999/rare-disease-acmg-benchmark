#!/usr/bin/env python3
"""Generate Exomiser job files and batch manifests for Arm A.

Run matrix
----------
positives  x {informed, clinvar_blind, pheno_blind, double_blind}
knockouts  x {informed, clinvar_blind}          (matched negative controls)
subset     x {imprecise, noisy, truncated}      (HPO perturbation, informed only)

The phenotype-blinded conditions ignore HPO entirely, so they are not repeated
across perturbation tiers.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import json, os, random, sys, math

REPO = _paths.s(_paths.REPO)
JOBS = _paths.s(_paths.JOBS)
OUT = _paths.s(_paths.ARMA)
NWORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 8
HPO_SUBSET = int(sys.argv[2]) if len(sys.argv) > 2 else 300
random.seed(20260902)

sys.path.insert(0, f"{REPO}/scripts")
from make_job import job

for d in (JOBS, OUT):
    os.makedirs(d, exist_ok=True)


def vcf(case, kind):
    """Per-case VCF path, derived rather than stored.

    The manifest used to carry absolute paths, which pinned it to one machine
    and put 2,068 of them in a public file. They are a function of the case id
    and RDACMG_HOME, so they are computed here instead.
    """
    return _paths.s(_paths.CASES / f"{case['case_id']}.{kind}.vcf.gz")

manifest = json.load(open(f"{REPO}/data/manifest.json"))
print(f"{len(manifest)} cases")

subset_ids = {c["case_id"] for c in random.sample(manifest, min(HPO_SUBSET, len(manifest)))}

specs = []   # (case_id, condition, hpo_tier, vcf, hpo_terms, phenotype_on)
for c in manifest:
    hpo = c["hpo_sets"]["as_reported"]
    for cond, pheno in (("informed", True), ("clinvar_blind", True),
                        ("pheno_blind", False), ("double_blind", False)):
        specs.append((c["case_id"], cond, "as_reported", vcf(c, "pos"), hpo, pheno))
    for cond in ("informed", "clinvar_blind"):
        specs.append((c["case_id"], cond, "as_reported", vcf(c, "ko"), hpo, True))
    if c["case_id"] in subset_ids:
        for tier in ("imprecise", "noisy", "truncated"):
            specs.append((c["case_id"], "informed", tier, vcf(c, "pos"),
                          c["hpo_sets"][tier], True))

print(f"{len(specs)} Exomiser runs")

# properties file per condition: the two blind-to-ClinVar conditions share one
PROPS = {"informed": "informed", "clinvar_blind": "clinvar_blind",
         "pheno_blind": "informed", "double_blind": "clinvar_blind"}

lines_by_prop = {"informed": [], "clinvar_blind": []}
for case_id, cond, tier, vcf, hpo, pheno in specs:
    kind = "ko" if vcf.endswith(".ko.vcf.gz") else "pos"
    name = f"{case_id}__{cond}__{tier}__{kind}"
    path = f"{JOBS}/{name}.yml"
    with open(path, "w") as fh:
        fh.write(job(vcf, hpo, name, OUT, pheno))
    lines_by_prop[PROPS[cond]].append(f"--analysis {path}")

os.makedirs(f"{REPO}/../batches", exist_ok=True)
BDIR = _paths.s(_paths.BATCHES)
os.makedirs(BDIR, exist_ok=True)
for f in os.listdir(BDIR):
    os.remove(os.path.join(BDIR, f))

total_files = 0
for prop, lines in lines_by_prop.items():
    random.shuffle(lines)
    per = math.ceil(len(lines) / NWORKERS)
    for w in range(NWORKERS):
        chunk = lines[w * per:(w + 1) * per]
        if not chunk:
            continue
        p = f"{BDIR}/{prop}_w{w}.batch"
        open(p, "w").write("\n".join(chunk) + "\n")
        total_files += 1
    print(f"  {prop:14s} {len(lines):6d} jobs -> {NWORKERS} batch files")

print(f"wrote {total_files} batch files to {BDIR}")
print(f"job YAMLs: {len(specs)} in {JOBS}")
