#!/usr/bin/env python3
"""Spike causal alleles into the trio, one VCF per case, plus matched knockouts.

FORMAT fields for spiked records are sampled from real trio variants matched on
variant type and genotype, so a spiked variant is not separable from a real one
by its quality fields. Without this the quality filter never touches a causal
variant and Arm A recall is inflated.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import json, os, random, subprocess, sys, gzip
from collections import defaultdict

WORK = _paths.s(_paths.WORK)
REPO = _paths.s(_paths.REPO)
OUT = _paths.s(_paths.CASES)
random.seed(20260902)
os.makedirs(OUT, exist_ok=True)

CHROM_ORDER = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}

# ------------------------------------------------ read background once
header, body = [], []
with gzip.open(f"{WORK}/trio.vcf.gz", "rt") as fh:
    for line in fh:
        (header if line.startswith("#") else body).append(line.rstrip("\n"))
print(f"background: {len(body)} variants, {len(header)} header lines")

FORMAT_KEYS = "GT:DP:GQ"

# ------------------------------------------------ FORMAT donor pool
# key: (is_indel, genotype) -> list of (DP, GQ) seen in real data
donors = defaultdict(list)
for line in body:
    f = line.split("\t")
    is_indel = len(f[3]) != 1 or len(f[4]) != 1
    fmt = f[8].split(":")
    for s in f[9:]:
        parts = dict(zip(fmt, s.split(":")))
        gt = (parts.get("GT") or "./.").replace("|", "/")
        dp, gq = parts.get("DP", "."), parts.get("GQ", ".")
        if dp not in (".", "") and gq not in (".", ""):
            donors[(is_indel, gt)].append((dp, gq))
for k in list(donors):
    random.shuffle(donors[k])
print("FORMAT donor pool:", {f"{'indel' if k[0] else 'snv'}/{k[1]}": len(v)
                             for k, v in sorted(donors.items())})

def sample_format(is_indel, gt):
    pool = donors.get((is_indel, gt)) or donors.get((False, gt)) or donors.get((is_indel, "0/1"))
    dp, gq = random.choice(pool)
    return f"{gt}:{dp}:{gq}"

# ------------------------------------------------ genotype assignment
def genotypes(model, index):
    """(proband, father, mother) for variant `index` of the case."""
    return {
        "denovo":       [("0/1", "0/0", "0/0")],
        "ad_inherited": [("0/1", "0/1", "0/0")],
        "ar_hom":       [("1/1", "0/1", "0/1")],
        "xlr":          [("1/1", "0/0", "0/1")],   # hemizygous male proband
        "comphet":      [("0/1", "0/1", "0/0"), ("0/1", "0/0", "0/1")],
    }[model][index]

cases = json.load(open(f"{REPO}/data/cases.json"))
print(f"{len(cases)} cases")

# a background record already at a spiked position must be replaced, not duplicated
bg_index = {}
for i, line in enumerate(body):
    f = line.split("\t")
    bg_index[(f[0], int(f[1]), f[3], f[4])] = i

manifest = []
for n, case in enumerate(cases, 1):
    spiked_lines, truth = [], []
    for idx, v in enumerate(case["variants"]):
        gt_p, gt_f, gt_m = genotypes(case["model"], idx if case["model"] == "comphet" else 0)
        is_indel = len(v["ref"]) != 1 or len(v["alt"]) != 1
        cols = [v["chrom"], str(v["pos"]), ".", v["ref"], v["alt"], "50", "PASS", "SPIKED=1",
                FORMAT_KEYS,
                sample_format(is_indel, gt_p),
                sample_format(is_indel, gt_f),
                sample_format(is_indel, gt_m)]
        spiked_lines.append("\t".join(cols))
        truth.append({"chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
                      "gene": v["gene"], "gt": [gt_p, gt_f, gt_m]})

    for tag, extra in (("pos", spiked_lines), ("ko", [])):
        recs = list(body) + extra
        recs.sort(key=lambda l: (CHROM_ORDER.get(l.split("\t", 1)[0], 99),
                                 int(l.split("\t")[1])))
        path = f"{OUT}/{case['case_id']}.{tag}.vcf.gz"
        with gzip.open(path, "wt", compresslevel=1) as fh:
            for h in header:
                if h.startswith("#CHROM"):
                    fh.write('##INFO=<ID=SPIKED,Number=0,Type=Flag,Description="Injected causal allele">\n')
                fh.write(h + "\n")
            fh.write("\n".join(recs) + "\n")

    manifest.append({**{k: case[k] for k in
                        ("case_id", "gene", "disease_id", "disease", "model", "hpo", "hpo_sets")},
                     "truth": truth,
                     "vcf_pos": f"{OUT}/{case['case_id']}.pos.vcf.gz",
                     "vcf_ko": f"{OUT}/{case['case_id']}.ko.vcf.gz"})
    if n % 200 == 0:
        print(f"  {n}/{len(cases)}", flush=True)

with open(f"{REPO}/data/manifest.json", "w") as fh:
    json.dump(manifest, fh, indent=1)
print(f"wrote {len(manifest)} case pairs -> {OUT}")
print(f"manifest -> {REPO}/data/manifest.json")
