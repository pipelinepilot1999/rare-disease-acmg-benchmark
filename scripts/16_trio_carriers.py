#!/usr/bin/env python3
"""Enumerate the real ClinVar P/LP alleles the trio already carries.

A knockout case is negative only for the spiked variant. HG002 is a healthy
individual but not a genome free of pathogenic alleles, and saying so with a
list is better than saying so in the abstract.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, gzip, os, subprocess
from collections import Counter

RES = _paths.s(_paths.RES)
WORK = _paths.s(_paths.WORK)
REPO = _paths.s(_paths.REPO)

trio = {}
p = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
                    f"{WORK}/trio.vcf.gz"], capture_output=True, text=True, check=True)
for line in p.stdout.splitlines():
    f = line.split("\t")
    trio[(f[0][3:], int(f[1]), f[2], f[3])] = f[4:7]
print(f"trio variants: {len(trio)}")

rows = []
with gzip.open(f"{RES}/clinvar/clinvar.vcf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        key = (f[0], int(f[1]), f[3], f[4])
        gts = trio.get(key)
        if not gts:
            continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        sig = (info.get("CLNSIG") or "").replace("_", " ")
        low = sig.lower()
        if "pathogenic" not in low or "conflicting" in low:
            continue
        if any("1" in g for g in gts):
            rows.append({
                "chrom": "chr" + f[0], "pos": f[1], "ref": f[3], "alt": f[4],
                "variation_id": f[2], "gene": (info.get("GENEINFO") or "").split(":")[0],
                "clnsig": sig, "revstat": (info.get("CLNREVSTAT") or "").replace("_", " "),
                "condition": (info.get("CLNDN") or "").replace("_", " ")[:80],
                "HG002": gts[0], "HG003": gts[1], "HG004": gts[2],
            })

print(f"\nClinVar P/LP alleles present in the trio: {len(rows)}")
print("  by proband genotype:", Counter(r["HG002"] for r in rows).most_common())
hom = [r for r in rows if r["HG002"] in ("1/1", "1|1")]
print(f"  homozygous in the proband: {len(hom)}")
for r in hom[:8]:
    print(f"    {r['gene']:<12s} {r['chrom']}:{r['pos']} {r['clnsig']:<22s} {r['condition'][:46]}")

with open(f"{REPO}/data/trio_carriers.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {REPO}/data/trio_carriers.tsv")
