#!/usr/bin/env python3
"""Protein-change index for PS1 and PM5, built from ClinVar variant_summary.

PS1 asks whether the same amino acid change has already been called pathogenic
by a *different* nucleotide variant; PM5 asks whether a different change at the
same residue has. Both need an index keyed on (gene, residue), and both are
knowledge-derived, which is why they are removed in the blinded condition.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import gzip, os, re, csv, json
from collections import defaultdict

RES = _paths.s(_paths.RES)
REPO = _paths.s(_paths.REPO)
WORK = _paths.s(_paths.WORK)

P_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\)")

PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}

# residue -> {(aa_ref, aa_alt)} seen as P/LP, with the variation ids that support it
index = defaultdict(lambda: defaultdict(set))
selfmap = {}          # variation_id -> (gene, pos, aa_ref, aa_alt)
n = kept = 0

with gzip.open(f"{RES}/clinvar/variant_summary.txt.gz", "rt") as fh:
    head = fh.readline().lstrip("#").rstrip("\n").split("\t")
    ix = {k: head.index(k) for k in
          ("Name", "GeneSymbol", "ClinicalSignificance", "Assembly", "VariationID")}
    for line in fh:
        n += 1
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(ix.values()) or f[ix["Assembly"]] != "GRCh38":
            continue
        m = P_RE.search(f[ix["Name"]])
        if not m:
            continue
        gene = f[ix["GeneSymbol"]]
        aa_ref, pos, aa_alt = m.group(1), int(m.group(2)), m.group(3)
        vid = f[ix["VariationID"]]
        selfmap[vid] = (gene, pos, aa_ref, aa_alt)
        if f[ix["ClinicalSignificance"]] in PATHOGENIC:
            index[(gene, pos)][(aa_ref, aa_alt)].add(vid)
            kept += 1

print(f"variant_summary rows: {n}")
print(f"protein changes parsed: {len(selfmap)}")
print(f"P/LP protein changes indexed: {kept} across {len(index)} residues")

# ---------------------------------------------------- annotate the Arm B set
rows = list(csv.DictReader(open(f"{REPO}/data/armb_annotated.tsv"), delimiter="\t"))
hit_ps1 = hit_pm5 = has_aa = 0
for r in rows:
    vid = r["variation_id"]
    aa = selfmap.get(vid)
    r["aa_pos"] = r["aa_ref"] = r["aa_alt"] = ""
    r["ps1_hit"] = r["pm5_hit"] = "0"
    if not aa:
        continue
    gene, pos, aa_ref, aa_alt = aa
    has_aa += 1
    r["aa_pos"], r["aa_ref"], r["aa_alt"] = pos, aa_ref, aa_alt
    residue = index.get((gene, pos), {})
    # PS1: same amino acid substitution, established by a different variant
    others = residue.get((aa_ref, aa_alt), set()) - {vid}
    if others:
        r["ps1_hit"] = "1"; hit_ps1 += 1
    # PM5: a different substitution at the same residue, established as P/LP
    for (ar, al), vids in residue.items():
        if al != aa_alt and (vids - {vid}):
            r["pm5_hit"] = "1"; hit_pm5 += 1
            break

print(f"\nArm B variants with a protein change: {has_aa}/{len(rows)}")
print(f"  PS1 candidates: {hit_ps1}")
print(f"  PM5 candidates: {hit_pm5}")

cols = list(rows[0].keys())
with open(f"{REPO}/data/armb_annotated.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"updated {REPO}/data/armb_annotated.tsv")
