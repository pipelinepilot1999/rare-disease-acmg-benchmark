#!/usr/bin/env python3
"""Merge the sharded gnomAD caches with the REVEL subset and write the
annotated Arm B table."""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, glob, json, os

REPO = _paths.s(_paths.REPO)
CACHE = _paths.s(_paths.CACHE)

faf = {}
for p in sorted(glob.glob(f"{CACHE}/gnomad_shard*.json")):
    d = json.load(open(p))
    faf.update(d)
    print(f"  {os.path.basename(p)}: {len(d)}")
legacy = f"{CACHE}/gnomad_faf.json"
if os.path.exists(legacy):
    for k, v in json.load(open(legacy)).items():
        faf.setdefault(k, v)
print(f"gnomAD entries: {len(faf)}")

revel = {}
for line in open(f"{CACHE}/revel_subset.tsv"):
    c, p, rf, al, s = line.rstrip("\n").split("\t")
    revel[(c, int(p), rf, al)] = float(s)
print(f"REVEL entries: {len(revel)}")

rows = list(csv.DictReader(open(f"{REPO}/data/armb_variants.tsv"), delimiter="\t"))
cols = list(rows[0].keys()) + ["gnomad_id", "revel", "faf95", "faf_group",
                               "gnomad_ac", "gnomad_an", "faf_source"]
missing_faf = 0
for r in rows:
    vid = f'{r["chrom"][3:]}-{r["pos"]}-{r["ref"]}-{r["alt"]}'
    r["gnomad_id"] = vid
    r["revel"] = revel.get((r["chrom"][3:], int(r["pos"]), r["ref"], r["alt"]), "")
    g = faf.get(vid)
    if g is None and vid in faf:
        # present in cache but absent from gnomAD: a truly absent variant has
        # filtering AF 0, which is what PM2 needs
        g = {"faf95": 0.0, "grp": None, "ac": 0, "an": None, "source": "absent"}
    if g is None:
        missing_faf += 1
        r.update(faf95="", faf_group="", gnomad_ac="", gnomad_an="", faf_source="not_queried")
    else:
        r.update(faf95=g["faf95"], faf_group=g["grp"], gnomad_ac=g["ac"],
                 gnomad_an=g["an"], faf_source=g["source"])

with open(f"{REPO}/data/armb_annotated.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader(); w.writerows(rows)

from collections import Counter
print(f"\nrows: {len(rows)}  never queried: {missing_faf}")
print("faf source:", Counter(r["faf_source"] for r in rows).most_common())
print(f"with REVEL: {sum(1 for r in rows if r['revel'] != '')}")
print(f"wrote {REPO}/data/armb_annotated.tsv")
