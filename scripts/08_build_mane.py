#!/usr/bin/env python3
"""MANE Select transcript structures, for the PVS1 decision tree.

PVS1 turns on whether a premature termination codon triggers nonsense-mediated
decay: a PTC in the last exon, or within the last 50 nucleotides of the
penultimate exon, escapes NMD and the criterion is downgraded. That question
cannot be answered without exon coordinates, which is why v1's plan to take
PVS1 from a tool that does not model them was inadequate.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import gzip, json, os, re
from collections import defaultdict

RES = _paths.s(_paths.RES)
WORK = _paths.s(_paths.WORK)

tx = {}
exons = defaultdict(list)
cds = defaultdict(list)

attr_re = re.compile(r'(\S+) "([^"]*)"')

with gzip.open(f"{RES}/gencode/gencode.v50.basic.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] not in ("transcript", "exon", "CDS"):
            continue
        if "MANE_Select" not in f[8]:
            continue
        a = dict(attr_re.findall(f[8]))
        tid = a.get("transcript_id")
        if not tid:
            continue
        if f[2] == "transcript":
            tx[tid] = {"gene": a.get("gene_name"), "chrom": f[0], "strand": f[6],
                       "start": int(f[3]), "end": int(f[4]), "tx": tid}
        elif f[2] == "exon":
            exons[tid].append((int(f[3]), int(f[4])))
        else:
            cds[tid].append((int(f[3]), int(f[4])))

out = {}
for tid, t in tx.items():
    ex = sorted(exons.get(tid, []))
    cd = sorted(cds.get(tid, []))
    if not ex or not cd:
        continue
    if t["strand"] == "-":
        ex = ex[::-1]        # store in transcription order
    t["exons"] = ex
    t["cds_start"] = min(c[0] for c in cd)
    t["cds_end"] = max(c[1] for c in cd)
    t["n_exons"] = len(ex)
    out[t["gene"]] = t

os.makedirs(WORK, exist_ok=True)
json.dump(out, open(f"{WORK}/mane_select.json", "w"))
print(f"MANE Select transcripts: {len(out)} genes")
multi = sum(1 for t in out.values() if t["n_exons"] > 1)
print(f"  multi-exon: {multi}  single-exon: {len(out)-multi}")
ex = out.get("PAH")
if ex:
    print(f"  example PAH: {ex['tx']} {ex['chrom']}{ex['strand']} {ex['n_exons']} exons "
          f"CDS {ex['cds_start']}-{ex['cds_end']}")
