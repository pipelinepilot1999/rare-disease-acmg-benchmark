#!/usr/bin/env python3
"""Split the single-threaded gnomAD cache into per-shard caches so the
parallel workers do not refetch what was already retrieved."""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import json, os, csv

N = 6
C = _paths.s(_paths.CACHE)
src_path = f"{C}/gnomad_faf.json"
src = json.load(open(src_path)) if os.path.exists(src_path) else {}
rows = list(csv.DictReader(open(_paths.s(_paths.REPO / "data/armb_variants.tsv")),
                           delimiter="\t"))
shards = [{} for _ in range(N)]
for i, r in enumerate(rows):
    vid = f'{r["chrom"][3:]}-{r["pos"]}-{r["ref"]}-{r["alt"]}'
    if vid in src:
        shards[i % N][vid] = src[vid]
for i, d in enumerate(shards):
    json.dump(d, open(f"{C}/gnomad_shard{i}.json", "w"))
print("seeded shards:", [len(d) for d in shards])
