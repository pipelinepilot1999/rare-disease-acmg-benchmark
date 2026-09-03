#!/usr/bin/env python3
"""Split the single-threaded gnomAD cache into per-shard caches so the
parallel workers do not refetch what was already retrieved."""
import json, os, csv

N = 6
C = os.path.expanduser("~/rdacmg/work/cache")
src_path = f"{C}/gnomad_faf.json"
src = json.load(open(src_path)) if os.path.exists(src_path) else {}
rows = list(csv.DictReader(open(os.path.expanduser("~/rdacmg/repo/data/armb_variants.tsv")),
                           delimiter="\t"))
shards = [{} for _ in range(N)]
for i, r in enumerate(rows):
    vid = f'{r["chrom"][3:]}-{r["pos"]}-{r["ref"]}-{r["alt"]}'
    if vid in src:
        shards[i % N][vid] = src[vid]
for i, d in enumerate(shards):
    json.dump(d, open(f"{C}/gnomad_shard{i}.json", "w"))
print("seeded shards:", [len(d) for d in shards])
