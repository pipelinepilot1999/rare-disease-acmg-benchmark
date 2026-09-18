#!/usr/bin/env python3
"""Fetch gnomAD v4 joint grpmax filtering AF for one shard of the Arm B set.

Sharded because the work is latency-bound: the API caps a query at 25 variant
aliases, so throughput comes from concurrent requests, not bigger batches.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, json, os, sys, time, urllib.request

shard, nshard = int(sys.argv[1]), int(sys.argv[2])
REPO = _paths.s(_paths.REPO)
CACHE = _paths.s(_paths.CACHE)
os.makedirs(CACHE, exist_ok=True)
OUT = f"{CACHE}/gnomad_shard{shard}.json"

rows = list(csv.DictReader(open(f"{REPO}/data/armb_variants.tsv"), delimiter="\t"))
ids = [f'{r["chrom"][3:]}-{r["pos"]}-{r["ref"]}-{r["alt"]}'
       for i, r in enumerate(rows) if i % nshard == shard]
ids = list(dict.fromkeys(ids))

faf = json.load(open(OUT)) if os.path.exists(OUT) else {}
todo = [v for v in ids if v not in faf]
print(f"shard {shard}/{nshard}: {len(ids)} ids, {len(todo)} to fetch", flush=True)

API = "https://gnomad.broadinstitute.org/api"
F = "{ac an faf95{popmax popmax_population}}"
BATCH = 20

def fetch(chunk):
    parts = [f'v{i}: variant(variantId: "{v}", dataset: gnomad_r4) '
             f'{{ variant_id joint{F} exome{F} genome{F} }}' for i, v in enumerate(chunk)]
    body = json.dumps({"query": "query {" + " ".join(parts) + "}"}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "rare-disease-acmg-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode())

for start in range(0, len(todo), BATCH):
    chunk = todo[start:start + BATCH]
    d = None
    for attempt in range(7):
        try:
            d = fetch(chunk); break
        except Exception as e:
            wait = min(180, 6 * (2 ** attempt))
            print(f"  shard{shard} retry {wait}s ({e})", flush=True)
            time.sleep(wait)
    if d is None:
        print(f"  shard{shard} GIVING UP at {start}", flush=True); continue
    data = d.get("data") or {}
    for i, vid in enumerate(chunk):
        v = data.get(f"v{i}")
        if not v:
            faf[vid] = None; continue
        best = None
        for src in ("joint", "exome", "genome"):
            blk = v.get(src)
            if blk and blk.get("faf95") and blk["faf95"].get("popmax") is not None:
                best = {"faf95": blk["faf95"]["popmax"],
                        "grp": blk["faf95"].get("popmax_population"),
                        "ac": blk.get("ac"), "an": blk.get("an"), "source": src}
                break
        if best is None:
            blk = v.get("joint") or v.get("exome") or v.get("genome") or {}
            best = {"faf95": 0.0, "grp": None, "ac": blk.get("ac"),
                    "an": blk.get("an"), "source": "present_no_faf"}
        faf[vid] = best
    if (start // BATCH) % 20 == 0:
        json.dump(faf, open(OUT, "w"))
        print(f"  shard{shard}: {start+len(chunk)}/{len(todo)}", flush=True)
    time.sleep(0.3)

json.dump(faf, open(OUT, "w"))
print(f"shard {shard} DONE: {len(faf)} entries", flush=True)
