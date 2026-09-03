#!/usr/bin/env python3
"""Annotate the Arm B variant set with the evidence the criteria engine needs.

* gnomAD v4 joint grpmax filtering AF (faf95.popmax) - the value ClinGen's
  March 2024 guidance names for BA1/BS1. Not raw popmax: the filtering AF is the
  95% CI lower bound, and it excludes bottlenecked genetic ancestry groups.
* REVEL v1.3 - the single calibrated missense predictor used for PP3/BP4.

Both are cached to disk so the step is restartable.
"""
import csv, json, os, sys, time, zipfile, io, gzip
import urllib.request
from collections import defaultdict

REPO = os.path.expanduser("~/rdacmg/repo")
RES = os.path.expanduser("~/rdacmg/resources")
CACHE = os.path.expanduser("~/rdacmg/work/cache")
os.makedirs(CACHE, exist_ok=True)

rows = list(csv.DictReader(open(f"{REPO}/data/armb_variants.tsv"), delimiter="\t"))
print(f"{len(rows)} Arm B variants")

# ------------------------------------------------------------- REVEL
REVEL_CACHE = f"{CACHE}/revel_subset.tsv"
wanted = {(r["chrom"][3:], int(r["pos"]), r["ref"], r["alt"]) for r in rows}
revel = {}
if os.path.exists(REVEL_CACHE):
    for line in open(REVEL_CACHE):
        c, p, rf, al, s = line.rstrip("\n").split("\t")
        revel[(c, int(p), rf, al)] = float(s)
    print(f"REVEL cache: {len(revel)} scores")
else:
    print("scanning REVEL v1.3 ...")
    zf = zipfile.ZipFile(f"{RES}/revel/revel-v1.3.zip")
    name = zf.namelist()[0]   # single member: revel_with_transcript_ids (CSV)
    n = 0
    with zf.open(name) as raw, open(REVEL_CACHE, "w") as out:
        head = raw.readline().decode().rstrip("\n").split(",")
        ix = {k: head.index(k) for k in ("chr", "grch38_pos", "ref", "alt", "REVEL")}
        for bline in raw:
            f = bline.decode().rstrip("\n").split(",")
            pos = f[ix["grch38_pos"]]
            if pos in (".", ""):
                continue
            key = (f[ix["chr"]], int(pos), f[ix["ref"]], f[ix["alt"]])
            if key in wanted:
                score = float(f[ix["REVEL"]])
                # a position can carry several transcript rows; keep the max
                if revel.get(key, -1) < score:
                    revel[key] = score
            n += 1
            if n % 20_000_000 == 0:
                print(f"  {n/1e6:.0f}M rows, {len(revel)} hits", flush=True)
        for (c, p, rf, al), s in sorted(revel.items()):
            out.write(f"{c}\t{p}\t{rf}\t{al}\t{s}\n")
    print(f"REVEL: {len(revel)} of {len(wanted)} variants scored")

# ------------------------------------------------------------- gnomAD
GNOMAD_CACHE = f"{CACHE}/gnomad_faf.json"
faf = json.load(open(GNOMAD_CACHE)) if os.path.exists(GNOMAD_CACHE) else {}
print(f"gnomAD cache: {len(faf)} variants")

API = "https://gnomad.broadinstitute.org/api"
FIELDS = "{ac an faf95{popmax popmax_population}}"

def query_batch(ids):
    parts = []
    for i, vid in enumerate(ids):
        parts.append(f'v{i}: variant(variantId: "{vid}", dataset: gnomad_r4) '
                     f'{{ variant_id joint{FIELDS} exome{FIELDS} genome{FIELDS} }}')
    body = json.dumps({"query": "query {" + " ".join(parts) + "}"}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "rare-disease-acmg-benchmark/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode())

todo = []
for r in rows:
    vid = f'{r["chrom"][3:]}-{r["pos"]}-{r["ref"]}-{r["alt"]}'
    r["gnomad_id"] = vid
    if vid not in faf:
        todo.append(vid)
todo = list(dict.fromkeys(todo))
print(f"gnomAD lookups needed: {len(todo)}")

BATCH = 20   # gnomAD caps query cost at 25 aliases per request
for start in range(0, len(todo), BATCH):
    chunk = todo[start:start + BATCH]
    for attempt in range(6):
        try:
            d = query_batch(chunk)
            break
        except Exception as e:
            wait = min(120, 8 * (2 ** attempt))
            print(f"  retry in {wait}s ({e})", flush=True)
            time.sleep(wait)
    else:
        print(f"  GIVING UP on batch at {start}"); continue

    data = d.get("data") or {}
    for i, vid in enumerate(chunk):
        v = data.get(f"v{i}")
        if not v:
            faf[vid] = None            # absent from gnomAD entirely
            continue
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
    if (start // BATCH) % 10 == 0:
        json.dump(faf, open(GNOMAD_CACHE, "w"))
        print(f"  {start+len(chunk)}/{len(todo)}", flush=True)
    time.sleep(1.0)

json.dump(faf, open(GNOMAD_CACHE, "w"))
print(f"gnomAD: {sum(1 for v in faf.values() if v)} with data, "
      f"{sum(1 for v in faf.values() if not v)} absent")

# ------------------------------------------------------------- write
out_cols = list(rows[0].keys()) + ["revel", "faf95", "faf_group", "gnomad_ac", "gnomad_an", "faf_source"]
with open(f"{REPO}/data/armb_annotated.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        key = (r["chrom"][3:], int(r["pos"]), r["ref"], r["alt"])
        r["revel"] = revel.get(key, "")
        g = faf.get(r["gnomad_id"])
        r["faf95"] = g["faf95"] if g else ""
        r["faf_group"] = g["grp"] if g else ""
        r["gnomad_ac"] = g["ac"] if g else ""
        r["gnomad_an"] = g["an"] if g else ""
        r["faf_source"] = g["source"] if g else "absent"
        w.writerow(r)
print(f"wrote {REPO}/data/armb_annotated.tsv")
