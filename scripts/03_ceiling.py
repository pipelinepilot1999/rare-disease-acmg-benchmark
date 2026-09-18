#!/usr/bin/env python3
"""Arm B' - the ceiling on automated ACMG classification.

Takes each ClinGen VCEP's own evidence codes, removes the ones no automated
system could assess, re-runs the combining rules, and records the tier that
survives. Aggregated, this bounds what any rule engine can achieve - a property
of the framework and the evidence distribution, not of any particular tool.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import json, os, sys, csv
from collections import Counter, defaultdict

sys.path.insert(0, _paths.s(_paths.REPO / "src"))
from acmgx.combine import combine_2015, combine_points, parse_code
from acmgx.assessability import bucket

RAW = _paths.s(_paths.RES / "erepo/erepo_raw.json")
OUTDIR = _paths.s(_paths.REPO / "results")
os.makedirs(OUTDIR, exist_ok=True)

TIER_MAP = {
    "Pathogenic": "P",
    "Likely Pathogenic": "LP",
    "Uncertain Significance": "VUS",
    "Uncertain Significance - Insufficient Evidence": "VUS",
    "Uncertain Significance - Conflicting Evidence": "VUS",
    "Likely Benign": "LB",
    "Benign": "B",
    "Benign - Stand Alone": "B",
}

records = json.load(open(RAW))
print(f"loaded {len(records)} ERepo interpretations")

rows, skipped = [], Counter()
for r in records:
    caid = r.get("caid")
    cond = (r.get("condition") or {}).get("label")
    gene = (r.get("gene") or {}).get("label") or (r.get("gene") or {}).get("@id")
    for g in r.get("guidelines") or []:
        # The final call sits on the guideline. Agents carry the evidence codes,
        # and only sometimes repeat the outcome.
        label = (g.get("outcome") or {}).get("label")
        tier = TIER_MAP.get(label)
        if tier is None:
            skipped[f"outcome:{label}"] += 1
            continue
        met, vceps = [], []
        for a in g.get("agents") or []:
            name = a.get("affiliation") or a.get("label")
            if name:
                vceps.append(name)
            for c in a.get("evidenceCodes") or []:
                if c.get("status") == "Met" and c.get("label"):
                    met.append(c["label"])
        good = []
        for code in dict.fromkeys(met):          # de-duplicate, preserve order
            try:
                parse_code(code); good.append(code)
            except ValueError:
                skipped[f"code:{code}"] += 1
        if not good:
            skipped["no_met_codes"] += 1
            continue
        rows.append({
            "caid": caid, "condition": cond, "gene": gene,
            "vcep": vceps[0] if vceps else None,
            "tier_vcep": tier, "codes": good,
        })

print(f"usable interpretations: {len(rows)}")
print("skipped outcomes/codes:", dict(skipped.most_common(12)))

# ---------------------------------------------------------------- analysis
def subset(codes, keep):
    out = []
    for c in codes:
        base, _, _ = parse_code(c)
        if bucket(base) in keep:
            out.append(c)
    return out

for row in rows:
    codes = row["codes"]
    row["n_codes"] = len(codes)
    row["tier_recomputed"], row["pts_recomputed"] = (combine_2015(codes), combine_points(codes)[1])
    ass = subset(codes, {"assessable", "knowledge"})
    bli = subset(codes, {"assessable"})
    row["tier_ceiling"], row["pts_ceiling"] = combine_2015(ass), combine_points(ass)[1]
    row["tier_ceiling_blind"], row["pts_ceiling_blind"] = combine_2015(bli), combine_points(bli)[1]
    row["n_unassessable"] = sum(1 for c in codes if bucket(parse_code(c)[0]) == "unassessable")
    row["n_knowledge"] = sum(1 for c in codes if bucket(parse_code(c)[0]) == "knowledge")

def agree(key, ref="tier_vcep"):
    n = sum(1 for r in rows if r[key] == r[ref])
    return n, n / len(rows)

print()
print("=" * 74)
print("VALIDATION - do our combining rules reproduce the VCEPs' own calls?")
n, f = agree("tier_recomputed")
print(f"  recomputed from ALL VCEP codes == VCEP tier : {n}/{len(rows)} ({f:.1%})")
print()
print("CEILING - re-running the VCEPs' evidence with unassessable criteria removed")
n, f = agree("tier_ceiling")
print(f"  machine-assessable + knowledge-derived      : {n}/{len(rows)} ({f:.1%})")
n, f = agree("tier_ceiling_blind")
print(f"  machine-assessable only (ClinVar-blinded)   : {n}/{len(rows)} ({f:.1%})")
print("=" * 74)

# dependence on unassessable evidence
dep = Counter()
for r in rows:
    dep[(r["tier_vcep"], r["tier_ceiling"])] += 1
print("\nTier flips when unassessable evidence is removed (VCEP -> ceiling), top 15:")
for (a, b), n in dep.most_common(15):
    mark = "" if a == b else "   <-- lost"
    print(f"  {a:4s} -> {b:4s}  {n:6d}{mark}")

print("\nFraction of interpretations depending on >=1 unassessable criterion: "
      f"{sum(1 for r in rows if r['n_unassessable']>0)/len(rows):.1%}")

by_tier = defaultdict(lambda: [0, 0])
for r in rows:
    by_tier[r["tier_vcep"]][0] += 1
    by_tier[r["tier_vcep"]][1] += int(r["tier_vcep"] == r["tier_ceiling"])
print("\nCeiling recall by expert tier:")
for t in ("P", "LP", "VUS", "LB", "B"):
    tot, ok = by_tier[t]
    if tot:
        print(f"  {t:4s} n={tot:6d}  retained {ok:6d} ({ok/tot:.1%})")

# criterion frequency and how often each is the load-bearing one
crit = Counter(); crit_bucket = Counter()
for r in rows:
    for c in r["codes"]:
        base = parse_code(c)[0]
        crit[base] += 1
        crit_bucket[bucket(base)] += 1
print("\nEvidence-code usage across all VCEP interpretations:")
tot_codes = sum(crit.values())
for b in ("assessable", "knowledge", "unassessable", "unknown"):
    print(f"  {b:14s} {crit_bucket[b]:7d} ({crit_bucket[b]/tot_codes:.1%})")
print("\n  most-used criteria:", ", ".join(f"{k}={v}" for k, v in crit.most_common(18)))

with open(f"{OUTDIR}/ceiling.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["caid","vcep","gene","condition","tier_vcep","tier_recomputed","tier_ceiling",
                "tier_ceiling_blind","pts_recomputed","pts_ceiling","pts_ceiling_blind",
                "n_codes","n_unassessable","n_knowledge","codes"])
    for r in rows:
        w.writerow([r["caid"],r["vcep"],r["gene"],r["condition"],r["tier_vcep"],r["tier_recomputed"],
                    r["tier_ceiling"],r["tier_ceiling_blind"],r["pts_recomputed"],
                    r["pts_ceiling"],r["pts_ceiling_blind"],r["n_codes"],
                    r["n_unassessable"],r["n_knowledge"],";".join(r["codes"])])
print(f"\nwrote {OUTDIR}/ceiling.tsv")
