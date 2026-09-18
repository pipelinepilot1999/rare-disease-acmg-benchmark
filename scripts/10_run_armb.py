#!/usr/bin/env python3
"""Arm B: run acmgx over the ERepo ground truth under three blinding conditions,
alongside a ClinVar lookup baseline that performs no reasoning at all.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, json, os, sys
from collections import Counter, defaultdict

sys.path.insert(0, _paths.s(_paths.REPO / "src"))
from acmgx.criteria import Variant, evaluate, CRITERIA, NOT_IMPLEMENTED
from acmgx.combine import combine_2015, combine_points, parse_code

REPO = _paths.s(_paths.REPO)
WORK = _paths.s(_paths.WORK)
RES = _paths.s(_paths.RES)

CONDITIONS = {
    "full":  set(),                                  # includes deprecated PP5/BP6
    "svi":   {"PP5", "BP6"},                         # current ClinGen guidance
    "blind": {"PP5", "BP6", "PS1", "PM5"},           # no ClinVar-derived evidence
}
TIERS = ["P", "LP", "VUS", "LB", "B"]

mane = json.load(open(f"{WORK}/mane_select.json"))

hi = {}
for line in open(f"{RES}/clingen/dosage.tsv"):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) > 4:
        hi[f[0]] = f[4].strip()
print(f"ClinGen dosage curations: {len(hi)} genes")

rows = list(csv.DictReader(open(f"{REPO}/data/armb_annotated.tsv"), delimiter="\t"))
print(f"Arm B variants: {len(rows)}")

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

results = []
for r in rows:
    v = Variant(
        chrom=r["chrom"], pos=int(r["pos"]), ref=r["ref"], alt=r["alt"],
        gene=r["gene"], consequence=r["consequence"],
        revel=num(r.get("revel")), faf95=num(r.get("faf95")),
        clinvar_tier=r.get("clinvar_tier") or None,
        clinvar_revstat=r.get("clnrevstat", ""),
        aa_pos=int(r["aa_pos"]) if r.get("aa_pos") else None,
        aa_ref=r.get("aa_ref") or None, aa_alt=r.get("aa_alt") or None,
        hi_score=hi.get(r["gene"]),
        mane=mane.get(r["gene"]),
        ps1_hit=r.get("ps1_hit") == "1", pm5_hit=r.get("pm5_hit") == "1",
    )
    out = {"variation_id": r["variation_id"], "gene": r["gene"], "vcep": r["vcep"],
           "tier_vcep": r["tier_vcep"], "consequence": r["consequence"],
           "codes_vcep": r["codes"], "tier_lookup": r.get("clinvar_tier") or "VUS"}
    for name, disabled in CONDITIONS.items():
        codes = evaluate(v, disabled)
        out[f"codes_{name}"] = ";".join(codes)
        out[f"tier_{name}"] = combine_2015(codes)
        out[f"pts_{name}"] = combine_points(codes)[1]
    results.append(out)

# ------------------------------------------------------------ assertions
for forbidden, cond in (("PP5", "svi"), ("BP6", "svi"), ("PP5", "blind"),
                        ("BP6", "blind"), ("PS1", "blind"), ("PM5", "blind")):
    n = sum(1 for r in results if forbidden in r[f"codes_{cond}"])
    assert n == 0, f"ASSERTION FAILED: {forbidden} fired {n} times under '{cond}'"
print("assertion OK: disabled criteria fire zero times")

# ------------------------------------------------------------ metrics
def concordance(key):
    n = sum(1 for r in results if r[key] == r["tier_vcep"])
    return n, n / len(results)

def per_class(key):
    out = {}
    for t in TIERS:
        tp = sum(1 for r in results if r[key] == t and r["tier_vcep"] == t)
        fp = sum(1 for r in results if r[key] == t and r["tier_vcep"] != t)
        fn = sum(1 for r in results if r[key] != t and r["tier_vcep"] == t)
        n = sum(1 for r in results if r["tier_vcep"] == t)
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        out[t] = (n, prec, rec, f1)
    return out

print("\n" + "=" * 78)
print("ARM B - tier concordance with ClinGen VCEP calls")
print("=" * 78)
order = [("tier_lookup", "ClinVar lookup baseline (no reasoning)"),
         ("tier_full", "acmgx full (incl. deprecated PP5/BP6)"),
         ("tier_svi", "acmgx SVI-compliant (no PP5/BP6)"),
         ("tier_blind", "acmgx ClinVar-blinded (no PS1/PM5 either)")]
for key, label in order:
    n, f = concordance(key)
    print(f"  {label:46s} {n:6d}/{len(results)}  {f:6.1%}")

print("\nper-class F1:")
hdr = "  {:<14s}" + "".join(f"{t:>10s}" for t in TIERS)
print(hdr.format("condition"))
for key, label in order:
    pc = per_class(key)
    print("  {:<14s}".format(key.replace("tier_", "")) +
          "".join(f"{pc[t][3]:10.2f}" for t in TIERS))
print("  {:<14s}".format("n") + "".join(f"{per_class('tier_svi')[t][0]:10d}" for t in TIERS))

print("\nVUS rate (fraction of variants called VUS):")
print(f"  expert panels           {sum(1 for r in results if r['tier_vcep']=='VUS')/len(results):6.1%}")
for key, label in order:
    print(f"  {key.replace('tier_',''):22s}  {sum(1 for r in results if r[key]=='VUS')/len(results):6.1%}")

# ------------------------------------------------- criterion-level agreement
print("\n" + "=" * 78)
print("CRITERION-LEVEL agreement with VCEP evidence codes (acmgx svi condition)")
print("=" * 78)
print(f"  {'criterion':<10s} {'VCEP':>7s} {'acmgx':>7s} {'TP':>6s} {'prec':>7s} {'rec':>7s} {'F1':>7s}")
crit_rows = []
# Metrics are computed against the `svi` calls, so a criterion disabled under svi
# can only ever score zero. Emitting it as "precision 0.00" reads as a failure of
# the engine rather than a condition of the experiment, so it is excluded and the
# condition is recorded in the file.
METRIC_CONDITION = "svi"
for c in sorted(set(CRITERIA) - CONDITIONS[METRIC_CONDITION]):
    vcep_set, ours = set(), set()
    for i, r in enumerate(results):
        vc = {parse_code(x)[0] for x in r["codes_vcep"].split(";") if x}
        oc = {parse_code(x)[0] for x in r["codes_svi"].split(";") if x}
        if c in vc: vcep_set.add(i)
        if c in oc: ours.add(i)
    tp = len(vcep_set & ours)
    prec = tp / len(ours) if ours else 0.0
    rec = tp / len(vcep_set) if vcep_set else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    print(f"  {c:<10s} {len(vcep_set):7d} {len(ours):7d} {tp:6d} {prec:7.2f} {rec:7.2f} {f1:7.2f}")
    crit_rows.append([METRIC_CONDITION, c, len(vcep_set), len(ours), tp,
                      round(prec, 4), round(rec, 4), round(f1, 4)])

print("\nCriteria acmgx does not implement (declared, not silently skipped):")
for c, why in sorted(NOT_IMPLEMENTED.items()):
    print(f"  {c:<6s} {why}")

# ------------------------------------------------------------ write
os.makedirs(f"{REPO}/results", exist_ok=True)
cols = list(results[0].keys())
with open(f"{REPO}/results/arm_b_calls.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader(); w.writerows(results)
with open(f"{REPO}/results/criterion_metrics.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["condition","criterion","n_vcep","n_acmgx","tp","precision","recall","f1"])
    w.writerows(crit_rows)
print(f"\nwrote {REPO}/results/arm_b_calls.tsv and criterion_metrics.tsv")
