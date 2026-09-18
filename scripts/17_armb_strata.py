#!/usr/bin/env python3
"""Arm B: stratify by VCEP, and place acmgx against the ceiling.

VCEP calls apply gene-specific criteria specifications. Scoring a generic engine
against them conflates "the rules are wrong" with "the rules were not specialised
for this gene". Stratifying separates the two; the aggregate number does not.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, os, sys
from collections import defaultdict, Counter
import numpy as np

sys.path.insert(0, _paths.s(_paths.REPO / "src"))
from acmgx.stats import wilson_ci

REPO = _paths.s(_paths.REPO)
calls = list(csv.DictReader(open(f"{REPO}/results/arm_b_calls.tsv"), delimiter="\t"))
ceil = {r["caid"]: r for r in csv.DictReader(open(f"{REPO}/results/ceiling.tsv"), delimiter="\t")}
print(f"{len(calls)} classified variants")

# ------------------------------------------------------- ceiling vs achieved
# join on the ERepo variant identity: ceiling.tsv is keyed by CAID, calls by
# ClinVar VariationID, so re-derive the ceiling per variant from its VCEP codes.
sys.path.insert(0, f"{REPO}/src")
from acmgx.combine import combine_2015, parse_code
from acmgx.assessability import bucket

def ceiling_tier(codes):
    keep = [c for c in codes if bucket(parse_code(c)[0]) in ("assessable", "knowledge")]
    return combine_2015(keep)

for r in calls:
    codes = [c for c in r["codes_vcep"].split(";") if c]
    r["tier_ceiling"] = ceiling_tier(codes)
    r["n_unassessable"] = sum(1 for c in codes if bucket(parse_code(c)[0]) == "unassessable")

def acc(rows, key):
    if not rows: return (0.0, 0.0, 0.0, 0)
    n = sum(1 for r in rows if r[key] == r["tier_vcep"])
    lo, hi = wilson_ci(n, len(rows))
    return (n / len(rows), lo, hi, len(rows))

print("\n" + "=" * 84)
print("CEILING vs ACHIEVED")
print("=" * 84)
for key, lab in (("tier_ceiling", "ceiling (VCEP evidence, unassessable removed)"),
                 ("tier_full", "acmgx full (reads ClinVar)"),
                 ("tier_svi", "acmgx SVI-compliant"),
                 ("tier_blind", "acmgx ClinVar-blinded")):
    a, lo, hi, n = acc(calls, key)
    print(f"  {lab:<48s} {a:6.1%} [{lo:.1%},{hi:.1%}]")
c_ceil = acc(calls, "tier_ceiling")[0]
c_svi = acc(calls, "tier_svi")[0]
print(f"\n  gap between the ceiling and a real SVI-compliant engine: {c_ceil - c_svi:.1%}")
print(f"  (the ceiling is a property of the evidence, not of any tool)")

# ------------------------------------------------------- restricted to reachable
reachable = [r for r in calls if r["n_unassessable"] == 0]
print(f"\nOn the {len(reachable)} variants whose VCEP call used NO unassessable evidence:")
for key, lab in (("tier_full", "acmgx full"), ("tier_svi", "acmgx SVI"),
                 ("tier_blind", "acmgx blinded")):
    a, lo, hi, n = acc(reachable, key)
    print(f"  {lab:<24s} {a:6.1%} [{lo:.1%},{hi:.1%}]")

# ------------------------------------------------------- by VCEP
print("\n" + "=" * 84)
print("BY EXPERT PANEL (panels with >= 150 variants)")
print("=" * 84)
by = defaultdict(list)
for r in calls:
    by[r["vcep"] or "unknown"].append(r)
print(f"  {'panel':<46s} {'n':>5s} {'lookup':>8s} {'svi':>8s} {'blind':>8s} {'ceiling':>8s}")
for v, rows in sorted(by.items(), key=lambda kv: -len(kv[1])):
    if len(rows) < 150:
        continue
    print(f"  {v[:44]:<46s} {len(rows):5d} "
          f"{acc(rows,'tier_lookup')[0]:8.1%} {acc(rows,'tier_svi')[0]:8.1%} "
          f"{acc(rows,'tier_blind')[0]:8.1%} {acc(rows,'tier_ceiling')[0]:8.1%}")
small = [r for v, rows in by.items() if len(rows) < 150 for r in rows]
if small:
    print(f"  {'(panels with < 150 variants, pooled)':<46s} {len(small):5d} "
          f"{acc(small,'tier_lookup')[0]:8.1%} {acc(small,'tier_svi')[0]:8.1%} "
          f"{acc(small,'tier_blind')[0]:8.1%} {acc(small,'tier_ceiling')[0]:8.1%}")

# ------------------------------------------------------- by consequence
print("\n" + "=" * 84)
print("BY MOLECULAR CONSEQUENCE (acmgx SVI condition)")
print("=" * 84)
byc = defaultdict(list)
for r in calls:
    byc[r["consequence"]].append(r)
print(f"  {'consequence':<26s} {'n':>5s} {'svi':>8s} {'ceiling':>8s}")
for c, rows in sorted(byc.items(), key=lambda kv: -len(kv[1]))[:9]:
    print(f"  {c:<26s} {len(rows):5d} {acc(rows,'tier_svi')[0]:8.1%} {acc(rows,'tier_ceiling')[0]:8.1%}")

with open(f"{REPO}/results/arm_b_strata.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["stratum_type", "stratum", "n", "lookup", "full", "svi", "blind", "ceiling"])
    for v, rows in sorted(by.items(), key=lambda kv: -len(kv[1])):
        w.writerow(["vcep", v, len(rows)] +
                   [round(acc(rows, k)[0], 4) for k in
                    ("tier_lookup", "tier_full", "tier_svi", "tier_blind", "tier_ceiling")])
    for c, rows in sorted(byc.items(), key=lambda kv: -len(kv[1])):
        w.writerow(["consequence", c, len(rows)] +
                   [round(acc(rows, k)[0], 4) for k in
                    ("tier_lookup", "tier_full", "tier_svi", "tier_blind", "tier_ceiling")])
print(f"\nwrote {REPO}/results/arm_b_strata.tsv")
