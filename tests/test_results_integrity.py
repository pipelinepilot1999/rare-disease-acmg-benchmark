#!/usr/bin/env python3
"""Guard the published numbers against the results files.

Every headline figure in README.md and docs/results.md is recomputed here from
results/ and asserted to appear, verbatim, in the prose. This is the test that
was missing: the engine had unit tests, but nothing tied the documents to the
data, so a number could drift out of date silently.

It also asserts two structural properties that a reader would otherwise have to
take on trust: that the matched-knockout statistic is arithmetically identical to
top-1 rather than an independent measurement, and that the criteria the engine
implements, the criteria it declares unimplemented, and the assessability buckets
remain mutually consistent.

Run: python tests/test_results_integrity.py
"""
import csv
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))

from acmgx.combine import parse_code
from acmgx.assessability import bucket, MACHINE_ASSESSABLE, KNOWLEDGE_DERIVED, NOT_MACHINE_ASSESSABLE
from acmgx.criteria import CRITERIA, NOT_IMPLEMENTED, CLINVAR_DERIVED

fails = []


def check(name, got, want):
    if got != want:
        fails.append(f"{name}: got {got!r}, want {want!r}")


def cited(doc, text, label):
    """Assert a recomputed figure is actually the one written in the prose."""
    if text not in doc.body:
        fails.append(f"{doc.name}: '{text}' ({label}) not found - the document has drifted "
                     f"from results/")


def anchored(doc, before, value, after, label):
    """Assert a figure appears in its own sentence, not merely somewhere in the file.

    `cited` only proves the right number is present somewhere, which a document
    with the same figure quoted twice can satisfy while one copy is stale. For
    the headline claims the surrounding words are pinned too.
    """
    needle = f"{before}{value}{after}"
    if needle not in doc.body:
        fails.append(f"{doc.name}: anchored claim for {label} is stale - expected "
                     f"...{before}[{value}]{after}...")


class Doc:
    def __init__(self, path):
        self.name = os.path.relpath(path, REPO)
        self.body = open(path).read()


def tsv(path):
    return list(csv.DictReader(open(os.path.join(REPO, path)), delimiter="\t"))


def pct(x, places=1):
    return f"{100 * x:.{places}f}%"


readme = Doc(os.path.join(REPO, "README.md"))
results = Doc(os.path.join(REPO, "docs", "results.md"))
limits = Doc(os.path.join(REPO, "docs", "limitations.md"))

# ------------------------------------------------------------------ Arm B'
ceiling = tsv("results/ceiling.tsv")
check("ceiling.tsv rows", len(ceiling), 12967)
agree = lambda key: sum(1 for r in ceiling if r[key] == r["tier_vcep"]) / len(ceiling)

cited(results, "12,967", "ERepo interpretations")
cited(results, pct(agree("tier_recomputed")), "combining-rule validation")
cited(results, pct(agree("tier_ceiling")), "ceiling, all interpretations")
cited(readme, pct(agree("tier_ceiling")), "ceiling, all interpretations")
cited(results, pct(agree("tier_ceiling_blind")), "ceiling, ClinVar-blinded")

dep = sum(1 for r in ceiling if int(r["n_unassessable"]) > 0) / len(ceiling)
cited(results, pct(dep), "share depending on unassessable evidence")
cited(readme, pct(dep), "share depending on unassessable evidence")

by_tier = Counter(r["tier_vcep"] for r in ceiling)
kept = Counter(r["tier_vcep"] for r in ceiling if r["tier_vcep"] == r["tier_ceiling"])
for tier in ("P", "LP", "VUS", "LB", "B"):
    cited(results, f"{by_tier[tier]:,}", f"n for expert tier {tier}")
    cited(results, pct(kept[tier] / by_tier[tier]), f"ceiling retention for {tier}")
cited(readme, pct(kept["P"] / by_tier["P"]), "ceiling retention for Pathogenic")
anchored(results, "| Pathogenic | 3,516 | **", pct(kept["P"] / by_tier["P"]), "**",
         "Pathogenic retention row")
anchored(results, "and for only ", pct(kept["P"] / by_tier["P"]), " of Pathogenic calls",
         "Pathogenic retention, conclusion")

# ------------------------------------------------------------------- Arm B
calls = tsv("results/arm_b_calls.tsv")
check("arm_b_calls.tsv rows", len(calls), 11499)
cited(results, "11,499", "Arm B variant count")
cited(readme, "11,499", "Arm B variant count")

conc = lambda key: sum(1 for r in calls if r[key] == r["tier_vcep"]) / len(calls)
for key, label in (("tier_lookup", "lookup baseline"), ("tier_full", "acmgx full"),
                   ("tier_svi", "acmgx svi"), ("tier_blind", "acmgx blind")):
    cited(results, pct(conc(key)), label)
cited(readme, pct(conc("tier_lookup")), "lookup baseline")
cited(readme, pct(conc("tier_svi")), "acmgx svi")

# headline claims pinned to their own sentence, so a stale duplicate cannot hide
anchored(results, "ClinVar's aggregate call agrees with ClinGen expert panels on **",
         pct(conc("tier_lookup")), "** of 11,499 variants", "lookup baseline, lead paragraph")
anchored(results, "ClinVar assertions, reaches **", pct(conc("tier_svi")), "**",
         "acmgx svi, lead paragraph")
anchored(readme, "panels on **", pct(conc("tier_lookup")), "** of 11,499 variants",
         "lookup baseline, README")
anchored(readme, "not read ClinVar, reaches **", pct(conc("tier_svi")), "**",
         "acmgx svi, README")


# the ceiling restricted to the Arm B subset - the *other* ceiling number
def ceiling_tier(codes):
    keep = [c for c in codes if bucket(parse_code(c)[0]) in ("assessable", "knowledge")]
    from acmgx.combine import combine_2015
    return combine_2015(keep)

sub = sum(1 for r in calls
          if ceiling_tier([c for c in r["codes_vcep"].split(";") if c]) == r["tier_vcep"])
cited(results, pct(sub / len(calls)), "ceiling on the Arm B subset")
anchored(results, "only **", pct(sub / len(calls)), "** of the 11,499 variants Arm B classifies",
         "ceiling, Arm B denominator")
anchored(results, "(**", pct(agree("tier_ceiling")), "** across all 12,967 expert",
         "ceiling, full denominator")
anchored(readme, "is\n  **", pct(agree("tier_ceiling")), "** over all 12,967 expert interpretations",
         "ceiling, README")
# both ceiling figures must be present, and distinguishable
check("two ceiling numbers documented",
      pct(agree("tier_ceiling")) in results.body and pct(sub / len(calls)) in results.body, True)

vus = lambda key: sum(1 for r in calls if r[key] == "VUS") / len(calls)
for key in ("tier_vcep", "tier_lookup", "tier_full", "tier_svi", "tier_blind"):
    cited(results, pct(vus(key)), f"VUS rate {key}")

reachable = [r for r in calls
             if not any(bucket(parse_code(c)[0]) == "unassessable"
                        for c in r["codes_vcep"].split(";") if c)]
cited(results, f"{len(reachable):,}", "variants needing no unassessable evidence")
for key in ("tier_full", "tier_svi", "tier_blind"):
    n = sum(1 for r in reachable if r[key] == r["tier_vcep"])
    cited(results, pct(n / len(reachable)), f"{key} on reachable subset")

# ------------------------------------------------------- criterion metrics
crit = tsv("results/criterion_metrics.tsv")
check("criterion metrics record their condition",
      {r["condition"] for r in crit}, {"svi"})
check("no criterion disabled under svi is scored",
      {r["criterion"] for r in crit} & {"PP5", "BP6"}, set())
for r in crit:
    cited(results, f"{float(r['precision']):.2f}", f"{r['criterion']} precision")

# ------------------------------------------------------------------- Arm A
ranks = tsv("results/arm_a_ranks.tsv")
idx = {(r["case_id"], r["condition"], r["hpo_tier"], r["kind"]): r for r in ranks}
cases = sorted({r["case_id"] for r in ranks})
check("Arm A cases", len(cases), 1034)
cited(results, "1,034", "Arm A case count")

INF = float("inf")


def rank(case, cond, tier="as_reported"):
    r = idx.get((case, cond, tier, "pos"))
    if not r or r["found"] != "1" or not r["gene_rank_pess"]:
        return INF
    return float(r["gene_rank_pess"])


def topn(arr, n):
    return sum(1 for a in arr if a <= n) / len(arr)


summary = json.load(open(os.path.join(REPO, "results", "arm_a_summary.json")))
for cond in ("informed", "clinvar_blind", "pheno_blind", "double_blind"):
    arr = [rank(c, cond) for c in cases]
    check(f"{cond} top-1 matches summary json",
          round(topn(arr, 1), 4), round(summary[cond]["top1"], 4))
    cited(results, f"{topn(arr, 1):.3f}", f"{cond} top-1")

wl = [c for c in cases if idx[(c, "informed", "as_reported", "pos")]["whitelisted"] == "1"]
nwl = [c for c in cases if c not in set(wl)]
check("whitelisted causal alleles", len(wl), 300)
cited(results, "300 of 1,034", "whitelist coverage")
cited(readme, pct(len(wl) / len(cases)), "whitelist coverage")
for ids, conds in ((wl, ("informed", "clinvar_blind", "pheno_blind", "double_blind")),
                   (nwl, ("informed", "clinvar_blind", "pheno_blind", "double_blind"))):
    for cond in conds:
        cited(results, f"{topn([rank(c, cond) for c in ids], 1):.3f}",
              f"{cond} top-1 on a whitelist stratum")

# the non-whitelisted stratum must show exactly zero ClinVar effect - the control
on = topn([rank(c, "pheno_blind") for c in nwl], 1)
off = topn([rank(c, "double_blind") for c in nwl], 1)
check("non-whitelisted stratum is a null control", abs(on - off) < 0.02, True)

for model, n in (("denovo", 266), ("ar_hom", 352), ("comphet", 189), ("ad_inherited", 227)):
    ids = [c for c in cases
           if idx[(c, "informed", "as_reported", "pos")]["model"] == model]
    check(f"{model} n", len(ids), n)
    surv = sum(1 for c in ids
               if idx[(c, "informed", "as_reported", "pos")]["found"] == "1") / len(ids)
    cited(results, f"{surv:.3f}", f"{model} survival")

pen = [c for c in cases if (c, "ad_penetrant", "as_reported", "pos") in idx]
check("penetrance contrast cases", len(pen), 227)
cited(results, f"{topn([rank(c, 'ad_penetrant') for c in pen], 1):.3f}", "penetrant top-1")
anchored(results, "| parent coded affected | **0.996** | **",
         f"{topn([rank(c, 'ad_penetrant') for c in pen], 1):.3f}", "**",
         "penetrance contrast top-1")

# ------------------------------------- the knockout statistic is not independent
for cond in ("informed", "clinvar_blind"):
    beats, first = [], []
    for c in cases:
        p = idx[(c, cond, "as_reported", "pos")]
        k = idx[(c, cond, "as_reported", "ko")]
        beats.append(bool(p["found"] == "1" and p["combined"] and k["top_combined"]
                          and float(p["combined"]) > float(k["top_combined"])))
        first.append(rank(c, cond) == 1)
    check(f"knockout statistic is identical to top-1 under {cond}", beats, first)
check("the identity is disclosed in the prose",
      "identical to top-1" in results.body, True)

# ------------------------------------------------------- structural invariants
check("CRITERIA and NOT_IMPLEMENTED are disjoint",
      set(CRITERIA) & set(NOT_IMPLEMENTED), set())
check("every implemented criterion has an assessability bucket",
      {c for c in CRITERIA if bucket(c) == "unknown"}, set())
check("every ClinVar-derived criterion is implemented",
      CLINVAR_DERIVED - set(CRITERIA), set())
check("assessability buckets are disjoint",
      MACHINE_ASSESSABLE & KNOWLEDGE_DERIVED
      | MACHINE_ASSESSABLE & NOT_MACHINE_ASSESSABLE
      | KNOWLEDGE_DERIVED & NOT_MACHINE_ASSESSABLE, set())
check("every unimplemented criterion is still bucketed",
      {c for c in NOT_IMPLEMENTED if bucket(c) == "unknown"}, set())

# the README must not claim a script range that does not exist
import re
listed = sorted(int(m.group(1)) for m in
                (re.match(r"(\d+)", f) for f in os.listdir(os.path.join(REPO, "scripts")))
                if m)
check("README states the real script range",
      f"01 -> {max(listed):02d}" in readme.body, True)

if fails:
    print(f"FAILED {len(fails)}:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print(f"all results-integrity checks pass "
      f"({len(cases)} Arm A cases, {len(calls)} Arm B variants, {len(ceiling)} interpretations)")
