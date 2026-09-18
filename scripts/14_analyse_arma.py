#!/usr/bin/env python3
"""Arm A analysis: blinding deltas, with intervals on everything."""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, json, os, sys
from collections import defaultdict, Counter
import numpy as np

sys.path.insert(0, _paths.s(_paths.REPO / "src"))
from acmgx.stats import (bootstrap_ci, paired_bootstrap_diff, mcnemar_exact,
                         wilson_ci, mrr, top_n, min_detectable_effect)

REPO = _paths.s(_paths.REPO)
rows = list(csv.DictReader(open(f"{REPO}/results/arm_a_ranks.tsv"), delimiter="\t"))
print(f"{len(rows)} runs")

INF = float("inf")

# Exomiser ranks GENES: every variant in a gene carries that gene's combined
# score, so a variant-level rank counts a gene's other variants against the
# causal one. For compound heterozygotes that makes rank 1 unreachable outright.
# Gene rank is therefore the primary metric; the variant-level figure is
# reported beside it.
RANK_FIELD = "gene_rank_pess"

def rank_of(r, field=None):
    """Rank of the causal gene. Unranked if it did not survive the filter cascade."""
    field = field or RANK_FIELD
    if r["found"] != "1" or not r.get(field):
        return INF
    return float(r[field])

# index: (case_id, condition, tier, kind) -> row
idx = {(r["case_id"], r["condition"], r["hpo_tier"], r["kind"]): r for r in rows}
cases = sorted({r["case_id"] for r in rows})
models = {r["case_id"]: r["model"] for r in rows}

CONDS = ["informed", "clinvar_blind", "pheno_blind", "double_blind"]
ALL_CONDS = CONDS + ["ad_penetrant"]

def ranks_for(cond, tier="as_reported", kind="pos", subset=None):
    out, ids = [], []
    for c in (subset or cases):
        r = idx.get((c, cond, tier, kind))
        if r is None:
            continue
        out.append(rank_of(r)); ids.append(c)
    return np.array(out, float), ids

print("\n" + "=" * 82)
print("ARM A - causal variant rank by blinding condition (pessimistic tie policy)")
print("=" * 82)
print(f"  {'condition':<16s} {'n':>5s} {'top-1':>18s} {'top-10':>18s} {'MRR':>8s} {'median':>8s} {'survived':>8s}")
summary = {}
for cond in CONDS:
    a, ids = ranks_for(cond)
    if a.size == 0:
        continue
    t1, t1lo, t1hi = bootstrap_ci(a, lambda x: top_n(x, 1))
    t10, t10lo, t10hi = bootstrap_ci(a, lambda x: top_n(x, 10))
    m = mrr(a)
    med = float(np.median(a[np.isfinite(a)])) if np.isfinite(a).any() else float("nan")
    nfound = sum(1 for c in ids if idx[(c, cond, "as_reported", "pos")]["found"] == "1")
    summary[cond] = dict(n=a.size, top1=t1, top10=t10, mrr=m, median=med,
                         survived=nfound / a.size, ranks=a, ids=ids)
    print(f"  {cond:<16s} {a.size:5d} "
          f"{t1:6.3f} [{t1lo:.3f},{t1hi:.3f}] "
          f"{t10:6.3f} [{t10lo:.3f},{t10hi:.3f}] "
          f"{m:8.3f} {med:8.1f} {nfound/a.size:8.3f}")

print(f"\n  minimum detectable difference in top-1 at n={summary['informed']['n']}: "
      f"{min_detectable_effect(summary['informed']['n']):.3f}")

print("\n  tie policy sensitivity and metric level (informed condition):")
ri = [idx[(c, "informed", "as_reported", "pos")] for c in cases
      if (c, "informed", "as_reported", "pos") in idx]
for lbl, fld in (("gene, pessimistic", "gene_rank_pess"), ("gene, optimistic", "gene_rank_opt"),
                 ("variant, pessimistic", "rank_pess"), ("variant, optimistic", "rank_opt")):
    a = np.array([rank_of(r, fld) for r in ri])
    print(f"    {lbl:<22s} top-1={top_n(a,1):.3f}  top-10={top_n(a,10):.3f}")

# --------------------------------------------------- paired blinding deltas
print("\n" + "=" * 82)
print("BLINDING DELTAS (paired: same cases, same background, one channel removed)")
print("=" * 82)
def paired(a_cond, b_cond, label):
    common = [c for c in cases
              if (c, a_cond, "as_reported", "pos") in idx and (c, b_cond, "as_reported", "pos") in idx]
    A = np.array([rank_of(idx[(c, a_cond, "as_reported", "pos")]) for c in common])
    B = np.array([rank_of(idx[(c, b_cond, "as_reported", "pos")]) for c in common])
    for n in (1, 10):
        d, lo, hi = paired_bootstrap_diff(A, B, lambda x: top_n(x, n))
        # McNemar on the discordant pairs
        b = int(np.sum((A <= n) & (B > n))); c_ = int(np.sum((A > n) & (B <= n)))
        p = mcnemar_exact(b, c_)
        print(f"  {label:<38s} top-{n:<2d} delta={d:+.3f} [{lo:+.3f},{hi:+.3f}]  "
              f"discordant {b}/{c_}  McNemar p={p:.2e}")
    return common

paired("informed", "clinvar_blind", "ClinVar informed - blinded")
paired("informed", "pheno_blind",  "phenotype informed - blinded")
paired("informed", "double_blind", "informed - double blinded")

# --------------------------------------------------- whitelist effect
print("\n" + "=" * 82)
print("CLINVAR WHITELIST - which causal variants were looked up, and what it bought")
print("=" * 82)
wl, nwl = [], []
for c in cases:
    ri = idx.get((c, "informed", "as_reported", "pos"))
    rb = idx.get((c, "clinvar_blind", "as_reported", "pos"))
    if not ri or not rb:
        continue
    (wl if ri.get("whitelisted") == "1" else nwl).append(
        (rank_of(ri), rank_of(rb)))
for label, grp in (("whitelisted causal variants", wl), ("not whitelisted", nwl)):
    if not grp:
        print(f"  {label:<32s} n=0")
        continue
    A = np.array([g[0] for g in grp]); B = np.array([g[1] for g in grp])
    d, lo, hi = paired_bootstrap_diff(A, B, lambda x: top_n(x, 1))
    print(f"  {label:<32s} n={len(grp):5d}  top-1 informed={top_n(A,1):.3f} "
          f"blinded={top_n(B,1):.3f}  delta={d:+.3f} [{lo:+.3f},{hi:+.3f}]")

# --------------------------------------------------- by inheritance model
print("\n" + "=" * 82)
print("BY INHERITANCE MODEL (informed condition)")
print("=" * 82)
print(f"  {'model':<14s} {'n':>5s} {'survived filters':>17s} {'top-1':>18s} {'top-10':>18s}")
by_model, surv = defaultdict(list), defaultdict(list)
for c in cases:
    r = idx.get((c, "informed", "as_reported", "pos"))
    if r:
        by_model[r["model"]].append(rank_of(r))
        surv[r["model"]].append(1.0 if r["found"] == "1" else 0.0)
for m, v in sorted(by_model.items()):
    a = np.array(v)
    sv = float(np.mean(surv[m]))
    slo, shi = wilson_ci(int(sum(surv[m])), len(surv[m]))
    t1, lo, hi = bootstrap_ci(a, lambda x: top_n(x, 1))
    t10, lo10, hi10 = bootstrap_ci(a, lambda x: top_n(x, 10))
    print(f"  {m:<14s} {a.size:5d} {sv:6.3f} [{slo:.2f},{shi:.2f}] "
          f"{t1:6.3f} [{lo:.3f},{hi:.3f}] {t10:6.3f} [{lo10:.3f},{hi10:.3f}]")

# ------------------------------- what the inheritance filter costs under
# incomplete penetrance: the same cases with the transmitting parent coded
# affected instead of unaffected.
pen = [c for c in cases if (c, "ad_penetrant", "as_reported", "pos") in idx]
if pen:
    A = np.array([rank_of(idx[(c, "ad_penetrant", "as_reported", "pos")]) for c in pen])
    B = np.array([rank_of(idx[(c, "informed", "as_reported", "pos")]) for c in pen])
    sa = np.mean([idx[(c, "ad_penetrant", "as_reported", "pos")]["found"] == "1" for c in pen])
    sb = np.mean([idx[(c, "informed", "as_reported", "pos")]["found"] == "1" for c in pen])
    print(f"\n  inherited-dominant cases, transmitting parent coded AFFECTED (n={len(pen)}):")
    print(f"    survived filters {sa:.3f} (vs {sb:.3f} when coded unaffected)")
    d, lo, hi = paired_bootstrap_diff(A, B, lambda x: top_n(x, 1))
    print(f"    top-1 {top_n(A,1):.3f} vs {top_n(B,1):.3f}   delta={d:+.3f} [{lo:+.3f},{hi:+.3f}]")
    d, lo, hi = paired_bootstrap_diff(A, B, lambda x: top_n(x, 10))
    print(f"    top-10 {top_n(A,10):.3f} vs {top_n(B,10):.3f}  delta={d:+.3f} [{lo:+.3f},{hi:+.3f}]")

# --------------------------------------------------- HPO perturbation
print("\n" + "=" * 82)
print("HPO PERTURBATION (informed condition, paired within case)")
print("=" * 82)
base_ids = [c for c in cases if (c, "informed", "imprecise", "pos") in idx]
if base_ids:
    A0 = np.array([rank_of(idx[(c, "informed", "as_reported", "pos")]) for c in base_ids])
    print(f"  {'tier':<14s} {'n':>5s} {'top-1':>7s} {'top-10':>7s} {'delta top-1 vs as_reported':>32s}")
    print(f"  {'as_reported':<14s} {A0.size:5d} {top_n(A0,1):7.3f} {top_n(A0,10):7.3f}")
    for tier in ("imprecise", "noisy", "truncated"):
        B = np.array([rank_of(idx[(c, "informed", tier, "pos")]) for c in base_ids
                      if (c, "informed", tier, "pos") in idx])
        if B.size != A0.size:
            A = np.array([rank_of(idx[(c, "informed", "as_reported", "pos")]) for c in base_ids
                          if (c, "informed", tier, "pos") in idx])
        else:
            A = A0
        d, lo, hi = paired_bootstrap_diff(B, A, lambda x: top_n(x, 1))
        print(f"  {tier:<14s} {B.size:5d} {top_n(B,1):7.3f} {top_n(B,10):7.3f} "
              f"{d:+18.3f} [{lo:+.3f},{hi:+.3f}]")

# --------------------------------------------------- negative controls
print("\n" + "=" * 82)
print("MATCHED KNOCKOUTS - what the pipeline does when there is nothing to find")
print("=" * 82)
for cond in ("informed", "clinvar_blind"):
    pos_top, ko_top, causal = [], [], []
    for c in cases:
        rp = idx.get((c, cond, "as_reported", "pos"))
        rk = idx.get((c, cond, "as_reported", "ko"))
        if not rp or not rk:
            continue
        try:
            ko_top.append(float(rk["top_combined"])); pos_top.append(float(rp["top_combined"]))
            causal.append(float(rp["combined"]) if rp["combined"] else 0.0)
        except ValueError:
            continue
    if not ko_top:
        continue
    ko = np.array(ko_top); cz = np.array(causal)
    # how often does the causal variant beat the best hit of its own knockout?
    beats = float(np.mean(cz > ko))
    lo, hi = wilson_ci(int(np.sum(cz > ko)), cz.size)
    print(f"  {cond:<14s} n={cz.size:4d}  median top score in knockout={np.median(ko):.4f}  "
          f"median causal score={np.median(cz):.4f}")
    print(f"  {'':14s}        causal beats its own knockout's best hit: "
          f"{beats:.3f} [{lo:.3f},{hi:.3f}]")

# --------------------------------------------------- candidate set size
cand = [int(r["n_candidates"]) for r in rows if r["kind"] == "pos" and r["condition"] == "informed"]
gen = [int(r["n_genes"]) for r in rows if r["kind"] == "pos" and r["condition"] == "informed" and r["n_genes"]]
if cand:
    print(f"\ncandidate set after the filter cascade: median {int(np.median(cand))} variants, "
          f"{int(np.median(gen))} genes (IQR {int(np.percentile(cand,25))}-{int(np.percentile(cand,75))})")

json.dump({k: {kk: (vv if not isinstance(vv, np.ndarray) else None)
               for kk, vv in v.items() if kk not in ("ranks", "ids")}
           for k, v in summary.items()},
          open(f"{REPO}/results/arm_a_summary.json", "w"), indent=1)
print(f"\nwrote {REPO}/results/arm_a_summary.json")
