#!/usr/bin/env python3
"""Extract causal-variant and causal-gene ranks from every Exomiser run.

Tie policy is explicit and pessimistic: where the causal gene ties with others
on the combined score, its rank is the worst position in the tie block. The
optimistic and mean-of-block variants are recorded alongside, so the effect of
the policy is visible rather than buried in a tie-break the tool chose.
"""
import csv, json, os, sys, glob
from collections import defaultdict

REPO = os.path.expanduser("~/rdacmg/repo")
ARMA = os.path.expanduser("~/rdacmg/arma")

manifest = {c["case_id"]: c for c in json.load(open(f"{REPO}/data/manifest.json"))}
print(f"{len(manifest)} cases in manifest")

files = sorted(glob.glob(f"{ARMA}/*.variants.tsv"))
print(f"{len(files)} Exomiser variant tables")

rows = []
for path in files:
    name = os.path.basename(path)[: -len(".variants.tsv")]
    try:
        case_id, cond, tier, kind = name.split("__")
    except ValueError:
        continue
    case = manifest.get(case_id)
    if case is None:
        continue

    with open(path) as fh:
        rd = list(csv.DictReader(fh, delimiter="\t"))
    if not rd:
        rows.append(dict(case_id=case_id, condition=cond, hpo_tier=tier, kind=kind,
                         model=case["model"], gene=case["gene"], n_candidates=0,
                         found=0, rank_pess="", rank_opt="", gene_rank_pess="",
                         combined="", variant_score="", pheno_score="", whitelisted="",
                         top_combined="", acmg=""))
        continue

    def f(r, k, default=0.0):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return default

    # gene-level scores, one per gene
    gene_score = {}
    for r in rd:
        g = r.get("GENE_SYMBOL") or ""
        gene_score[g] = max(gene_score.get(g, -1.0), f(r, "EXOMISER_GENE_COMBINED_SCORE"))
    top_combined = max(gene_score.values()) if gene_score else 0.0

    truth = {(t["chrom"], t["pos"], t["ref"], t["alt"]) for t in case["truth"]}
    hit = None
    for r in rd:
        key = ("chr" + r["CONTIG"] if not r["CONTIG"].startswith("chr") else r["CONTIG"],
               int(r["START"]), r["REF"], r["ALT"])
        if key in truth:
            if hit is None or f(r, "EXOMISER_GENE_COMBINED_SCORE") > f(hit, "EXOMISER_GENE_COMBINED_SCORE"):
                hit = r

    rec = dict(case_id=case_id, condition=cond, hpo_tier=tier, kind=kind,
               model=case["model"], gene=case["gene"], n_candidates=len(rd),
               n_genes=len(gene_score), top_combined=round(top_combined, 6))

    if hit is None:
        rec.update(found=0, rank_pess="", rank_opt="", gene_rank_pess="",
                   combined="", variant_score="", pheno_score="", whitelisted="", acmg="")
    else:
        cs = f(hit, "EXOMISER_GENE_COMBINED_SCORE")
        # A compound-heterozygous case has two causal alleles that necessarily
        # share the gene's combined score. Counting the partner allele against
        # the variant makes rank 1 unreachable for every comp-het case, which is
        # a property of the metric and not of the tool. Exclude a case's own
        # other causal alleles from the tie block.
        def is_partner(r):
            key = ("chr" + r["CONTIG"] if not r["CONTIG"].startswith("chr") else r["CONTIG"],
                   int(r["START"]), r["REF"], r["ALT"])
            return key in truth and r is not hit
        vs_all = [f(r, "EXOMISER_GENE_COMBINED_SCORE") for r in rd if not is_partner(r)]
        rank_pess = sum(1 for s in vs_all if s >= cs)
        rank_opt = sum(1 for s in vs_all if s > cs) + 1
        gs = sorted(gene_score.values(), reverse=True)
        grank_pess = sum(1 for s in gs if s >= cs)
        grank_opt = sum(1 for s in gs if s > cs) + 1
        rec.update(found=1, rank_pess=rank_pess, rank_opt=rank_opt,
                   gene_rank_pess=grank_pess, gene_rank_opt=grank_opt,
                   combined=round(cs, 6),
                   variant_score=round(f(hit, "EXOMISER_VARIANT_SCORE"), 6),
                   pheno_score=round(f(hit, "EXOMISER_GENE_PHENO_SCORE"), 6),
                   whitelisted=hit.get("WHITELIST_VARIANT", ""),
                   acmg=hit.get("EXOMISER_ACMG_CLASSIFICATION", ""))
    rows.append(rec)

cols = ["case_id", "condition", "hpo_tier", "kind", "model", "gene", "found",
        "rank_pess", "rank_opt", "gene_rank_pess", "gene_rank_opt",
        "combined", "variant_score", "pheno_score", "whitelisted", "acmg",
        "n_candidates", "n_genes", "top_combined"]
os.makedirs(f"{REPO}/results", exist_ok=True)
with open(f"{REPO}/results/arm_a_ranks.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})

from collections import Counter
print("\nruns by condition/kind:", Counter((r["condition"], r["kind"]) for r in rows).most_common())
pos = [r for r in rows if r["kind"] == "pos"]
print(f"positives with causal variant found: "
      f"{sum(r['found'] for r in pos)}/{len(pos)}")
print(f"\nwrote {REPO}/results/arm_a_ranks.tsv ({len(rows)} rows)")
