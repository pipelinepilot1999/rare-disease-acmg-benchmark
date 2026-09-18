#!/usr/bin/env python3
"""Figures for docs/results.md.

Palette is the validated categorical set (slots 1-4), which passes the CVD and
normal-vision floors on the adjacent pairlist. Two of the four sit below 3:1
contrast on the light surface, so every chart carries visible direct labels -
the relief rule, not an optional nicety.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import csv, json, os, sys
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, _paths.s(_paths.REPO / "src"))
from acmgx.stats import bootstrap_ci, top_n, wilson_ci

REPO = _paths.s(_paths.REPO)
FIG = f"{REPO}/results/figures"
os.makedirs(FIG, exist_ok=True)

S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURFACE = "#fcfcfb"
GRID = "#e5e4e0"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f4f8fd", "#2a78d6", "#123356"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 9,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.6,
})

def style(ax, title=None, sub=None):
    """Title above subtitle above the plot, with enough room that they never touch."""
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=22 if sub else 8)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, color=MUTED, fontsize=8, va="bottom")
    ax.tick_params(length=0)

INF = float("inf")
def rank_of(r):
    # gene-level, pessimistic tie policy - see docs/decisions.md D-019
    return INF if r["found"] != "1" or not r["gene_rank_pess"] else float(r["gene_rank_pess"])

# ---------------------------------------------------------------- Arm A
ranks_path = f"{REPO}/results/arm_a_ranks.tsv"
have_a = os.path.exists(ranks_path)
if have_a:
    rows = list(csv.DictReader(open(ranks_path), delimiter="\t"))
    idx = {(r["case_id"], r["condition"], r["hpo_tier"], r["kind"]): r for r in rows}
    cases = sorted({r["case_id"] for r in rows})
    CONDS = ["informed", "clinvar_blind", "pheno_blind", "double_blind"]
    LABELS = {"informed": "Informed", "clinvar_blind": "ClinVar blinded",
              "pheno_blind": "Phenotype blinded", "double_blind": "Double blinded"}

    def arr(cond, tier="as_reported", kind="pos", ids=None):  # noqa: E306
        out = []
        for c in (ids or cases):
            r = idx.get((c, cond, tier, kind))
            if r: out.append(rank_of(r))
        return np.array(out, float)

    # --- Figure 1: blinding deltas
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    x = np.arange(len(CONDS)); w = 0.36
    for j, (n, col, lab) in enumerate([(1, S1, "top-1"), (10, S3, "top-10")]):
        vals, los, his = [], [], []
        for cond in CONDS:
            a = arr(cond)
            v, lo, hi = bootstrap_ci(a, lambda z: top_n(z, n)) if a.size else (0, 0, 0)
            vals.append(v); los.append(v - lo); his.append(hi - v)
        b = ax.bar(x + (j - 0.5) * w, vals, w * 0.92, color=col, label=lab,
                   edgecolor=SURFACE, linewidth=2)
        ax.errorbar(x + (j - 0.5) * w, vals, yerr=[los, his], fmt="none",
                    ecolor=INK2, elinewidth=1.2, capsize=3)
        for xi, v in zip(x + (j - 0.5) * w, vals):
            ax.text(xi, v + 0.035, f"{v:.2f}", ha="center", color=INK, fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels([LABELS[c] for c in CONDS])
    ax.set_ylim(0, 1.08); ax.set_ylabel("recall of causal gene")
    ax.yaxis.grid(True); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")
    style(ax, "What survives blinding",
          f"causal-gene recall, {len(cases)} matched cases, pessimistic tie policy, 95% bootstrap CI")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig1_blinding.png", dpi=170); plt.close(fig)

    # --- Figure 2: rank CDF
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for cond, col in zip(CONDS, [S1, S2, S3, S4]):
        a = arr(cond); a = a[np.isfinite(a)]
        if not a.size: continue
        xs = np.sort(a); ys = np.arange(1, xs.size + 1) / len(arr(cond))
        ax.step(xs, ys, where="post", color=col, linewidth=2, label=LABELS[cond])
        ax.text(xs[-1], ys[-1], f" {LABELS[cond]}", color=INK2, fontsize=8, va="center")
    ax.set_xscale("log"); ax.set_xlabel("rank of causal gene (log)")
    ax.set_ylabel("fraction of cases at or above rank")
    ax.set_ylim(0, 1.02); ax.yaxis.grid(True); ax.set_axisbelow(True)
    style(ax, "Rank distribution", "the whole curve, not just the top-N points")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig2_rank_cdf.png", dpi=170); plt.close(fig)

    # --- Figure 3: HPO tiers
    tiers = ["as_reported", "imprecise", "noisy", "truncated"]
    sub = [c for c in cases if (c, "informed", "imprecise", "pos") in idx]
    if sub:
        fig, ax = plt.subplots(figsize=(7.2, 3.8))
        vals, los, his = [], [], []
        for t in tiers:
            a = arr("informed", t, "pos", sub)
            v, lo, hi = bootstrap_ci(a, lambda z: top_n(z, 1)) if a.size else (0, 0, 0)
            vals.append(v); los.append(v - lo); his.append(hi - v)
        xb = np.arange(len(tiers))
        ax.bar(xb, vals, 0.55, color=S1, edgecolor=SURFACE, linewidth=2)
        ax.errorbar(xb, vals, yerr=[los, his], fmt="none", ecolor=INK2, elinewidth=1.2, capsize=3)
        for xi, v in zip(xb, vals):
            ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", color=INK, fontsize=9)
        ax.set_xticks(xb); ax.set_xticklabels(["as reported", "imprecise", "noisy", "truncated"])
        ax.set_ylabel("top-1 recall"); ax.yaxis.grid(True); ax.set_axisbelow(True)
        style(ax, "Phenotype quality", f"same {len(sub)} cases under four perturbations of the reported HPO terms")
        fig.tight_layout(); fig.savefig(f"{FIG}/fig3_hpo_tiers.png", dpi=170); plt.close(fig)

    # --- Figure 4: knockouts
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    cz, ko = [], []
    for c in cases:
        rp, rk = idx.get((c, "informed", "as_reported", "pos")), idx.get((c, "informed", "as_reported", "ko"))
        if not rp or not rk or not rp["combined"]: continue
        try:
            cz.append(float(rp["combined"])); ko.append(float(rk["top_combined"]))
        except ValueError: pass
    if cz:
        cz, ko = np.array(cz), np.array(ko)
        ax.scatter(ko, cz, s=14, color=S1, alpha=0.45, edgecolor="none")
        lim = [0, 1.02]
        ax.plot(lim, lim, color=MUTED, linewidth=1, linestyle="--")
        frac = float(np.mean(cz > ko))
        ax.text(0.04, 0.95, f"causal beats its own knockout's best hit\nin {frac:.1%} of cases",
                transform=ax.transAxes, color=INK, fontsize=9, va="top")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("best score in the matched knockout"); ax.set_ylabel("causal variant score")
        ax.grid(True); ax.set_axisbelow(True)
        style(ax, "Is the top hit meaningful?", "each point is one case and its own variant-removed control")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig4_knockout.png", dpi=170); plt.close(fig)

    # --- Figure 8: the 2x2, split by whether ClinVar could look the variant up
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=True)
    wl_ids = [c for c in cases
              if idx.get((c, "informed", "as_reported", "pos"), {}).get("whitelisted") == "1"]
    nwl_ids = [c for c in cases if c not in set(wl_ids)]
    for ax, ids, name in ((axes[0], wl_ids, "on Exomiser's ClinVar whitelist"),
                          (axes[1], nwl_ids, "not on the whitelist")):
        xb = np.arange(2); w = 0.36
        for j, (cv_on, cv_off, col, lab) in enumerate(
                [("informed", "clinvar_blind", S1, "phenotype supplied"),
                 ("pheno_blind", "double_blind", S2, "variant-only")]):
            vals = [top_n(arr(cv_on, ids=ids), 1), top_n(arr(cv_off, ids=ids), 1)]
            ax.bar(xb + (j - 0.5) * w, vals, w * 0.92, color=col,
                   edgecolor=SURFACE, linewidth=2, label=lab if ax is axes[0] else None)
            for xi, v in zip(xb + (j - 0.5) * w, vals):
                ax.text(xi, v + 0.025, f"{v:.2f}", ha="center", color=INK, fontsize=9)
        ax.set_xticks(xb); ax.set_xticklabels(["ClinVar on", "ClinVar blinded"])
        ax.set_ylim(0, 1.06); ax.yaxis.grid(True); ax.set_axisbelow(True)
        if ax is axes[0]:
            ax.set_ylabel("top-1 recall of causal gene")
            ax.legend(frameon=False, loc="upper right", fontsize=8.5)
        style(ax, f"{name}  (n={len(ids)})")
    fig.suptitle("Where the lookup effect actually lives", x=0.008, ha="left",
                 color=INK, fontsize=12, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(f"{FIG}/fig8_interaction.png", dpi=170); plt.close(fig)

# ---------------------------------------------------------------- Arm B
calls_path = f"{REPO}/results/arm_b_calls.tsv"
if os.path.exists(calls_path):
    B = list(csv.DictReader(open(calls_path), delimiter="\t"))
    TIERS = ["P", "LP", "VUS", "LB", "B"]
    keys = [("tier_lookup", "ClinVar lookup"), ("tier_full", "acmgx full"),
            ("tier_svi", "acmgx SVI"), ("tier_blind", "acmgx blinded")]
    fig, axes = plt.subplots(1, 4, figsize=(13.4, 3.7))
    for ax, (k, lab) in zip(axes, keys):
        M = np.zeros((5, 5))
        for r in B:
            if r["tier_vcep"] in TIERS and r[k] in TIERS:
                M[TIERS.index(r["tier_vcep"]), TIERS.index(r[k])] += 1
        Mn = M / np.maximum(M.sum(1, keepdims=True), 1)
        ax.imshow(Mn, cmap=SEQ, vmin=0, vmax=1)
        for i in range(5):
            for j in range(5):
                if M[i, j]:
                    ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=7.5,
                            color="#ffffff" if Mn[i, j] > 0.55 else INK)
        acc = np.trace(M) / max(M.sum(), 1)
        ax.set_xticks(range(5)); ax.set_xticklabels(TIERS, fontsize=8)
        ax.set_yticks(range(5)); ax.set_yticklabels(TIERS, fontsize=8)
        ax.set_xlabel("predicted"); 
        if ax is axes[0]: ax.set_ylabel("ClinGen VCEP")
        style(ax, f"{lab}", f"concordance {acc:.1%}")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig5_confusion.png", dpi=170); plt.close(fig)

crit_path = f"{REPO}/results/criterion_metrics.tsv"
if os.path.exists(crit_path):
    C = [r for r in csv.DictReader(open(crit_path), delimiter="\t") if int(r["n_vcep"]) > 0]
    C.sort(key=lambda r: -int(r["n_vcep"]))
    fig, ax = plt.subplots(figsize=(7.4, 0.34 * len(C) + 1.6))
    y = np.arange(len(C))
    ax.barh(y - 0.19, [float(r["precision"]) for r in C], 0.34, color=S1,
            edgecolor=SURFACE, linewidth=2, label="precision")
    ax.barh(y + 0.19, [float(r["recall"]) for r in C], 0.34, color=S2,
            edgecolor=SURFACE, linewidth=2, label="recall")
    for i, r in enumerate(C):
        ax.text(float(r["precision"]) + 0.012, i - 0.19, f'{float(r["precision"]):.2f}',
                va="center", fontsize=7.5, color=INK)
        ax.text(float(r["recall"]) + 0.012, i + 0.19, f'{float(r["recall"]):.2f}',
                va="center", fontsize=7.5, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{r["criterion"]}  (n={r["n_vcep"]})' for r in C], fontsize=8)
    ax.invert_yaxis(); ax.set_xlim(0, 1.16); ax.xaxis.grid(True); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    style(ax, "Per-criterion agreement with expert panels",
          "acmgx (SVI condition) against the VCEPs' own Met evidence codes")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig6_criterion.png", dpi=170); plt.close(fig)

ceil_path = f"{REPO}/results/ceiling.tsv"
if os.path.exists(ceil_path):
    Z = list(csv.DictReader(open(ceil_path), delimiter="\t"))
    TIERS = ["P", "LP", "VUS", "LB", "B"]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xs = np.arange(len(TIERS))
    keep, ns = [], []
    for t in TIERS:
        g = [r for r in Z if r["tier_vcep"] == t]
        ns.append(len(g))
        keep.append(sum(1 for r in g if r["tier_ceiling"] == t) / max(len(g), 1))
    bars = ax.bar(xs, keep, 0.55, color=S1, edgecolor=SURFACE, linewidth=2)
    for xi, v, n in zip(xs, keep, ns):
        ax.text(xi, v + 0.02, f"{v:.0%}", ha="center", color=INK, fontsize=9.5)
        ax.text(xi, -0.06, f"n={n}", ha="center", color=MUTED, fontsize=8)
    overall = sum(1 for r in Z if r["tier_vcep"] == r["tier_ceiling"]) / len(Z)
    ax.axhline(overall, color=S2, linewidth=2, linestyle="--")
    # park the reference label over the low bars, where there is whitespace
    ax.text(0.5, overall + 0.03, f"overall ceiling {overall:.1%}", color=S2,
            fontsize=8.5, ha="center",
            bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
    ax.set_xticks(xs); ax.set_xticklabels(TIERS)
    ax.set_ylim(-0.1, 1.08); ax.set_ylabel("share of calls a machine could still reach")
    ax.yaxis.grid(True); ax.set_axisbelow(True)
    style(ax, "The ceiling on automated ACMG classification",
          "VCEP calls re-derived from their own evidence, minus the criteria no machine can assess")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig7_ceiling.png", dpi=170); plt.close(fig)

print("figures written to", FIG)
for f in sorted(os.listdir(FIG)):
    print("  ", f)
