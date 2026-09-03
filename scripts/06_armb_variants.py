#!/usr/bin/env python3
"""Resolve the ERepo ground truth to GRCh38 coordinates and basic annotation.

ERepo records carry a ClinVar VariationID, which is the ID column of the ClinVar
VCF, so the join gives coordinates, molecular consequence, gene and ClinVar's own
assertion in one pass. ClinVar's assertion is loaded deliberately: it is both the
evidence for PP5/BP6 and the "lookup baseline" classifier of section 6.3.
"""
import gzip, json, os, sys, csv
from collections import Counter, defaultdict

sys.path.insert(0, os.path.expanduser("~/rdacmg/repo/src"))
from acmgx.combine import parse_code

RES = os.path.expanduser("~/rdacmg/resources")
REPO = os.path.expanduser("~/rdacmg/repo")

TIER_MAP = {
    "Pathogenic": "P", "Likely Pathogenic": "LP", "Uncertain Significance": "VUS",
    "Uncertain Significance - Insufficient Evidence": "VUS",
    "Uncertain Significance - Conflicting Evidence": "VUS",
    "Likely Benign": "LB", "Benign": "B", "Benign - Stand Alone": "B",
}

# ClinVar aggregate classification -> five tiers (the lookup baseline)
def clinvar_tier(clnsig):
    if not clnsig:
        return None
    s = clnsig.replace("_", " ").lower()
    if "conflicting" in s:
        return "VUS"
    for needle, tier in (("likely pathogenic", "LP"), ("likely benign", "LB"),
                         ("pathogenic", "P"), ("benign", "B"),
                         ("uncertain significance", "VUS")):
        if needle in s:
            return tier
    return None

# ------------------------------------------------------- ERepo truth
records = json.load(open(f"{RES}/erepo/erepo_raw.json"))
truth = {}
for r in records:
    vid = r.get("variationId")
    if not vid:
        continue
    for g in r.get("guidelines") or []:
        tier = TIER_MAP.get((g.get("outcome") or {}).get("label"))
        if tier is None:
            continue
        codes, vceps = [], []
        for a in g.get("agents") or []:
            if a.get("affiliation") or a.get("label"):
                vceps.append(a.get("affiliation") or a.get("label"))
            for c in a.get("evidenceCodes") or []:
                if c.get("status") == "Met" and c.get("label"):
                    codes.append(c["label"])
        good = []
        for c in dict.fromkeys(codes):
            try:
                parse_code(c); good.append(c)
            except ValueError:
                pass
        if not good:
            continue
        truth[str(vid)] = {
            "variation_id": str(vid), "caid": r.get("caid"),
            "gene": (r.get("gene") or {}).get("label"),
            "vcep": vceps[0] if vceps else None,
            "condition": (r.get("condition") or {}).get("label"),
            "tier_vcep": tier, "codes": good,
        }
print(f"ERepo interpretations with a ClinVar VariationID: {len(truth)}")

# ------------------------------------------------------- ClinVar join
def info_dict(s):
    d = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d

found = 0
stats = Counter()
with gzip.open(f"{RES}/clinvar/clinvar.vcf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        vid = f[2]
        t = truth.get(vid)
        if t is None:
            continue
        info = info_dict(f[7])
        if f[4] in (".", "") or f[3] in (".", ""):
            stats["no_alt_allele"] += 1
            continue
        t.update({
            "chrom": "chr" + f[0], "pos": int(f[1]), "ref": f[3], "alt": f[4],
            "mc": info.get("MC", ""), "geneinfo": info.get("GENEINFO", ""),
            "clnsig": info.get("CLNSIG", ""), "clnrevstat": info.get("CLNREVSTAT", ""),
            "clnhgvs": info.get("CLNHGVS", ""), "clnvc": info.get("CLNVC", ""),
        })
        found += 1
print(f"joined to ClinVar VCF: {found}  (unmatched: {len(truth)-found})")
print("join issues:", dict(stats))

rows = [t for t in truth.values() if "chrom" in t]

# molecular consequence -> a compact term for the criteria engine
def main_consequence(mc):
    # MC looks like "SO:0001583|missense_variant,SO:0001587|nonsense"
    terms = [p.split("|")[-1] for p in mc.split(",") if "|" in p]
    order = ["nonsense", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
             "initiator_codon_variant", "stop_lost", "missense_variant", "inframe_deletion",
             "inframe_insertion", "splice_site_variant", "synonymous_variant",
             "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
             "non-coding_transcript_variant", "genic_upstream_transcript_variant"]
    for o in order:
        if o in terms:
            return o
    return terms[0] if terms else "unknown"

for r in rows:
    r["consequence"] = main_consequence(r["mc"])
    r["clinvar_tier"] = clinvar_tier(r["clnsig"])
    r["is_expert_panel"] = "reviewed_by_expert_panel" in (r["clnrevstat"] or "")

print("\nconsequence distribution:", Counter(r["consequence"] for r in rows).most_common(12))
print("\nexpert tier distribution:", Counter(r["tier_vcep"] for r in rows).most_common())
print("\nClinVar review status  :", Counter(r["clnrevstat"] for r in rows).most_common(6))

os.makedirs(f"{REPO}/data", exist_ok=True)
cols = ["variation_id","caid","chrom","pos","ref","alt","gene","vcep","condition",
        "tier_vcep","consequence","mc","clnsig","clnrevstat","clinvar_tier",
        "is_expert_panel","codes"]
with open(f"{REPO}/data/armb_variants.tsv","w",newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in sorted(rows, key=lambda x: (x["chrom"], x["pos"])):
        w.writerow([r.get(c) if c != "codes" else ";".join(r["codes"]) for c in cols])
print(f"\nwrote {len(rows)} rows -> {REPO}/data/armb_variants.tsv")
