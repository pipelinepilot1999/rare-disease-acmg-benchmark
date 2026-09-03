#!/usr/bin/env python3
"""Build Arm A cases from phenopacket-store.

Each case is a real published patient's HPO profile plus that patient's causal
allele, spiked into the GIAB Ashkenazi trio under the inheritance model implied
by the reported zygosity. Phenotypes are never derived from gene-to-phenotype
annotations - that is the circularity this project exists to avoid.
"""
import json, glob, os, random, sys, csv
from collections import Counter, defaultdict

RES = os.path.expanduser("~/rdacmg/resources")
WORK = os.path.expanduser("~/rdacmg/work")
REPO = os.path.expanduser("~/rdacmg/repo")
random.seed(20260902)

# ---------------------------------------------------------------- regions
def load_bed(path):
    d = defaultdict(list)
    for line in open(path):
        f = line.split()
        d[f[0]].append((int(f[1]), int(f[2])))
    return {k: sorted(v) for k, v in d.items()}

import bisect
class Regions:
    def __init__(self, path):
        self.d = load_bed(path)
        self.starts = {c: [s for s, _ in v] for c, v in self.d.items()}
    def contains(self, chrom, pos):        # pos is 1-based
        v = self.d.get(chrom)
        if not v: return False
        i = bisect.bisect_right(self.starts[chrom], pos - 1) - 1
        return i >= 0 and v[i][0] <= pos - 1 < v[i][1]

regions = Regions(f"{WORK}/trio_hiconf.bed")
print(f"testable regions loaded: {sum(len(v) for v in regions.d.values())} intervals")

# ------------------------------------------------- background trio variants
background = set()
import gzip, subprocess
p = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", f"{WORK}/trio.vcf.gz"],
                   capture_output=True, text=True, check=True)
for line in p.stdout.splitlines():
    c, pos, ref, alt = line.split("\t")
    background.add((c, int(pos), ref, alt))
print(f"background trio variants: {len(background)}")

# ---------------------------------------------------------------- HPO
hp = json.load(open(f"{RES}/hpo/hp.json"))["graphs"][0]
parents = defaultdict(set)
labels = {}
for n in hp["nodes"]:
    if n.get("id", "").startswith("http://purl.obolibrary.org/obo/HP_"):
        labels[n["id"].rsplit("/", 1)[-1].replace("_", ":")] = (n.get("lbl") or "")
for e in hp["edges"]:
    if e.get("pred") == "is_a":
        s = e["sub"].rsplit("/", 1)[-1].replace("_", ":")
        o = e["obj"].rsplit("/", 1)[-1].replace("_", ":")
        if s.startswith("HP:") and o.startswith("HP:"):
            parents[s].add(o)
ALL_HPO = [t for t in labels if t.startswith("HP:")]
PHENO_ROOT = "HP:0000118"

def generalize(term, levels=1):
    """Walk up the ontology, stopping before the phenotypic-abnormality root."""
    cur = term
    for _ in range(levels):
        ps = [p for p in parents.get(cur, ()) if p != PHENO_ROOT]
        if not ps:
            break
        cur = sorted(ps)[0]
    return cur

# ---------------------------------------------------------------- packets
files = sorted(glob.glob(f"{RES}/phenopackets/store/**/*.json", recursive=True))
stats = Counter()
candidates = []

for fp in files:
    d = json.load(open(fp))
    feats = [f["type"]["id"] for f in d.get("phenotypicFeatures", [])
             if not f.get("excluded") and f.get("type", {}).get("id", "").startswith("HP:")]
    feats = [t for t in dict.fromkeys(feats) if t in labels]
    if len(feats) < 3:
        stats["too_few_hpo"] += 1; continue

    variants, disease, gene = [], None, None
    for it in d.get("interpretations") or []:
        diag = it.get("diagnosis") or {}
        disease = disease or diag.get("disease")
        for g in diag.get("genomicInterpretations") or []:
            vd = (g.get("variantInterpretation") or {}).get("variationDescriptor") or {}
            rec = vd.get("vcfRecord")
            if not rec or rec.get("genomeAssembly") not in ("hg38", "GRCh38"):
                continue
            gene = gene or (vd.get("geneContext") or {}).get("symbol")
            variants.append({
                "chrom": rec["chrom"] if str(rec["chrom"]).startswith("chr") else "chr" + str(rec["chrom"]),
                "pos": int(rec["pos"]), "ref": rec["ref"], "alt": rec["alt"],
                "zyg": ((vd.get("allelicState") or {}).get("label") or "").lower(),
                "gene": (vd.get("geneContext") or {}).get("symbol"),
            })
    if not variants:
        stats["no_variant"] += 1; continue
    if any(len(v["ref"]) > 50 or len(v["alt"]) > 50 for v in variants):
        stats["large_indel"] += 1; continue
    if not all(regions.contains(v["chrom"], v["pos"]) for v in variants):
        stats["outside_testable_region"] += 1; continue
    if any((v["chrom"], v["pos"], v["ref"], v["alt"]) in background for v in variants):
        stats["already_in_trio"] += 1; continue

    # ------------------------------------------------ inheritance model
    zygs = [v["zyg"] for v in variants]
    genes = {v["gene"] for v in variants}
    if len(variants) == 2 and len(genes) == 1 and all(z == "heterozygous" for z in zygs):
        model = "comphet"
    elif len(variants) == 1 and zygs[0] == "homozygous":
        model = "ar_hom"
    elif len(variants) == 1 and zygs[0] == "hemizygous" and variants[0]["chrom"] == "chrX":
        model = "xlr"
    elif len(variants) == 1 and zygs[0] == "heterozygous":
        model = "denovo"        # a subset is reassigned to inherited dominant below
    else:
        stats["unmodellable_zygosity"] += 1; continue

    candidates.append({
        "case_id": os.path.splitext(os.path.basename(fp))[0],
        "source": os.path.relpath(fp, f"{RES}/phenopackets/store"),
        "gene": gene, "disease_id": (disease or {}).get("id"),
        "disease": (disease or {}).get("label"),
        "model": model, "variants": variants, "hpo": feats,
    })
    stats["usable"] += 1

print("\ncase filtering:")
for k, v in stats.most_common():
    print(f"  {k:28s} {v:6d}")

# Reassign a third of the de novo pool to inherited-dominant, so the cost of the
# de novo inheritance filter is measurable rather than assumed.
dn = [c for c in candidates if c["model"] == "denovo"]
random.shuffle(dn)
for c in dn[: len(dn) // 3]:
    c["model"] = "ad_inherited"

# ------------------------------------------------- balance and sub-sample
by_model = defaultdict(list)
for c in candidates:
    by_model[c["model"]].append(c)
print("\nby inheritance model (all usable):")
for m, v in sorted(by_model.items()):
    print(f"  {m:14s} {len(v):5d}  ({len({c['gene'] for c in v})} genes)")

# One case per (gene, model) keeps the set from being dominated by heavily
# published genes, then cap per model.
CAP = 400
selected = []
for m, v in sorted(by_model.items()):
    seen, keep = set(), []
    random.shuffle(v)
    for c in v:
        if c["gene"] in seen:
            continue
        seen.add(c["gene"]); keep.append(c)
    selected.extend(keep[:CAP])
print(f"\nselected {len(selected)} cases "
      f"({len({c['gene'] for c in selected})} distinct genes, "
      f"{len({c['disease_id'] for c in selected})} diseases)")
for m, v in sorted(Counter(c['model'] for c in selected).items()):
    print(f"  {m:14s} {v:5d}")

# ------------------------------------------------- HPO perturbation tiers
pool = sorted(ALL_HPO)
for c in selected:
    terms = c["hpo"]
    c["hpo_sets"] = {
        "as_reported": terms,
        "imprecise": sorted({generalize(t, 2) for t in terms}),
        "noisy": terms + random.sample(pool, max(3, len(terms) // 2)),
        "truncated": terms[: max(2, len(terms) // 3)],
    }

os.makedirs(f"{REPO}/data", exist_ok=True)
with open(f"{REPO}/data/cases.json", "w") as fh:
    json.dump(selected, fh, indent=1)

with open(f"{REPO}/data/cases.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["case_id", "gene", "disease_id", "disease", "model", "n_variants",
                "n_hpo", "variants", "source"])
    for c in selected:
        w.writerow([c["case_id"], c["gene"], c["disease_id"], c["disease"], c["model"],
                    len(c["variants"]), len(c["hpo"]),
                    ";".join(f"{v['chrom']}:{v['pos']}{v['ref']}>{v['alt']}" for v in c["variants"]),
                    c["source"]])
print(f"\nwrote {REPO}/data/cases.json and cases.tsv")
