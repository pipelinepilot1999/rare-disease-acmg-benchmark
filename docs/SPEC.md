# rare-disease-acmg-benchmark — build specification (v2)

**Supersedes:** `rare-disease-acmg-benchmark-spec.md` (v1).
**Execution mode:** Autonomous. Decide, proceed, log.

---

## 0. What changed from v1, and why

v1 proposed two loosely related arms — phenotype-driven prioritization, and rule-based ACMG
classification — and correctly identified that benchmarking an ACMG classifier against ClinVar
while letting it read ClinVar measures nothing. It then applied that insight to one arm only.

Five defects made v1 unbuildable as written:

1. **Arm A had the same circularity as Arm B, uncorrected.** Exomiser's ClinVar whitelist gives
   P/LP variants a maximal variant score that bypasses all filtering. v1 spiked ClinVar P/LP
   variants. Recall would have measured ClinVar lookup, not prioritization. Compounding it, v1
   derived case HPO terms from gene→phenotype annotations and then ranked with hiPhive, which
   scores against those same annotations.
2. **Arm B's headline analysis had no ground truth.** v1 called criterion-level disagreement "the
   analytical core" but sourced ground truth from ClinVar, which records tiers, not criteria.
3. **Negative controls were n=1.** One unmodified trio replicated 20–30 times with different
   phenotypes is one genome, not 20–30 negatives.
4. **The rule set was 2015-vintage.** PP5/BP6 are deprecated by ClinGen SVI, not merely circular;
   PVS1 is a four-strength decision tree; PP3/BP4 require a single calibrated predictor; PM2 is
   Supporting; combining should also be reported in points.
5. **The two arms answered different questions,** so the project had no thesis.

v2 fixes all five, and the fix to (5) subsumes the rest.

---

## 1. The question

> **How much of automated rare-disease variant interpretation is knowledge-base recall rather
> than inference?**

Automated tools report strong performance on retrospective benchmarks. Those benchmarks are built
from curated databases, and the tools read those same databases. This project measures the size of
that effect in the two places it matters, using one experimental design in both:

**Every headline number is a delta between an informed and a blinded condition, never an absolute.**

| | Instrument | Informed | Blinded | The delta measures |
|---|---|---|---|---|
| **Arm A** | Exomiser 15.1.0 | ClinVar whitelist + ClinVar pathogenicity on | both off | how much ranking is lookup |
| **Arm A′** | Exomiser 15.1.0 | HPO terms supplied | variant-only | how much ranking is phenotype |
| **Arm B** | `acmgx` (this repo) | PP5/BP6/PS1/PM5 on | all ClinVar-derived criteria off | how much classification is lookup |
| **Arm B′** | VCEP criteria themselves | all evidence | machine-assessable evidence only | the **ceiling** for any automated classifier |

Arm B′ is the result v1 could not have produced and is the most interesting one: it bounds what
*any* rule engine can achieve, independent of implementation quality, by re-running the combining
rules on the expert panels' own evidence with the un-automatable criteria removed.

---

## 2. Environment (corrected)

v1 specified a 200 GB volume and spent most of its section 2 rationing against it. That
constraint drove real design decisions - cutting predictors, skipping splice prediction - and
none of them were necessary.

Requirements, as actually used for the published run:

- Linux x86-64, **8 vCPU, 62 GB RAM**, CPU-only. Arm A took ~3.5 h at that size; it scales
  roughly linearly with core count.
- **~130 GB of free disk** outside the repository: ~92 GB of downloaded resources (the Exomiser
  2512 bundles dominate at 37 GB compressed), ~7 GB of per-case VCFs, ~1 GB of Exomiser output.
  Set `RDACMG_HOME` to point at it; the default is `~/rdacmg`.
- **OpenJDK 21** for Exomiser, Python 3.11+, and `bcftools`/`samtools` 1.21.

**Consequence:** v1's entire disk-budget section is void. The instructions to cut REMM, to cut
full CADD, to prefer targeted plugin datasets over dbNSFP, and the 160 GB re-plan trigger were
all responses to a disk constraint that a modern volume does not have. Take the full resource set.

Long jobs run under `tmux`. Free space is checked before each acquisition step and logged.

---

## 3. Scope

Exome-equivalent regions, trio, GRCh38, germline SNV/small indel. No SV, no repeat expansions,
no mitochondrial, no non-coding beyond splice region. Stated plainly in the README.

**Exome regions are defined as GENCODE v50 basic CDS ± 20 bp**, not a vendor capture kit. Vendor
BEDs are license-encumbered and kit-specific; a GENCODE-derived region set is reproducible by
anyone and is the honest description of what is being tested.

---

## 4. Data

### 4.1 Genomic background

GIAB Ashkenazi trio (HG002 son, HG003 father, HG004 mother), NIST v4.2.1 GRCh38 benchmark VCFs,
restricted to the intersection of the three high-confidence BEDs and the exome region set.

Two properties of this background must be measured and reported, not assumed:

- **HG002/3/4 carry real ClinVar P/LP alleles.** A "negative" case containing a genuine pathogenic
  carrier allele is not negative. Enumerate every P/LP allele in the trio before building cases and
  publish the list in `docs/limitations.md`.
- **The v4.2.1 benchmark excludes difficult regions.** Causal variants cannot be placed there, so
  recall is systematically optimistic with respect to a real exome. Quantify the excluded fraction
  of the exome region set and report it.

### 4.2 Cases — real phenotypes, not derived ones

Cases come from **phenopacket-store 0.1.27** (GA4GH Phenopackets curated from published case
reports). Each supplies a real patient's HPO profile and that patient's causal allele.

For each usable phenopacket: map the causal allele to GRCh38, and spike it into the trio under the
inheritance model implied by the phenopacket's genotype and the disease MOI (de novo, AR
homozygous, compound heterozygous, X-linked hemizygous).

This replaces v1 §4.2/§4.4. It removes the HPO circularity outright, supplies real phenotypic
noise and incompleteness rather than simulated noise, and raises usable *n* by an order of
magnitude. Cases are stored **as phenopackets**, which Exomiser consumes natively, not as a TSV.

Perturbation tiers are applied to the real terms, not used to manufacture them:
`as_reported` (unmodified) · `imprecise` (each term generalized up the ontology) ·
`noisy` (unrelated terms added) · `truncated` (term set reduced).

### 4.3 Negative controls — matched, not replicated

For every positive case, a **matched knockout**: identical background, identical HPO terms, causal
variant removed. This is a paired design. It fixes v1's n=1 defect, and it licenses paired
statistics (McNemar, paired bootstrap) that are far better powered than the unpaired comparison
v1 implied.

Additionally, **phenotype-swapped** cases: causal variant present, HPO terms from an unrelated
disease. These separate "the variant looks bad" from "the variant matches the patient."

### 4.4 Arm B ground truth

**ClinGen Evidence Repository, all 13,243 published interpretations**, harvested via its API.
Each record carries the VCEP's final tier, the individual evidence codes with strength modifiers
(e.g. `PM3_Strong`, `PP4_Moderate`), each code's Met/Not-Met status, the VCEP affiliation, and a
ClinGen Allele Registry CAID.

This is the ground truth v1 needed and did not have. ClinVar 3-star records are downstream
summaries of these; using ERepo directly gives criterion-level truth and the VCEP identity needed
for stratification.

Held out entirely from Arm A.

---

## 5. Arm A — prioritization

### 5.1 Conditions

Full factorial over two blinding axes, plus baselines:

| Condition | ClinVar | Phenotype |
|---|---|---|
| `informed` | on | on |
| `clinvar_blind` | **off** | on |
| `pheno_blind` | on | **off** |
| `double_blind` | **off** | **off** |

Baselines: frequency+consequence filter only (rank by REVEL); random ranking within the surviving
set. Baselines bound the bottom of the scale, exactly as the blinded conditions bound the top.

ClinVar blinding means: Exomiser's ClinVar whitelist disabled **and** ClinVar removed as a
pathogenicity source. Verified by assertion, not by configuration intent (§10).

### 5.2 Filtering cascade

Applied in order with survivor counts recorded per step (the attrition table is a deliverable):

1. Quality — GIAB benchmark calls are already high-confidence; the spiked variants' FORMAT fields
   are sampled from matched real variants so they are not trivially separable (§10).
2. Frequency — **gnomAD v4 grpmax filtering AF** (95% CI lower bound), not raw popmax. Thresholds
   per inheritance model, justified from prevalence, logged as decisions.
3. Consequence — protein-altering and splice-region on **MANE Select** transcripts.
4. Inheritance — de novo / AR-hom / comp-het / XLR as separate tracks. Comp-het phase is inferred
   from parental genotypes, not read-backed; stated as such.

**The Ashkenazi interaction is a named sub-analysis, not a footnote.** gnomAD's grpmax FAF
deliberately excludes Ashkenazi Jewish, Finnish, Amish and Middle Eastern groups because
pathogenic alleles drift high in bottlenecked populations. The proband background *is* Ashkenazi.
Report causal-variant survival at step 2 separately for AJ founder alleles.

### 5.3 Metrics

Rank of the causal **gene** is primary; the causal **variant** rank is reported beside it.
Exomiser scores per gene and every variant inherits that score, so a variant-level rank counts a
gene's other alleles against the causal one - for compound heterozygotes it makes rank 1
unreachable outright.

**Tie policy is explicit and pessimistic**: where the causal variant ties with others, its rank is
the worst position in the tie block. Optimistic and mean-rank variants are reported alongside so
the effect of the policy is visible. Undeclared tie-breaking is a common source of inflated top-1.

Top-1/5/10 recall, MRR, full rank distribution. Bootstrap 95% CIs on everything. Paired tests
across conditions and across HPO tiers. Minimum detectable effect stated before the run.

Breakouts by inheritance model and HPO tier are reported **with CIs**, and where a cell is too
small to support a claim, that is said rather than the point estimate being discussed.

---

## 6. Arm B — classification

### 6.1 `acmgx`

A transparent ACMG/AMP engine written for this project, because the analysis requires
criterion-level toggling that no off-the-shelf tool exposes. Current SVI guidance, not 2015:

| Criterion | Implementation |
|---|---|
| PVS1 | ClinGen decision tree; emits PVS1 / _Strong / _Moderate / _Supporting |
| PM2 | **Supporting** strength; gnomAD v4 grpmax FAF |
| BA1, BS1 | grpmax FAF thresholds |
| PP3, BP4 | **single** calibrated predictor (REVEL) with graded strengths per Pejaver et al.; SpliceAI for splice, calibrated; never stacked |
| PS1, PM5 | amino-acid identity against ClinVar — *ClinVar-derived* |
| PP5, BP6 | ClinVar assertion — *ClinVar-derived, and SVI-deprecated* |
| PM1, PM4, BP3, BP7 | domain/hotspot, in-frame length change, repeat region, synonymous-no-splice |
| Unassessable | PS2/PM6, PS3/BS3, PS4, PM3, PP1/BS4, PP4, BS2, BP2/BP5 — declared, never silently skipped |

Combining reports **both** the 2015 categorical rules and the **Tavtigian points** total.

### 6.2 Conditions

- `full` — everything on, including PP5/BP6. Non-compliant with current guidance; included solely
  to measure the leakage. This is the condition most published "high concordance" numbers resemble.
- `svi` — PP5/BP6 off (current best practice), PS1/PM5 on.
- `blind` — every ClinVar-derived criterion off. Classification from sequence, frequency and
  calibrated predictors alone.

### 6.3 The lookup baseline

A "classifier" that ignores the variant entirely and returns ClinVar's aggregate classification.
It performs no reasoning. Its concordance with ERepo is the number every other result must be read
against, and it is reported first in `docs/results.md`.

### 6.4 The ceiling (Arm B′)

Take each ERepo record's **own** VCEP-assigned criteria, delete the codes no automated system can
assess (functional, segregation, case-control, de novo confirmation, phenotype specificity,
allelic data), re-run the combining rules on what remains, and record the tier.

This yields, per variant, the best tier any perfect rule engine could reach. Aggregated it gives
the **ceiling on automated concordance** — a property of the ACMG framework and the evidence
distribution, not of any tool. Reported as: ceiling, `acmgx` performance, and the gap.

### 6.5 Metrics

Per-tier precision/recall/F1 with class counts always shown; five-tier confusion matrices per
condition; **per-criterion** precision/recall/F1 against VCEP Met/Not-Met; points-scale error
(signed distance, which is more informative than tier agreement); VUS rate vs expert.

**Stratified by VCEP and by gene**, because VCEP calls apply gene-specific criteria
specifications. A generic engine scored against gene-specific expert rules conflates "wrong rules"
with "rules not specialized" — the stratification separates them.

---

## 7. Implementation

Arm A fans out with **Exomiser's own batch mode**, one JVM per worker, sharded across workers by
a shuffled job list. Nextflow was the v1 plan and is the wrong tool here: it would start a JVM per
task, and JVM startup plus loading the 287,098-variant ClinVar whitelist and the MVStores is most
of the per-case cost. Batch mode amortises all of it across a worker's whole shard, and measured
throughput is roughly twice what per-task invocation gives. Resumption is by re-running the
generator, which skips jobs whose output already exists.

Arm B is a Python package. It is embarrassingly parallel over variants and needs no workflow engine.

Versions are pinned in `versions.yml`. Smoke test on one case end-to-end, plus the section 10
blinding assertion, before any full run.

---

## 8. Deliverables

```
README.md                     docs/SPEC.md          docs/decisions.md
docs/results.md               docs/limitations.md   docs/ceiling.md
src/acmgx/                    scripts/              config/
data/cases/*.json             data/erepo_truth.tsv  data/trio_carriers.tsv
results/arm_a_ranks.tsv       results/arm_b_calls.tsv
results/criterion_metrics.tsv results/ceiling.tsv   results/figures/
```

### Figures

1. Causal-variant recall across the four blinding conditions (the headline)
2. Full rank distribution by condition, not just the top-N points
3. Recall by HPO perturbation tier, paired within case
4. Causal-variant score against its own matched knockout's best hit
5. Tier confusion matrices: lookup baseline / `full` / `svi` / `blind`
6. Per-criterion precision and recall against VCEP evidence codes
7. The ceiling, by expert tier
8. The ClinVar x phenotype 2x2, split by whitelist membership

Filter attrition is reported as a table rather than a figure: with a single shared background the
per-stage survivor counts are near-constant across cases, so a plot of them carries one number.

---

## 9. Reporting

`docs/results.md` answers §1 in its first paragraph, with a number and a CI.

Where a blinded number is much worse than an informed one, that is the finding — not a
disappointment to be softened. The likely honest outcome is that both arms lose a large fraction
of their apparent performance under blinding, and that Arm B′'s ceiling sits well below 100%.
Report it that way.

---

## 10. Assertions

**Setup** — reference build and contig naming identical across trio VCF, GENCODE, ClinVar, REVEL,
Exomiser bundle. Trio relationships confirmed from genotypes (Mendelian error rate), not filenames.
All variants left-aligned and decomposed; representation reconciled against ClinVar.

**Cases** — each spiked variant present at intended zygosity; parental genotypes match the model,
de novo checked individually against both parents; comp-het pairs come from different parents;
knockouts contain no spiked variant; spiked FORMAT fields statistically indistinguishable from
matched real variants (two-sample test on DP/GQ, reported).

**Arm A** — ClinVar blinding verified by confirming a known ClinVar P/LP variant scores *lower*
under `clinvar_blind` than under `informed`; a configuration flag alone is not evidence.
Causal-variant survival through the cascade reported per condition; a low value is investigated
before it is reported, but is not tuned away.

**Arm B** — criterion-level output captured for every variant; `blind` verified by confirming
PP5/BP6/PS1/PM5 fire zero times; per-class counts always beside aggregate figures; ceiling
computed from VCEP codes without reference to `acmgx` output.

---

## 11. Decision log

`docs/decisions.md`, sequential entries: **Chose / Alternatives / Reasoning / Source / Effect on
results.** Every threshold, every resource choice, every point where the spec left it open.

---

## 12. Scope boundary

No real patient data beyond published, de-identified case reports already in the public
phenopacket corpus. No clinical or diagnostic activity, no accredited-laboratory context, no
variant sign-out, no functional validation. A methods and evaluation project on public data.
