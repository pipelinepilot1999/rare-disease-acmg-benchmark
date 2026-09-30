# rare-disease-acmg-benchmark

**How much of automated rare-disease variant interpretation is knowledge-base recall rather than
inference?**

Automated tools for rare disease report strong retrospective performance. Those benchmarks are
built from curated databases, and the tools read the same databases. This project measures the size
of that effect in the two places it matters — phenotype-driven prioritization and ACMG/AMP
classification — using one experimental design in both.

**Every headline number here is a delta between an informed and a blinded condition, never an
absolute.**

| Arm | Instrument | Informed | Blinded | The delta measures |
|---|---|---|---|---|
| A | Exomiser 15.1.0 | ClinVar whitelist on | whitelist off | how much ranking is lookup |
| A′ | Exomiser 15.1.0 | HPO terms supplied | variant-only | how much ranking is phenotype |
| B | `acmgx` (this repo) | PP5/BP6/PS1/PM5 on | all ClinVar-derived criteria off | how much classification is lookup |
| B′ | the VCEPs' own criteria | all evidence | machine-assessable evidence only | **the ceiling** for any automated classifier |

Arm B′ is the one that bounds the problem rather than a tool: it takes each ClinGen expert panel's
own evidence codes, deletes the ones no automated system can assess — functional assays,
segregation, case-control, de novo confirmation, phenotype specificity — re-runs the combining
rules, and reports the tier that survives.

### What it found

- A "classifier" that returns ClinVar's call and reasons about nothing agrees with ClinGen expert
  panels on **97.5%** of 11,499 variants. An engine run under current ClinGen guidance, which may
  not read ClinVar, reaches **36.4%**.
- The **ceiling** — the panels' own evidence minus the criteria no machine can assess — is
  **49.2%** over all 12,967 expert interpretations (**51.1%** over the 11,499 that Arm B
  classifies), and only **9.3% of Pathogenic** calls survive it. 62.8% of expert calls depend on
  evidence no pipeline can generate.
- In prioritization Exomiser's ClinVar whitelist is worth **+0.400 top-1** on the 29.0% of causal
  alleles it actually covers, and exactly zero elsewhere. The marginal average across all cases is
  +0.007, which would have been a wrong conclusion from a correct calculation.
- The largest single effect measured was not a knowledge-base effect: coding a transmitting parent
  unaffected makes the trio inheritance filter discard **96%** of inherited-dominant causal alleles,
  a **70-point** top-1 swing.

See [`docs/results.md`](docs/results.md) for the findings and
[`docs/SPEC.md`](docs/SPEC.md) for the design and how it differs from the v1 specification.

---

## Scope — read this before quoting a number

- **Exome-equivalent regions, trio, GRCh38, germline SNV and small indel.** No structural variants,
  no repeat expansions, no mitochondrial, no non-coding beyond splice region.
- **Cases are synthetic at the VCF level.** Real published patients' HPO profiles and real causal
  alleles from [phenopacket-store](https://github.com/monarch-initiative/phenopacket-store) 0.1.27,
  injected into the GIAB Ashkenazi trio (HG002/3/4). **Nothing here tests variant calling** —
  sensitivity is inherited from the GIAB benchmark, not measured.
- **Phenotypes are never derived from gene-to-phenotype annotations.** Exomiser's hiPhive scores
  against those same annotations; deriving case phenotypes from them would guarantee the match.
  This is the single most important design choice in Arm A.
- **X-linked recessive is not tested.** Only HG002 has a published GIAB chrXY benchmark; the parents
  do not, so no trio can be built on chrX.
- **No real patient data**, no clinical or diagnostic activity, no accredited-laboratory context,
  no variant sign-out, no functional validation. A methods and evaluation project on public data.

Full list: [`docs/limitations.md`](docs/limitations.md).

## Which ACMG criteria are assessed

`acmgx` implements PVS1 (as the ClinGen decision tree over MANE Select exon structure),
PM2/BA1/BS1 (gnomAD v4 grpmax filtering AF), PP3/BP4 (REVEL with the calibrated intervals of
Pejaver et al., single predictor, graded strengths), PM4, BP7, and the ClinVar-derived
PS1/PM5/PP5/BP6.

**Declared not implemented**, never silently skipped: PM1, PP2, BP1, BP3 (need curated hotspot,
domain and constraint annotation not assembled here).

**Not assessable by any automated system**, and therefore the subject of the ceiling analysis:
PS2/PM6 (de novo), PS3/BS3 (functional), PS4 (case-control), PM3/BP2 (allelic data), PP1/BS4
(segregation), PP4 (phenotype specificity), BS2, BP5.

PP5 and BP6 are implemented only so they can be removed: ClinGen SVI recommends they not be used
at all. The `full` condition that includes them is not a baseline — it is the leakage measurement.

## Versions

| Resource | Version |
|---|---|
| Exomiser | 15.1.0, data bundle 2512 (hg38 + phenotype) |
| ClinVar | VCF release 2026-08-31; `variant_summary` same release |
| ClinGen Evidence Repository | harvested via API 2026-09-02, 13,243 interpretations |
| phenopacket-store | 0.1.27 |
| gnomAD | v4 joint, via the gnomAD GraphQL API |
| REVEL | v1.3 |
| GENCODE | v50 (MANE Select transcripts) |
| HPO | 2026-06-23 |
| GIAB benchmark | NIST v4.2.1, GRCh38, chr1–22 |
| ClinGen dosage sensitivity | 2026-09-02 |
| Reference | GCA_000001405.15 GRCh38 no-alt analysis set |

## Layout

```
docs/SPEC.md          the design, and what changed from the v1 specification
docs/decisions.md     every decision, with alternatives and effect on results
docs/results.md       findings
docs/limitations.md   what these numbers do not support
docs/interview-guide.pdf  21-page study guide, regenerated by `make guide`
src/acmgx/            the ACMG engine: criteria, combining rules, assessability, statistics
scripts/              numbered pipeline steps, 01 -> 19 (07 is sharded: 07a/07b/07c)
data/                 case definitions, Arm B variant table, trio carrier alleles
results/              metrics tables and figures
config/               Exomiser properties per condition, trio pedigree
```

## Building

```bash
python3 -m venv .venv && .venv/bin/pip install numpy matplotlib reportlab
make test          # engine unit tests + results-integrity checks
make analyse       # Arm A statistics, Arm B stratifications
make figures       # results/figures/*.png
make guide         # docs/interview-guide.pdf
```

Those four targets run on a fresh clone with no network and no external resources. The numbered
pipeline steps below do not.

## Reproducing

Steps run in order; each is idempotent and logs to `~/rdacmg`.

```bash
export RDACMG_HOME=/path/to/large/scratch     # default ~/rdacmg
bash scripts/render_config.sh                 # Exomiser properties for this machine
bash scripts/01_reference_and_regions.sh      # reference index, exome regions, testable intersection
bash scripts/02_trio_background.sh            # joint trio VCF + Mendelian assertion
python scripts/03_ceiling.py                  # Arm B' ceiling, from ERepo alone
python scripts/04_build_cases.py              # cases from phenopackets
python scripts/05_spike.py                    # per-case VCFs + matched knockouts
python scripts/06_armb_variants.py            # ERepo -> GRCh38 coordinates
python scripts/07b_gnomad_shard.py I N        # gnomAD FAF (sharded); then 07c to merge
python scripts/08_build_mane.py               # MANE Select exon structures for PVS1
python scripts/09_protein_index.py            # PS1/PM5 index from ClinVar
python scripts/11_build_batches.py            # Exomiser job matrix
bash   scripts/12_run_arma.sh                 # Arm A
python scripts/13_extract_ranks.py            # ranks, with an explicit tie policy
python scripts/14_analyse_arma.py             # Arm A statistics
python scripts/16_trio_carriers.py            # P/LP alleles already in the trio background
python scripts/18_ad_penetrant.py 6           # penetrance contrast job matrix
bash   scripts/19_run_penetrant.sh            # ...and run it; then re-run 13 and 14
python scripts/10_run_armb.py                 # Arm B, three conditions + lookup baseline
python scripts/17_armb_strata.py              # Arm B by VCEP and by consequence
python scripts/15_figures.py                  # figures
python tests/test_acmgx.py                    # engine unit tests
python tests/test_results_integrity.py        # published numbers match results/
make guide                                    # rebuild docs/interview-guide.pdf
```

### Paths

Code, `data/` and `results/` are resolved relative to the repository, so a clone runs where it is
checked out. Everything large - downloaded resources, Exomiser working directories, the per-case
VCFs - lives outside the repository under `~/rdacmg`, overridable with `RDACMG_HOME`. Steps 14, 15
and 17 need nothing but the repository and will run on a fresh clone; every other step needs the
external resources.
