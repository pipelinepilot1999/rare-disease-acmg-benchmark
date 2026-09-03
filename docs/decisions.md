# Decision log

Format: **Chose / Alternatives / Reasoning / Source / Effect on results.**

---

### D-001: Reframe the project around a single question
**Chose:** "How much of automated rare-disease interpretation is knowledge-base recall rather than
inference?", with both arms reporting deltas between informed and blinded conditions.
**Alternatives:** keep v1's two independent arms; keep absolute concordance as the headline.
**Reasoning:** v1 identified ClinVar circularity in Arm B and corrected for it, then built an Arm A
that had the same defect uncorrected. One question makes the correction structural rather than
optional, and makes the two arms answer the same thing with different instruments.
**Source:** own judgment, prompted by the ClinVar leakage in Exomiser's whitelist (D-004).
**Effect:** every headline number is now a paired delta with a confidence interval.

### D-002: Cases from phenopacket-store, not gene-to-phenotype annotations
**Chose:** phenopacket-store 0.1.27 (10,377 packets), real published patients' HPO profiles paired
with their causal alleles.
**Alternatives:** v1's plan - derive HPO terms from HPO's gene-to-phenotype file.
**Reasoning:** Exomiser's hiPhive scores phenotype-gene similarity against the same gene-to-phenotype
annotations. Deriving the case phenotype from them and then ranking with hiPhive guarantees the
match. Real case-report phenotypes also carry real incompleteness - median 7 terms - which is what
a referral actually looks like.
**Source:** Danis et al., *A corpus of GA4GH Phenopackets* (2024).
**Effect:** 1,034 usable cases across 628 genes and 660 diseases, with no HPO circularity.

### D-003: X-linked recessive dropped
**Chose:** report four inheritance models (de novo, inherited dominant, AR homozygous, compound
heterozygous); no XLR arm.
**Alternatives:** add the GIAB chrXY v1.0 benchmark.
**Reasoning:** only HG002 has a published chrXY benchmark; HG003 and HG004 do not, so no trio can
be constructed on chrX. 351 chrX variants in the corpus (3.2%) are therefore untestable. This is a
data limitation, not a design choice.
**Source:** GIAB FTP listing, checked 2026-09-02.
**Effect:** chrX excluded from the region set; 1.7 Mb (3.7% of the exome) not testable.

### D-004: ClinVar blinding for Arm A is the whitelist switch, not database removal
**Chose:** `exomiser.hg38.use-clinvar-white-list=false` against the normal data directory.
**Alternatives:** a data directory with the ClinVar MVStore removed; an emptied ClinVar database.
**Reasoning:** Exomiser will not start without `2512_hg38_clinvar.mv.db`, and the file is an
MVStore, not an H2 SQL database, so it cannot be emptied with the H2 tooling shipped alongside it.
The whitelist is the channel through which ClinVar reaches the *ranking* score: 287,098 whitelisted
variants receive a maximal variant score and bypass the filter cascade. ClinVar remains available
to the annotation-only `EXOMISER_ACMG_*` columns, which do not feed the score.
**Source:** Exomiser 15.1.0 startup log; Exomiser documentation on the ClinVar whitelist.
**Effect:** the blind removes the score channel, not the annotation channel. Stated in
`docs/limitations.md` as a partial blind.

### D-005: 63.6% of causal variants are ClinVar P/LP
**Chose:** measure and report this before running anything.
**Reasoning:** it bounds how much of Arm A could be lookup. 933 of 1,153 causal alleles (80.9%) are
in ClinVar at all; 733 (63.6%) are P/LP and therefore whitelist candidates.
**Source:** ClinVar VCF release 2026-08-31, joined to the case manifest.
**Effect:** frames the interpretation of the ClinVar blinding delta.

### D-006: Exome regions from GENCODE, not a vendor capture kit
**Chose:** GENCODE v50 basic CDS +/- 20 bp, autosomes.
**Alternatives:** Twist/Agilent/IDT capture BEDs.
**Reasoning:** vendor BEDs are licence-encumbered and kit-specific; a GENCODE-derived set is
reproducible by anyone and describes exactly what was tested.
**Effect:** 45.6 Mb region set; 43.9 Mb autosomal; 41.3 Mb testable after intersecting the three
GIAB high-confidence BEDs. 6.0% of the autosomal exome is lost to GIAB's difficult-region
exclusions and is reported as an optimism bias.

### D-007: Missing genotypes inside the high-confidence intersection are hom-ref
**Chose:** `bcftools merge -0`.
**Reasoning:** the GIAB benchmark asserts high confidence over the BED; absence of a call inside it
is a confident reference call, not missingness.
**Effect:** 42,778-variant trio VCF; Mendelian violation rate 0.035% over 42,778 informative sites,
which confirms the trio relationship from genotypes rather than filenames.

### D-008: Spiked FORMAT fields sampled from matched real variants
**Chose:** DP and GQ for each injected genotype drawn from real trio variants matched on variant
type (SNV/indel) and genotype.
**Alternatives:** fixed high-quality values, as VCF-level spiking usually does.
**Reasoning:** constant perfect quality values make a spiked variant trivially separable and mean
the quality filter can never touch a causal variant, inflating recall and making the first row of
the attrition table meaningless.
**Effect:** causal variants are subject to the same quality filtering as real ones.

### D-009: Matched knockouts as negative controls
**Chose:** for every positive case, an otherwise identical case with the causal variant removed.
**Alternatives:** v1's plan - 20-30 unmodified trios with HPO attached.
**Reasoning:** v1's negatives were one genome replicated, so they varied only in phenotype. A
matched knockout holds background and phenotype fixed and varies only the thing under test, which
also licenses paired statistics.
**Effect:** 1,034 matched negatives instead of 20-30 pseudo-replicates.

### D-010: gnomAD v4 joint grpmax filtering AF, not raw popmax
**Chose:** `joint.faf95.popmax` from the gnomAD v4 API for BA1/BS1/PM2.
**Alternatives:** raw population-maximum AF, as v1 specified.
**Reasoning:** ClinGen's March 2024 guidance names the joint exome+genome filtering AF - the 95% CI
lower bound - as the value for BA1 and BS1. It also excludes Ashkenazi Jewish, Finnish, Amish and
Middle Eastern groups, because pathogenic alleles drift high in bottlenecked populations.
**Source:** ClinGen, *Guidance to VCEPs regarding the use of gnomAD v4*, v1.0, March 2024.
**Effect:** the proband background is Ashkenazi, so this exclusion is a named sub-analysis rather
than a footnote. Exomiser's default frequency sources independently omit `GNOMAD_E_ASJ` and
`GNOMAD_G_ASJ`, and that default was kept.

### D-011: PP3/BP4 from one calibrated predictor
**Chose:** REVEL alone, with the score intervals of Pejaver et al. 2022 mapped to graded strengths
(PP3 at 0.644/0.773/0.932 for supporting/moderate/strong; BP4 at 0.290/0.183/0.016/0.003 for
supporting/moderate/strong/very strong).
**Alternatives:** v1's "REVEL, CADD or equivalent"; ANNOVAR/dbNSFP majority voting across predictors.
**Reasoning:** ClinGen SVI recommends a single calibrated tool and explicitly warns against
combining correlated predictors. Majority voting is not compliant with current guidance.
**Source:** Pejaver et al., *Calibration of computational tools...*, AJHG 2022, Table 2 (thresholds
read from the published table, not from tool defaults).
**Effect:** PP3/BP4 fire with calibrated strengths rather than a binary call.

### D-012: PM2 at Supporting
**Chose:** emit `PM2_Supporting`, never `PM2`.
**Source:** ClinGen SVI, *Recommendation for Absence/Rarity Criterion PM2*, v1.0.
**Effect:** matches observed VCEP practice - `PM2_Supporting` is the single most-used evidence code
in the Evidence Repository (7,692 uses), against 1,592 for unmodified `PM2`.

### D-013: PVS1 as a decision tree over MANE Select exon structure
**Chose:** implement NMD escape (last exon, or last 50 nt of the penultimate exon) against MANE
Select transcripts, and modulate strength by ClinGen dosage-sensitivity curation.
**Alternatives:** flat PVS1 for any LoF consequence, as InterVar does.
**Reasoning:** since 2018 PVS1 has four strength levels; a flat implementation both over- and
under-calls. Genes with no ClinGen dosage curation get one strength level less, reflecting the
uncertainty about whether LoF is the disease mechanism.
**Source:** Abou Tayoun et al. 2018; ClinGen dosage sensitivity curation list, 2026-09-02.
**Effect:** MANE Select structures built for 19,229 genes.

### D-014: PP5 and BP6 implemented only so they can be removed
**Chose:** implement them, then disable them in every condition except `full`.
**Reasoning:** ClinGen SVI recommends they not be used at all. v1 framed the PP5/BP6-enabled run as
the baseline and the disabled run as "corrected"; that has it backwards. The compliant run is the
baseline, and the enabled run measures how much apparent performance is circular.
**Source:** Biesecker and Harrison, *Recommendation for reputable source PP5 and BP6*.
**Effect:** three conditions - `full`, `svi`, `blind` - plus a ClinVar lookup baseline.

### D-015: ClinGen Evidence Repository as Arm B ground truth
**Chose:** all 13,243 ERepo interpretations, harvested via its API (`matchLimit`/`matchSkip`).
**Alternatives:** v1's plan - ClinVar records with expert-panel review status.
**Reasoning:** v1 named criterion-level disagreement "the analytical core" but sourced truth from
ClinVar, which records tiers and not criteria. ERepo publishes the VCEP's individual evidence codes
with strength modifiers and Met/Not-Met status, which is what that analysis needs.
**Effect:** 12,967 interpretations usable for the ceiling analysis; 11,499 resolved to GRCh38
coordinates by joining VariationID to the ClinVar VCF.

### D-016: Combining reported under both rule sets
**Chose:** report the 2015 categorical rules and the Tavtigian point total for every variant.
**Reasoning:** VCEP specifications are migrating to points; the two disagree on a small set of
evidence combinations and reporting one hides that.
**Effect:** our 2015 implementation reproduces the VCEPs' own published tier for 92.2% of
interpretations when given all of their evidence, which is what makes the ceiling analysis credible.

### D-017: PM1, PP2, BP1, BP3 declared not implemented
**Chose:** leave them unimplemented and declare them, rather than approximate them.
**Reasoning:** defensible versions need curated hotspot/domain annotation and gene-level missense
constraint that were not assembled. A poor approximation adds noise to the criterion-level
comparison, which is the arm's main analysis.
**Effect:** acmgx systematically under-calls relative to VCEPs on these; quantified in
`results/criterion_metrics.tsv`.

### D-018: SpliceAI run inside Exomiser, not standalone for Arm B
**Chose:** Arm A uses Exomiser's bundled SPLICE_AI pathogenicity source; Arm B's BP7 uses a
distance-to-exon-boundary heuristic instead.
**Alternatives:** run SpliceAI standalone over the 11,499 Arm B variants.
**Reasoning:** the standalone run costs several CPU-hours for a criterion that accounts for 1,107
of ~41,000 VCEP evidence codes. Arm A, where splice prediction actually affects ranking, does get
real SpliceAI.
**Effect:** BP7 is conservative in Arm B - a synonymous variant within 3 nt of a boundary does not
receive it. Listed in `docs/limitations.md`.

### D-019: Gene-level rank is the primary metric, with a pessimistic tie policy
**Chose:** rank the causal **gene**, taking the worst position in its tie block; report the
optimistic policy and the variant-level figure alongside.
**Alternatives:** variant-level rank, which was the first implementation.
**Reasoning:** Exomiser assigns a combined score per gene, and every variant in that gene inherits
it. A variant-level rank therefore counts a gene's *other* variants against the causal one. For
compound heterozygotes this is not a bias but an impossibility: the two causal alleles necessarily
tie with each other, so variant-level rank 1 is unreachable and comp-het top-1 was reported as
0.000 for all 189 cases. That was a defect in the metric, not a property of the tool - gene-level
top-1 for the same cases is 0.460. The extractor now also excludes a case's own partner allele from
the tie block.
**Source:** own diagnosis, after comp-het top-1 came out identically zero across every condition -
a number too clean to be real.
**Effect:** gene-level informed top-1 is 0.418 against 0.333 variant-level. The pessimistic policy
still costs real recall relative to an undeclared tie-break, and both are reported.

### D-023: Inherited-dominant cases re-run with the transmitting parent coded affected
**Chose:** keep the original run (both parents unaffected) and add a matched contrast in which the
carrier father is coded affected.
**Reasoning:** the original construction asserts a heterozygous parent who carries a dominant
pathogenic allele and is unaffected - that is a non-penetrant carrier, and Exomiser's inheritance
filter correctly rejects such a variant under AUTOSOMAL_DOMINANT. Only 4.0% of causal alleles
survived the cascade in that arm, so the near-zero recall was the filter working as designed, not
the ranker failing. Coding the father affected turns a construction error into a controlled
measurement of what the trio inheritance filter costs under incomplete penetrance.
**Source:** own diagnosis, after `ad_inherited` returned survived=0.040 against 0.989 for de novo.
**Effect:** reported as a paired contrast rather than as an inheritance-model recall figure.

### D-024: The ClinVar effect is reported conditionally, not as a headline average
**Chose:** report the ClinVar blinding delta split by whether the causal allele is actually on
Exomiser's whitelist, and crossed with the phenotype axis, rather than as a single average.
**Alternatives:** report the marginal delta (+0.007 top-1) and conclude the whitelist does not
matter - which is what the first pass of the analysis said.
**Reasoning:** the marginal number is close to zero and would have supported "ClinVar lookup does
not drive prioritization". Conditioning shows the opposite. Only 300 of 1,034 causal alleles are
whitelisted, and phenotype scoring masks the effect on those. Split out:

| n=300, whitelisted | ClinVar on | ClinVar blinded |
|---|---|---|
| phenotype supplied | 0.743 | 0.720 |
| variant-only | **0.907** | **0.507** |

| n=734, not whitelisted | ClinVar on | ClinVar blinded |
|---|---|---|
| phenotype supplied | 0.285 | 0.285 |
| variant-only | 0.338 | 0.330 |

The lookup effect is **+0.400 top-1** where it applies, and exactly zero where it does not - the
lower table is an internal control confirming that blinding touches only the whitelisted variants.
**Effect:** the marginal delta would have been a wrong conclusion drawn from a correct calculation.
The conditional split is the reported result.

### D-025: The phenotype effect is reported with its interaction, not as a main effect
**Chose:** report that supplying phenotypes *lowers* top-1 when ClinVar is informing the variant
score (-0.085, p=6e-8) and slightly raises it when ClinVar is blinded (+0.030, p=0.10).
**Reasoning:** stating "phenotype ranking makes results worse" unconditionally would be false. The
mechanism is visible: the whitelist pins a causal variant's score at maximum, ranking on that alone
puts the causal gene first 90.7% of the time, and blending in a phenotype score below 1.0 dilutes a
signal that was already saturated. Real case-report phenotypes are sparse - median 7 terms - so the
causal gene's phenotype score frequently is below 1.
**Effect:** the honest statement is about an interaction, and the 2x2 is figure 8.
