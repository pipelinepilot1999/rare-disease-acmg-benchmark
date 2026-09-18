# Limitations

Ordered by how much they should change your reading of the results.

---

## 1. Cases are synthetic at the VCF level

Causal alleles are injected into a healthy trio's VCF. Nothing here tests variant
calling: sensitivity, coverage, mapping quality and caller behaviour are inherited from the GIAB
benchmark, not measured. A pipeline that never calls the causal variant would score perfectly in
this benchmark. Read every Arm A number as conditional on the variant having been called.

The phenotypes and the causal alleles are real - they come from published case reports curated into
GA4GH Phenopackets - but the genomic background is one healthy individual's, and the pairing
between phenotype and background is artificial.

## 2. The ClinVar blind in Arm A is partial

Exomiser cannot start without its ClinVar database, and that database is an MVStore rather than an
H2 SQL database, so it cannot be emptied with the tooling shipped alongside it. The blinded
condition therefore disables the ClinVar **whitelist** - the channel that lets a variant take a
maximal score and bypass the filter cascade - but ClinVar remains available to the annotation-only
`EXOMISER_ACMG_*` columns.

Those columns do not feed the ranking score, so the blind is complete for the quantity being
measured. It is not a claim that no ClinVar-derived information exists anywhere in the run.

## 3. One genomic background, and it is Ashkenazi

Every case uses HG002/HG003/HG004. Two consequences:

- **Background diversity is zero.** The matched-knockout design controls for this within the
  experiment - each positive is compared against its own variant-removed twin - but between-case
  variation in background is not represented. A cohort with different ancestries and different rare
  variant burdens would rank differently.
- **The proband is Ashkenazi Jewish, and the frequency machinery deliberately excludes that group.**
  gnomAD's grpmax filtering AF omits Ashkenazi, Finnish, Amish and Middle Eastern populations
  because pathogenic alleles drift high in bottlenecked groups; Exomiser's default frequency
  sources independently omit `GNOMAD_E_ASJ` and `GNOMAD_G_ASJ`. The proband's own population is
  therefore invisible to the frequency filter that is supposed to remove his common variation.
  Reported as a sub-analysis rather than buried.

## 4. A knockout is negative only for the spiked variant

A knockout case is "negative" in the sense that the *spiked* variant is absent, not that the genome
is free of pathogenic alleles. We checked rather than assumed: within the 41.3 Mb testable region
the trio carries **5 ClinVar P/LP alleles**, of which **one is heterozygous in the proband and none
are homozygous**. That is small enough not to distort the negative-control comparison, but it is a
measured number rather than an assumption, and the alleles are listed in `data/trio_carriers.tsv`.

The restriction to high-confidence exome regions is doing most of the work here; a whole-genome or
less strictly filtered background would carry more.

## 5. X-linked recessive is untestable here

Only HG002 has a published GIAB chrXY benchmark; the parents do not, so no trio can be built on
chrX. 351 chrX variants in the phenopacket corpus (3.2%) are excluded, and the X-linked recessive
inheritance model - one of the four v1 set out to test - has no arm. chrX accounts for 1.7 Mb
(3.7%) of the exome region set.

## 6. 6% of the autosomal exome is not testable

The GIAB v4.2.1 benchmark excludes difficult regions - segmental duplications, low-complexity
sequence, regions where callers disagree. Intersecting the three samples' high-confidence BEDs with
the GENCODE CDS +/- 20 bp region set leaves 41.3 Mb of 43.9 Mb autosomal. Causal variants cannot be
placed in the missing 2.6 Mb, so Arm A is optimistic relative to a real exome by an amount this
design cannot measure. 783 otherwise-usable phenopacket cases were dropped for exactly this reason.

## 7. Arm B's ground truth is not a random sample of variants

ClinGen Evidence Repository interpretations are concentrated in the genes VCEPs have chosen to
curate: the five largest panels account for a large share of the corpus. Results are stratified by
VCEP and by gene for this reason, and aggregate figures should not be read as performance on
"variants in general."

More importantly, **VCEP calls apply gene-specific criteria specifications**, not the generic 2015
rules. A generic engine scored against gene-specific expert rules conflates "the rules are wrong"
with "the rules were not specialised for this gene." The stratification separates the two; the
aggregate number does not.

## 8. acmgx does not implement every criterion

PM1, PP2, BP1 and BP3 are declared unimplemented rather than approximated, because defensible
versions need curated hotspot/domain annotation and gene-level missense constraint that were not
assembled. acmgx therefore systematically under-calls relative to VCEPs on those criteria; the size
of the effect is in `results/criterion_metrics.tsv`.

BP7 uses a distance-to-exon-boundary heuristic rather than SpliceAI. Arm A does use real SpliceAI,
through Exomiser's bundled pathogenicity source.

## 9. The ceiling analysis assumes the combining rules are the VCEPs' own

The ceiling is computed by deleting un-automatable evidence codes from each VCEP's own criteria set
and re-running the combining rules. Our 2015 implementation reproduces the VCEPs' published tier
for 92.2% of interpretations when given all of their evidence, which is what makes this credible -
but it is not 100%. The residual 7.8% is mostly gene-specific combining rules that differ from the
2015 defaults. The ceiling should be read as accurate to a few percentage points, not exactly.

## 10. The inherited-dominant arm models a non-penetrant carrier

Cases in that arm were built with the transmitting parent heterozygous and coded unaffected, which
for a dominant condition asserts non-penetrance. Exomiser's `AUTOSOMAL_DOMINANT` filter discards
such variants by design, and 96% of those causal alleles never reached the ranker. The arm is
therefore reported as a **penetrance contrast** - the same cases re-run with the parent coded
affected - and not as an inheritance-model recall figure. The contrast is informative, but it was
not the experiment originally intended.

## 11. The whitelist strata are defined by an observed flag, and 222 cases have no observation

Arm A's ClinVar stratification uses Exomiser's own `WHITELIST_VARIANT` output column, which is only
emitted for variants that survived the filter cascade. Under `informed`, 812 of 1,034 causal alleles
survived; 300 of those carry the flag. The "not whitelisted" group of 734 is therefore 512 alleles
observed not to be whitelisted plus 222 whose status was never observed. Restricting the comparison
to the 512 observed cases leaves the conclusion unchanged - the ClinVar delta is 0.000 either way -
but the group label is looser than it looks.

Relatedly, 63.6% of causal alleles are ClinVar P/LP (D-005) while only 29.0% carried the whitelist
flag. The gap is unexplained here. The likely cause is that the Exomiser 2512 data bundle's ClinVar
snapshot predates the 2026-08-31 ClinVar VCF used for the P/LP count, and that the whitelist applies
a review-status floor, but neither was confirmed. **Anyone quoting these numbers should treat 29.0%
as the measured whitelist coverage and 63.6% as an upper bound on it.**

## 12. Single run, no replication

Exomiser and acmgx are deterministic given the same inputs and versions, so there is no run-to-run
variance to report. All intervals here are sampling intervals over cases, not over runs. A change
of ClinVar release, gnomAD version, or Exomiser data bundle would move these numbers, and none of
that uncertainty is captured.
