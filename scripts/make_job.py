#!/usr/bin/env python3
"""Emit an Exomiser job YAML for one case under one condition."""
import sys, json, os

def job(vcf, hpo, out_name, out_dir, phenotype=True):
    hpo_block = "[ " + ", ".join(f"'{t}'" for t in hpo) + " ]" if (phenotype and hpo) else "[ ]"
    steps = []
    if phenotype:
        steps.append("      hiPhivePrioritiser: { }")
    steps += [
        "      failedVariantFilter: { }",
        "      regulatoryFeatureFilter: { }",
        "      frequencyFilter: { maxFrequency: 2.0 }",
        # keepNonPathogenic: nothing is discarded, so every case yields a complete
        # ranking and the causal variant always has a defined rank.
        "      pathogenicityFilter: { keepNonPathogenic: true, target: NON_CODING }",
        "      inheritanceFilter: { }",
    ]
    if phenotype:
        steps.append("      omimPrioritiser: { }")
    return f"""---
sample:
  genomeAssembly: hg38
  vcf: {vcf}
  ped: {os.path.expanduser('~/rdacmg/repo/config/trio.ped')}
  proband: HG002
  hpoIds: {hpo_block}
analysis:
  inheritanceModes: {{
    AUTOSOMAL_DOMINANT: 0.1,
    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
    X_DOMINANT: 0.1,
    X_RECESSIVE_HOM_ALT: 0.1,
    X_RECESSIVE_COMP_HET: 2.0
  }}
  analysisMode: PASS_ONLY
  frequencySources: [
      GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_EAS, GNOMAD_E_NFE, GNOMAD_E_SAS,
      GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_EAS, GNOMAD_G_NFE, GNOMAD_G_SAS
  ]
  pathogenicitySources: [ REVEL, MVP, ALPHA_MISSENSE, SPLICE_AI ]
  steps: [
{chr(10).join(s + ',' for s in steps[:-1])}
{steps[-1]}
  ]
outputOptions:
  outputContributingVariantsOnly: false
  numGenes: 0
  outputDirectory: {out_dir}
  outputFileName: {out_name}
  outputFormats: [ TSV_VARIANT, TSV_GENE ]
"""

if __name__ == "__main__":
    vcf, hpo_csv, out_name, out_dir, pheno = sys.argv[1:6]
    hpo = [h for h in hpo_csv.split(",") if h]
    print(job(vcf, hpo, out_name, out_dir, pheno == "1"))
